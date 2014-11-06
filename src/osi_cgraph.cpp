#include <cstring>
#include <cmath>
#include <limits>
#include <cfloat>
#include <algorithm>
#include <climits>
#include <OsiSolverInterface.hpp>
#include <CoinPackedMatrix.hpp>
#include "osi_cgraph.h"
extern "C"
{
#include "memory.h"
#include "node_heap.h"
}

#define oo  (INT_MAX/2)
#define WORST_PRIORITY_ROW 10000   // to be considered

using namespace std;

#define EPS 1e-6

/* mininum size for a row to be considered a clique row */
#define MIN_CLIQUE_ROW 250

#define IS_NEGATIVE_DBL( v ) ( v < -1e-5 )

#define DBL_EQUAL( v1, v2 ) ( fabs(v1-v2) < EPS )

// conflict vector dimensions
#define CVEC_CAP   1800000
#define CVEC_FLUSH 1700000

#define MAX_NEGATIVE_COEFS 200

#define MAX_DIFFERENT_COEFS 100

const double LARGE_CONST = std::min( DBL_MAX/10.0, 1e20 );

double maxRowTime = 0.0;
int maxRowTimeIdx = 0;
int success = 1;
clock_t start, end;
double cpuTime;
int cliqueComp = 0, clique = 0, activePairwise = 0, inactivePairwise = 0, mixedPairwise = 0;

vector< pair< int, int > > cvec;/* conflict vector */
vector< int > neighs;/* used i fetch conflicts */
double *colLb;
double *colUb;
int nCols;
int nRows;

#ifdef DEBUG_CONF
    int newConflicts;
    vector< pair<int,int> > confS;
#endif /* DEBUG_CONF */

struct sort_sec_pair
{
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right)
    {
        if ( left.second != right.second )
            return left.second < right.second;

        return left.first < right.first;
    }
};

struct sort_columns
{
    bool operator()(const std::pair<int, double> &left, const std::pair<int,double> &right)
    {
        if ( fabs(left.second - right.second) > EPS )
            return ( (left.second - right.second) < EPS );

        return left.first < right.first;
    }
};

bool cliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs);

bool pairwiseAnalysisByGrouping(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs);

void pairwiseAnalysis(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs);

void processClique( const int n, const int *idx, CGraph *cgraph, const double *colLb, const double *colUb );

/* if size is large enough fech conflicts */
void fetchConflicts( const bool lastTime, CGraph *cgraph );

/* returns how much a variable can negatively contribute */
double mostNegativeContribution( const double coef, const char s, const double colLb, const double colUb );

/* returns how much a variable can positively contribute */
double unitaryContribution( const double coef, const char s, const double colLb, const double colUb );

CGraph *osi_build_cgraph( void *_lp )
{
    start = clock();

    OsiSolverInterface *lp = (OsiSolverInterface *)_lp;

    if (lp->getNumIntegers()<2)
        return 0;

    int cgraphSize = lp->getNumCols() * 2; //considering binary complement
    									  //using extra memory to facilitate indexing

    CGraph *cgraph = cgraph_create( cgraphSize );
    const char *ctype = lp->getColType();
    const CoinPackedMatrix *M = lp->getMatrixByRow();
    const double *rhs = lp->getRightHandSide();
    colLb = (double*)lp->getColLower();
    colUb = (double*)lp->getColUpper();
    const char *sense = lp->getRowSense();
    nCols = lp->getNumCols();
    nRows = lp->getNumRows();
    int idxRow;
    cvec.reserve( CVEC_CAP );
    neighs.reserve( 8192 );

    for(idxRow = 0; idxRow < nRows; idxRow++)
    {
        clock_t rowStart = clock();
#define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )
        const CoinShallowPackedVector &row = M->getVector(idxRow);
        const int nElements = row.getNumElements();
        const int *idx = row.getIndices();
        const double *coefs = row.getElements();
        vector< pair<int, double> > columns(nElements);
        int nBools = 0; // number of binary variables
        double sumNegCoefs = 0.0; //sum of all negative coefficients
        
        if ( (nElements<2) || (fabs(rhs[idxRow])>=LARGE_CONST) )
            continue;

        end = clock();
        cpuTime = ((double) (end - start)) / CLOCKS_PER_SEC;
        if(cpuTime > 60.0)
        {
            success = 0;
            break;
        }

        if ( sense[idxRow] == 'R' )  // lets not consider ranged constraints by now
        {
            printf("TODO: CHECK FOR RANGED CONSTRAINT (%s) rhs is %g\n", lp->getRowName(idxRow).c_str(), rhs[idxRow] );
            continue;
        }

#ifdef DEBUG_CONF
        newConflicts = 0;
        confS.clear();
#endif /* DEBUG_CONF */

		double mult = (sense[idxRow] == 'G') ? mult = -1.0 : mult = 1.0;
        for(int i = 0; i < nElements; i++)
        {
        	const int cidx = idx[i];

        	columns[i].first = idx[i];
        	columns[i].second = coefs[i] * mult;

        	if(ctype[cidx] == 1)
            	nBools++;

            if (FIXED_IN_ZERO(cidx))
                continue;

			if(columns[i].second <= -EPS)
            	sumNegCoefs += columns[i].second;

        	/* inserting trivial conflicts: variable-complement */
            if(ctype[cidx] != 1) //consider only binary variables
            	continue;
            cvec.push_back( pair<int, int>(cidx, cidx + nCols) );
#ifdef DEBUG_CONF
            newConflicts++;
            confS.push_back( pair<int,int>(cidx, cidx + nCols) );
#endif /* DEBUG_CONF */
        }

        /* considering just constraints which have only binary variables */
        if(nBools < nElements)
        	continue;

        sort(columns.begin(), columns.end(), sort_columns());
        bool foundClique = cliqueDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);

        if(!foundClique)/* normal case, checking pairwise conflicts
             * (these are computed increasing RHS in thisRhs) */
        {
        	bool groupingWorked = pairwiseAnalysisByGrouping(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);

	    	if(!groupingWorked)
	    		pairwiseAnalysis(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);
        }

        /*equality constraints are converted into two inequality constraints (<=).
        the first one is analyzed above and the second (multiplying the constraint by -1) is analyzed below.
        Example: x + y + z = 2 ==>  (x + y + z <= 2) and (- x - y - z <= -2)*/
        if(sense[idxRow] == 'E')
        {
            vector<pair<int, double> > newColumns(nElements);
            sumNegCoefs = 0.0;
            for(int i = 0; i < nElements; i++)
            {
                newColumns[i].first = columns[nElements-i-1].first;
                newColumns[i].second = -1.0 * columns[nElements-i-1].second;
                if(newColumns[i].second <= -EPS)
                    sumNegCoefs += newColumns[i].second;
            }

            foundClique = cliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);

            if(!foundClique)/* normal case, checking pairwise conflicts
                 * (these are computed increasing RHS in thisRhs) */
            {
                bool groupingWorked = pairwiseAnalysisByGrouping(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);

                if(!groupingWorked)
                    pairwiseAnalysis(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);
            }            
        }

#ifdef DEBUG_CONF
        if (newConflicts)
        {
            printf("%s,%d,%c,%g : ", lp->getRowName(idxRow).c_str(), nNeg, sense[idxRow], rhs[idxRow]);
            for ( int cc=0 ; (cc<confS.size()) ; ++cc )
            {
                const int var1 = confS[cc].first, var2 = confS[cc].second;
                printf("(%s, %s) ", var1 < nCols ? lp->getColName(var1).c_str() : ("¬" + lp->getColName(var1 - nCols)).c_str(),
                                    var2 < nCols ? lp->getColName(var2).c_str() : ("¬" + lp->getColName(var2 - nCols)).c_str() );
            }
            printf("\n");
        }
#endif
#undef FIXED_IN_ZERO
        clock_t rowEnd = clock();
        double rowTime = ((double) (rowEnd - rowStart)) / CLOCKS_PER_SEC;

        if(rowTime > maxRowTime)
        {
            maxRowTime = rowTime;
            maxRowTimeIdx = idxRow;
        }
    }

    fetchConflicts(true, cgraph);
    cgraph_update_min_max_degree( cgraph );

    printf("Clique: %d \t CliqueComp: %d \t Active: %d \t Inactive: %d \t Mixed: %d\n",
             clique, cliqueComp, activePairwise, inactivePairwise, mixedPairwise);

    return cgraph;
}

void processClique( const int n, const int *idx, CGraph *cgraph, const double *colLb, const double *colUb )
{
    if ( n >= MIN_CLIQUE_ROW )
    {
        vector< int > vidx;
        vidx.resize(n);

        memcpy( &(vidx[0]), idx, sizeof(int)*n );

        cgraph_add_clique( cgraph, &(vidx[0]), n );
    }
    else
    {
        int i1, i2, nm1 = n-1;
        int nCols = cgraph_size(cgraph) / 2; //divide by 2 to eliminate the complement of the variables

        for ( i1=0 ; (i1<nm1) ; ++i1 )
        {
        	if ( fabs( colUb[idx[i1]%nCols] <= EPS ) && fabs( colLb[idx[i1]%nCols] <= EPS ) )
                continue;
            for ( i2=i1+1 ; (i2<n) ; ++i2 )
            {
            	if ( fabs( colUb[idx[i2]%nCols] <= EPS ) && fabs( colLb[idx[i2]%nCols] <= EPS ) )
                    continue;
                cvec.push_back( pair<int,int>(idx[i1],idx[i2]) );
            }
        }
    }
}

void fetchConflicts( const bool lastTime, CGraph *cgraph )
{
    if ( (!lastTime) && (cvec.size()<CVEC_FLUSH) )
        return;

    if ( cvec.size() == 0 )
        return;
    pair< int, int > last;
    int currNode;
    vector< pair<int,int> >::const_iterator vIt;

    cvec.push_back( pair<int,int>( INT_MAX, INT_MAX ) );

    neighs.clear();
    sort( cvec.begin(), cvec.end() );
    currNode = cvec.begin()->first;
    last = pair< int,int >( -1, -1 );

    for ( vIt=cvec.begin() ; (vIt!=cvec.end()) ; last=*vIt,++vIt )
    {
        if (last == *vIt)  // skipping repetitions
            continue;

        if ( vIt->first!=currNode )
        {
            if (neighs.size())
            {
                cgraph_add_node_conflicts_no_sim( cgraph, currNode, &(neighs[0]), neighs.size() );
                neighs.clear();
            }
            currNode = vIt->first;
        }
        neighs.push_back( vIt->second );
    }

    neighs.clear();
    sort( cvec.begin(), cvec.end(), sort_sec_pair() );
    currNode = cvec.begin()->second;
    last = pair< int,int >( -1, -1 );
    for ( vIt=cvec.begin() ; (vIt!=cvec.end()) ; last=*vIt,++vIt )
    {
        if (last == *vIt)  // skipping repetitions
            continue;

        if ( vIt->second!=currNode )
        {
            if (neighs.size())
            {
                cgraph_add_node_conflicts_no_sim( cgraph, currNode, &(neighs[0]), neighs.size() );
                neighs.clear();
            }
            currNode = vIt->second;
        }
        neighs.push_back( vIt->first );
    }

    cvec.clear();
}


double mostNegativeContribution( const double coef, const char s, const double colLb, const double colUb )
{
    double mult = 1.0;

    if (s=='G')
        mult = -1.0;

    double c = coef*mult;

    return min( c*colLb, c*colUb );
}

double unitaryContribution( const double coef, const char s, const double colLb, const double colUb )
{
    if ( colUb <= EPS )
        return 0.0;

    if (s=='G')
        return coef*-1.0;
    else
        return coef;

    // some dumb compilers print warnings if
    // function does not have a return outside ifs
    return coef;
}

void pairwiseAnalysis(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs)
{

	int nElements = (int)columns.size();
    for(int j1 = 0; j1 < nElements; j1++)
    {
        const int cidx1 = columns[j1].first;
        const double coef1 = columns[j1].second;

        #define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )

        if (FIXED_IN_ZERO(cidx1))
            continue;

        for(int j2 = j1+1; j2 < nElements; j2++)
        {
            const int cidx2 = columns[j2].first;
            const double coef2 = columns[j2].second;

            if(FIXED_IN_ZERO(cidx2))
                continue;

            const double negDiscount = sumNegCoefs - min(0.0, coef1) - min(0.0, coef2);

            if(coef1 + coef2 + negDiscount > rhs + 0.001)
            {
                cvec.push_back( pair<int,int>(cidx1, cidx2) );
                activePairwise++;
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1, cidx2) );
#endif
            }

            if(coef1 + negDiscount > rhs + 0.001) /* cidx1 = 1 and cidx2 = 0 */
            {
                cvec.push_back( pair<int,int>(cidx1, cidx2+nCols) );
                mixedPairwise++;
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1, cidx2+nCols) );
#endif
            }

            if(coef2 + negDiscount > rhs + 0.001) /* cidx1 = 0 and cidx2 = 1 */
            {
                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2) );
                mixedPairwise++;
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1+nCols, cidx2) );
#endif
            }

            if(negDiscount > rhs + 0.001) /* cidx1 = 0 and cidx2 = 0 */
            {
                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2+nCols) );
                inactivePairwise++;
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1+nCols, cidx2+nCols) );
#endif
            }
        }
        fetchConflicts( false, cgraph );
    }
	#undef FIXED_IN_ZERO
}

bool pairwiseAnalysisByGrouping(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs)
{
	int nElements = (int)columns.size();
	vector<vector<int> > coefsGroup; //coefficients with same value will be grouped
	vector<double> diffCoefs; //different coefficients in the current row
	double prevCoef = -LARGE_CONST;
	for(int i = 0; i < nElements; i++)
	{
    	if( ((int)(diffCoefs.size())) >= MAX_DIFFERENT_COEFS)
    		return false;

    	if(columns[i].second > prevCoef)
    	{
    		vector<int> tmp(1, columns[i].first);
    		coefsGroup.push_back(tmp);
    		diffCoefs.push_back(columns[i].second);
    		prevCoef = columns[i].second;
    	}
    	else coefsGroup.rbegin()->push_back(columns[i].first);
	}

	#define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )

	for(int i = 0; i < (int)diffCoefs.size(); i++)
	{
		const double coef1 = diffCoefs[i];
		double negDiscount = sumNegCoefs - min(0.0, coef1) - min(0.0, coef1);

		if(coef1 + coef1 + negDiscount > rhs + 0.001)
		{
			for(int k = 0; k < (int)coefsGroup[i].size(); k++)
        	{
        		if (FIXED_IN_ZERO(coefsGroup[i][k]))
        			continue;
        		for(int l = k+1; l < (int)coefsGroup[i].size(); l++)
        		{
        			if (FIXED_IN_ZERO(coefsGroup[i][l]))
        				continue;
            		cvec.push_back(pair<int,int>(coefsGroup[i][k], coefsGroup[i][l]));
                    activePairwise++;
#ifdef DEBUG_CONF
            		newConflicts++;
            		confS.push_back(pair<int,int>(coefsGroup[i][k], coefsGroup[i][l]));
#endif
            	}
            }	
		}

		if(coef1 + negDiscount > rhs + 0.001)
		{
			for(int k = 0; k < (int)coefsGroup[i].size(); k++)
        	{
        		if (FIXED_IN_ZERO(coefsGroup[i][k]))
        			continue;
        		for(int l = k+1; l < (int)coefsGroup[i].size(); l++)
        		{
        			if (FIXED_IN_ZERO(coefsGroup[i][l]))
        				continue;
            		cvec.push_back( pair<int,int>(coefsGroup[i][k], coefsGroup[i][l]+nCols));
            		cvec.push_back( pair<int,int>(coefsGroup[i][k]+nCols, coefsGroup[i][l]));
                    mixedPairwise+=2;
#ifdef DEBUG_CONF
            		newConflicts+=2;
            		confS.push_back( pair<int,int>(coefsGroup[i][k], coefsGroup[i][l]+nCols) );
            		confS.push_back( pair<int,int>(coefsGroup[i][k]+nCols, coefsGroup[i][l]) );
#endif
            	}
            }
		}

		if(negDiscount > rhs + 0.001)
		{
			for(int k = 0; k < (int)coefsGroup[i].size(); k++)
        	{
        		if (FIXED_IN_ZERO(coefsGroup[i][k]))
        			continue;
        		for(int l = k+1; l < (int)coefsGroup[i].size(); l++)
        		{
        			if (FIXED_IN_ZERO(coefsGroup[i][l]))
        				continue;
            		cvec.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[i][l]+nCols));
                    inactivePairwise++;
#ifdef DEBUG_CONF
            		newConflicts++;
            		confS.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[i][l]+nCols) );
#endif
            	}
            }
		}

		for(int j = i + 1; j < (int)diffCoefs.size(); j++)
		{
			const double coef2 = diffCoefs[j];
			negDiscount = sumNegCoefs - min(0.0, coef1) - min(0.0, coef2);

            if(coef1 + coef2 + negDiscount > rhs + 0.001)
            {
            	for(int k = 0; k < (int)coefsGroup[i].size(); k++)
            	{
            		if (FIXED_IN_ZERO(coefsGroup[i][k]))
            			continue;
            		for(int l = 0; l < (int)coefsGroup[j].size(); l++)
            		{
            			if (FIXED_IN_ZERO(coefsGroup[j][l]))
            				continue;
                		cvec.push_back( pair<int,int>(coefsGroup[i][k], coefsGroup[j][l]));
                        activePairwise++;
#ifdef DEBUG_CONF
                		newConflicts++;
                		confS.push_back( pair<int,int>(coefsGroup[i][k], coefsGroup[j][l]) );
#endif
                	}
                }
            }
            if(coef1 + negDiscount > rhs + 0.001)
            {
            	for(int k = 0; k < (int)coefsGroup[i].size(); k++)
            	{
            		if (FIXED_IN_ZERO(coefsGroup[i][k]))
            			continue;
            		for(int l = 0; l < (int)coefsGroup[j].size(); l++)
            		{
            			if (FIXED_IN_ZERO(coefsGroup[j][l]))
            				continue;
                		cvec.push_back( pair<int,int>(coefsGroup[i][k], coefsGroup[j][l]+nCols) );
                        mixedPairwise++;
#ifdef DEBUG_CONF
                		newConflicts++;
                		confS.push_back( pair<int,int>(coefsGroup[i][k], coefsGroup[j][l]+nCols) );
#endif
                	}
                }
            }
            if(coef2 + negDiscount > rhs + 0.001)
            {
            	for(int k = 0; k < (int)coefsGroup[i].size(); k++)
            	{
            		if (FIXED_IN_ZERO(coefsGroup[i][k]))
            			continue;
            		for(int l = 0; l < (int)coefsGroup[j].size(); l++)
            		{
            			if (FIXED_IN_ZERO(coefsGroup[j][l]))
            				continue;
                		cvec.push_back( pair<int,int>(coefsGroup[i][k]+nCols, coefsGroup[j][l]) );
                        mixedPairwise++;
#ifdef DEBUG_CONF
                		newConflicts++;
                		confS.push_back( pair<int,int>(coefsGroup[i][k]+nCols, coefsGroup[j][l]) );
#endif
                	}
                }
            }
            if(negDiscount > rhs + 0.001)
            {
            	for(int k = 0; k < (int)coefsGroup[i].size(); k++)
            	{
            		if (FIXED_IN_ZERO(coefsGroup[i][k]))
            			continue;
            		for(int l = 0; l < (int)coefsGroup[j].size(); l++)
            		{
            			if (FIXED_IN_ZERO(coefsGroup[j][l]))
            				continue;
                		cvec.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[j][l]+nCols) );
                        inactivePairwise++;
#ifdef DEBUG_CONF
                		newConflicts++;
                		confS.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[j][l]+nCols) );
#endif
                	}
                }
            }
		}
		fetchConflicts(false, cgraph);
	}
	#undef FIXED_IN_ZERO
	return true;
}

bool cliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs)
{
	int nElements = (int)columns.size(), cliqueStart = -1, cliqueCompStart = -1;
    double L00, L11; //L00 = lower bound for LHS when the two variaveis with smallest coefficients are deactivated.
    				 //L11 = lower bound for LHS when the two variaveis with highest coefficients are activated.

    L00 = sumNegCoefs - min(0.0, columns[1].second) - min(0.0, columns[0].second);
    L11 = sumNegCoefs - min(0.0, columns[nElements-2].second) - min(0.0, columns[nElements-1].second)
    	  + columns[nElements-2].second + columns[nElements-1].second;

    if(L00 <= rhs && L11 <= rhs)
    	return false;

    if(L00 > rhs)
    {
        for(int i = nElements - 1; i > 0; i--)
        {
            double D = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[i-1].second);

            if(D > rhs)
            {
                cliqueCompStart = i;
                break;
            }
        }
        assert(cliqueCompStart > 0 && cliqueCompStart < nElements);
        int n = cliqueCompStart + 1, idxs[n];
        for(int i = 0; i < n; i++)
            idxs[i] = columns[i].first + nCols; //binary complement
        processClique( n, (const int *)idxs, cgraph, colLb, colUb );
        cliqueComp++;
    }

    if(L11 > rhs)
    {
        for(int i = 0; i < nElements - 1; i++)
        {
            double D = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[i+1].second);
            double L = D + columns[i].second + columns[i+1].second;

            if(L > rhs)
            {
                cliqueStart = i;
                break;
            }
        }
        assert(cliqueStart >= 0 && cliqueStart < nElements - 1);
        int n = nElements - cliqueStart, idxs[n];
        for(int i = cliqueStart, j = 0; i < nElements; i++)
            idxs[j++] = columns[i].first;
        processClique( n, (const int *)idxs, cgraph, colLb, colUb );
        clique++;
    }
    return true;
}