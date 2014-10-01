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
unsigned long int completeClique = 0L, incompleteClique = 0L, completeCliqueComplement = 0L, incompleteCliqueComplement = 0L,
					activePairwise = 0L, inactivePairwise = 0L, mixedPairwise = 0L;

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

void processClique( const int n, const int *idx, CGraph *cgraph, vector< pair<int, int> > &cvec, const double *colLb, const double *colUb );

/* if size is large enough fech conflicts */
void fetchConflicts( vector< pair<int,int> > &cvec, const bool lastTime, CGraph *cgraph, vector<int> &neighs );

/* returns how much a variable can negatively contribute */
double mostNegativeContribution( const double coef, const char s, const double colLb, const double colUb );

/* returns how much a variable can positively contribute */
double unitaryContribution( const double coef, const char s, const double colLb, const double colUb );

CGraph *osi_build_cgraph( void *_lp )
{
    clock_t start = clock(), end;
    double cpuTime;

    OsiSolverInterface *lp = (OsiSolverInterface *)_lp;

    if (lp->getNumIntegers()<2)
        return 0;

    int cgraphSize = lp->getNumCols() * 2; //considering binary complement
    									  //using extra memory to facilitate indexing

    CGraph *cgraph = cgraph_create( cgraphSize );
    const char *ctype = lp->getColType();
    const CoinPackedMatrix *M = lp->getMatrixByRow();
    const double *rhs = lp->getRightHandSide();
    const double *colLb = lp->getColLower();
    const double *colUb = lp->getColUpper();
    const char *sense = lp->getRowSense();
    const int nCols = lp->getNumCols();
    const int nRows = lp->getNumRows();
    int idxRow;

    /* conflict vector */
    vector< pair< int, int > > cvec;
    cvec.reserve( CVEC_CAP );

    /* used i fetch conflicts */
    vector< int > neighs;
    neighs.reserve( 8192 );

    for(idxRow = 0; idxRow < nRows; idxRow++)
    {
        clock_t rowStart = clock();
#define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )
        const CoinShallowPackedVector &row = M->getVector(idxRow);
        const int nElements = row.getNumElements();
        const int *idx = row.getIndices();
        const double *coefs = row.getElements();
        pair<int, double> columns[nElements];
        double minCoef = numeric_limits<double>::max();
        double maxCoef = -1.0 * numeric_limits<double>::max();
        int nBools = 0; // number of binary variables
        int posCoefs = 0; //number of positive coefficients
        int negCoefs = 0; //number of negative coefficients
        double sumPosCoefs = 0.0; //sum of all positive coefficients
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
            //printf("TODO: CHECK FOR RANGED CONSTRAINT (%s) rhs is %g\n", lp->getRowName(idxRow).c_str(), rhs[idxRow] );
            continue;
        }

#ifdef DEBUG_CONF
        int newConflicts = 0;
        vector< pair<int,int> > confS;
        confS.clear();
#endif /* DEBUG_CONF */

		double mult = 1.0;
      	if ( sense[idxRow] == 'G' )
        	mult = -1.0;
      	double thisRhs = rhs[idxRow] * mult;
        for(int i = 0; i < nElements; i++)
        {
        	const int cidx = idx[i];

        	columns[i].first = idx[i];
        	columns[i].second = coefs[i];

        	if(ctype[cidx] == 1)
            	nBools++;
            else break;

            if (FIXED_IN_ZERO(cidx))
                continue;

            const double neg = mostNegativeContribution( coefs[i], sense[idxRow], colLb[cidx], colUb[cidx] );

            if(coefs[i] >= -EPS)
            {
            	posCoefs++;
            	sumPosCoefs += coefs[i];
            }
            else if(coefs[i] <= -EPS)
            {
            	negCoefs++;
            	sumNegCoefs += coefs[i];
            }

            minCoef = min( minCoef, coefs[i] );
            maxCoef = max( maxCoef, coefs[i] );

            thisRhs -= neg;

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

        /* special case (to ignore), when all variables are fixed to zero */
        if ( (fabs(thisRhs)<=EPS) && DBL_EQUAL(maxCoef,minCoef)
                && DBL_EQUAL(maxCoef,1.0) )
            continue;

        sort(columns, columns + nElements, sort_columns());
        int C1 = -1, CC = -1; //k1: index where the formula Sk + Sk+1 > bi is satisfied (named C1 in document)
        					  //k2: index where the formula SA - Sk - Sk+1 < bi is satisfied (named CC in document)
        double Sminus, Splus;
        Sminus = sumNegCoefs - min(0.0, columns[nElements-2].second) - min(0.0, columns[nElements-1].second);
        Splus = sumPosCoefs - max(0.0, columns[nElements-2].second) - max(0.0, columns[nElements-1].second);

        if( (rhs[idxRow] == 'E' || rhs[idxRow] == 'L') && (Sminus + columns[nElements-2].second + columns[nElements-1].second <= rhs[idxRow]) )
        	continue; //there are no conflicts in this row

        if( (rhs[idxRow] == 'E' || rhs[idxRow] == 'G') && (Splus  >= rhs[idxRow]) )
        	continue; //there are no conflicts in this row

        for(int i = 0; i < nElements - 1; i++)
        {
        	double sk1, sk2;
        	sk1 = columns[i].second;
        	sk2 = columns[i+1].second;
        	Sminus = sumNegCoefs - min(0.0, sk1) - min(0.0, sk2);
        	Splus = sumPosCoefs - max(0.0, sk1) - max(0.0, sk2);
        	
        	if(C1 != -1 && CC != -1)
        		break;

        	if( (C1 == -1) && (sense[idxRow] ==  'E' || sense[idxRow] == 'L') && (Sminus + sk1 + sk2 > rhs[idxRow]) )
        		C1 = i;
        	if( (CC == -1) && (sense[idxRow] ==  'E' || sense[idxRow] == 'G') && (Splus < rhs[idxRow]) )
        		CC = i;
        }

        if(C1 != -1 || CC != -1)
        {
	        //special case: ready to use cliques
	        if(C1 != -1)
	        {
	        	int n = nElements - C1;
	        	int idxs[n];
	        	for(int i = C1, j = 0; i < nElements; i++)
	        		idxs[j++] = columns[i].first;
	        	processClique( n, (const int *)idxs, cgraph, cvec, colLb, colUb );	
	        	if(C1 == 0)
	        		completeClique++;
	        	else incompleteClique++;
	        }

	        //another special case: ready to use cliques (complement of binary variables)
	        if(CC != -1)
	        {
	        	int n = nElements - CC;
	        	int idxs[n];
	        	for(int i = CC, j = 0; i < nElements; i++)
	        		idxs[j++] = columns[i].first + nCols; //binary complement
	        	processClique( n, (const int *)idxs, cgraph, cvec, colLb, colUb );
	        	if(CC == 0)
	        		completeCliqueComplement++;
	        	else incompleteCliqueComplement++;
	        }
	    }

        else /* normal case, checking pairwise conflicts
             * (these are computed increasing RHS in thisRhs) */
        {
        	vector<vector<int> > coefsGroup; //coefficients with same value will be grouped
        	vector<double> diffCoefs; //different coefficients in the current row
        	double prevCoef = -LARGE_CONST;
        	for(int i = 0; i < nElements; i++)
        	{
	        	if( ((int)(diffCoefs.size())) >= MAX_DIFFERENT_COEFS)
	        		break;

	        	if(columns[i].second > prevCoef)
	        	{
	        		vector<int> tmp(1, columns[i].first);
	        		coefsGroup.push_back(tmp);
	        		diffCoefs.push_back(columns[i].second);
	        		prevCoef = columns[i].second;
	        	}
	        	else coefsGroup.rbegin()->push_back(columns[i].first);
	    	}

	    	if( ((int)(diffCoefs.size())) < MAX_DIFFERENT_COEFS )
	    	{
	    		for(int i = 0; i < (int)diffCoefs.size(); i++)
	    		{
	    			const double pos1 = unitaryContribution( diffCoefs[i], sense[idxRow], 0.0, 1.0 );
	    			const double neg1 = mostNegativeContribution( diffCoefs[i], sense[idxRow], 0.0, 1.0 );
                    double newThisRhs = thisRhs + neg1 + neg1;

	    			/* variable cannot be activated anyway, this is the case of integer variables
	                 * with strictly positive lower bounds, e.g. miplib 2010 30n20b8, variable s4 */
	    			if ( pos1-EPS >= thisRhs )
	    				continue;

	    			if (pos1+pos1>newThisRhs+0.001)
	    			{
	    				for(int k = 0; k < (int)coefsGroup[i].size(); k++)
                    	{
                    		if (FIXED_IN_ZERO(coefsGroup[i][k]))
                    			continue;
                    		for(int l = k+1; l < (int)coefsGroup[i].size(); l++)
                    		{
                    			if (FIXED_IN_ZERO(coefsGroup[i][l]))
                    				continue;
                    			activePairwise++;
                        		cvec.push_back( pair<int,int>(coefsGroup[i][k],coefsGroup[i][l]));
#ifdef DEBUG_CONF
                        		newConflicts++;
                        		confS.push_back( pair<int,int>(coefsGroup[i][k],coefsGroup[i][l]) );
#endif
                        	}
                        }	
	    			}

	    			if (pos1>newThisRhs+0.001)
	    			{
	    				for(int k = 0; k < (int)coefsGroup[i].size(); k++)
                    	{
                    		if (FIXED_IN_ZERO(coefsGroup[i][k]))
                    			continue;
                    		for(int l = k+1; l < (int)coefsGroup[i].size(); l++)
                    		{
                    			if (FIXED_IN_ZERO(coefsGroup[i][l]))
                    				continue;
                    			mixedPairwise+=2;
                        		cvec.push_back( pair<int,int>(coefsGroup[i][k],coefsGroup[i][l]+nCols));
                        		cvec.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[i][l]));
#ifdef DEBUG_CONF
                        		newConflicts+=2;
                        		confS.push_back( pair<int,int>(coefsGroup[i][k],coefsGroup[i][l]+nCols) );
                        		confS.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[i][l]) );
#endif
                        	}
                        }
	    			}
	    			
	    			if (newThisRhs<=-EPS)
	    			{
	    				for(int k = 0; k < (int)coefsGroup[i].size(); k++)
                    	{
                    		if (FIXED_IN_ZERO(coefsGroup[i][k]))
                    			continue;
                    		for(int l = k+1; l < (int)coefsGroup[i].size(); l++)
                    		{
                    			if (FIXED_IN_ZERO(coefsGroup[i][l]))
                    				continue;
                    			inactivePairwise++;
                        		cvec.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[i][l]+nCols));
#ifdef DEBUG_CONF
                        		newConflicts++;
                        		confS.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[i][l]+nCols) );
#endif
                        	}
                        }
	    			}

	    			for(int j = i + 1; j < (int)diffCoefs.size(); j++)
	    			{
	    				const double pos2 = unitaryContribution( diffCoefs[j], sense[idxRow], 0.0, 1.0 );
	    				
	    				if ( pos2-EPS >= thisRhs )
	    					continue;

	                    const double neg2 = mostNegativeContribution( diffCoefs[j], sense[idxRow], 0.0, 1.0 );
	                    newThisRhs = thisRhs + neg1 + neg2;

	                    if (pos1+pos2>newThisRhs+0.001)
	                    {
	                    	for(int k = 0; k < (int)coefsGroup[i].size(); k++)
	                    	{
	                    		if (FIXED_IN_ZERO(coefsGroup[i][k]))
	                    			continue;
	                    		for(int l = 0; l < (int)coefsGroup[j].size(); l++)
	                    		{
	                    			if (FIXED_IN_ZERO(coefsGroup[j][l]))
	                    				continue;
	                    			activePairwise++;
	                        		cvec.push_back( pair<int,int>(coefsGroup[i][k],coefsGroup[j][l]));
	#ifdef DEBUG_CONF
	                        		newConflicts++;
	                        		confS.push_back( pair<int,int>(coefsGroup[i][k],coefsGroup[j][l]) );
	#endif
	                        	}
	                        }
	                    }
	                    if (pos1>newThisRhs+0.001)
	                    {
	                    	for(int k = 0; k < (int)coefsGroup[i].size(); k++)
	                    	{
	                    		if (FIXED_IN_ZERO(coefsGroup[i][k]))
	                    			continue;
	                    		for(int l = 0; l < (int)coefsGroup[j].size(); l++)
	                    		{
	                    			if (FIXED_IN_ZERO(coefsGroup[j][l]))
	                    				continue;
	                    			mixedPairwise++;
	                        		cvec.push_back( pair<int,int>(coefsGroup[i][k],coefsGroup[j][l]+nCols) );
	#ifdef DEBUG_CONF
	                        		newConflicts++;
	                        		confS.push_back( pair<int,int>(coefsGroup[i][k],coefsGroup[j][l]+nCols) );
	#endif
	                        	}
	                        }
	                    }
	                    if (pos2>newThisRhs+0.001)
	                    {
	                    	for(int k = 0; k < (int)coefsGroup[i].size(); k++)
	                    	{
	                    		if (FIXED_IN_ZERO(coefsGroup[i][k]))
	                    			continue;
	                    		for(int l = 0; l < (int)coefsGroup[j].size(); l++)
	                    		{
	                    			if (FIXED_IN_ZERO(coefsGroup[j][l]))
	                    				continue;
	                    			mixedPairwise++;
	                        		cvec.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[j][l]) );
	#ifdef DEBUG_CONF
	                        		newConflicts++;
	                        		confS.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[j][l]) );
	#endif
	                        	}
	                        }
	                    }
	                    if (newThisRhs<=-EPS)
	                    {
	                    	for(int k = 0; k < (int)coefsGroup[i].size(); k++)
	                    	{
	                    		if (FIXED_IN_ZERO(coefsGroup[i][k]))
	                    			continue;
	                    		for(int l = 0; l < (int)coefsGroup[j].size(); l++)
	                    		{
	                    			if (FIXED_IN_ZERO(coefsGroup[j][l]))
	                    				continue;
	                    			inactivePairwise++;
	                        		cvec.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[j][l]+nCols) );
	#ifdef DEBUG_CONF
	                        		newConflicts++;
	                        		confS.push_back( pair<int,int>(coefsGroup[i][k]+nCols,coefsGroup[j][l]+nCols) );
	#endif
	                        	}
	                        }
	                    }
	    			}
	    			fetchConflicts( cvec, false, cgraph, neighs );
	    		}
	    	}
	    	else
	    	{
	            int j1, j2, cm1 = nElements-1;
	            for ( j1=0 ; (j1<cm1) ; ++j1 )
	            {
	                const int cidx1 = idx[j1];

	                if (FIXED_IN_ZERO(cidx1))
	                    continue;

	                const double pos1 = unitaryContribution( coefs[j1], sense[idxRow], colLb[cidx1], colUb[cidx1] );

	                end = clock();
	                cpuTime = ((double) (end - start)) / CLOCKS_PER_SEC;
	                if(cpuTime > 60.0)
	                {
	                    success = 0;
	                    break;
	                }

	                /* variable cannot be activated anyway, this is the case of integer variables
	                 * with strictly positive lower bounds, e.g. miplib 2010 30n20b8, variable s4 */
	                if ( pos1-EPS >= thisRhs )
	                    continue;

	                for ( j2=j1+1 ; (j2<nElements) ; ++j2 )
	                {
	                    const int cidx2 = idx[j2];

	                    if (FIXED_IN_ZERO(cidx2))
	                        continue;

	                    const double pos2 = unitaryContribution( coefs[j2], sense[idxRow], colLb[cidx2], colUb[cidx2] );

	                    if ( pos2-EPS >= thisRhs )
	                        continue;

	                    end = clock();
	                    cpuTime = ((double) (end - start)) / CLOCKS_PER_SEC;
	                    if(cpuTime > 60.0)
	                    {
	                        success = 0;
	                        break;
	                    }

	                    /* recalculating thisRhs, excluding fixed variables */
	                    const double neg1 = mostNegativeContribution( coefs[j1], sense[idxRow], colLb[cidx1], colUb[cidx1] );
	                    const double neg2 = mostNegativeContribution( coefs[j2], sense[idxRow], colLb[cidx2], colUb[cidx2] );
	                    const double newThisRhs = thisRhs + neg1 + neg2;

	                    if (pos1+pos2>newThisRhs+0.001)
	                    {
	                    	activePairwise++;
	                        cvec.push_back( pair<int,int>(cidx1,cidx2) );
	#ifdef DEBUG_CONF
	                        newConflicts++;
	                        confS.push_back( pair<int,int>(cidx1,cidx2) );
	#endif
	                    }

	                    if (pos1>newThisRhs+0.001) /* cidx1 = 1 and cidx2 = 0 */
	                    {
	                    	mixedPairwise++;
	                        cvec.push_back( pair<int,int>(cidx1,cidx2+nCols) );
	#ifdef DEBUG_CONF
	                        newConflicts++;
	                        confS.push_back( pair<int,int>(cidx1,cidx2+nCols) );
	#endif
	                    }

	                    if (pos2>newThisRhs+0.001) /* cidx1 = 0 and cidx2 = 1 */
	                    {
	                    	mixedPairwise++;
	                        cvec.push_back( pair<int,int>(cidx1+nCols,cidx2) );
	#ifdef DEBUG_CONF
	                        newConflicts++;
	                        confS.push_back( pair<int,int>(cidx1+nCols,cidx2) );
	#endif
	                    }

	                    if (newThisRhs<=-EPS) /* cidx1 = 0 and cidx2 = 0 */
	                    {
	                    	inactivePairwise++;
	                        cvec.push_back( pair<int,int>(cidx1+nCols,cidx2+nCols) );
	#ifdef DEBUG_CONF
	                        newConflicts++;
	                        confS.push_back( pair<int,int>(cidx1+nCols,cidx2+nCols) );
	#endif
	                    }
	                }
	                fetchConflicts( cvec, false, cgraph, neighs );
	            }
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

    fetchConflicts( cvec, true, cgraph, neighs );

    return cgraph;
}

void processClique( const int n, const int *idx, CGraph *cgraph, vector< pair<int, int> > &cvec, const double *colLb, const double *colUb )
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

void fetchConflicts( vector< pair<int,int> > &cvec, const bool lastTime, CGraph *cgraph, vector<int> &neighs )
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