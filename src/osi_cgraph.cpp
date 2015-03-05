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

#define MIN_PAIRWISE_ANALYSIS 100

#define MAX_NONZEROS 128

const double LARGE_CONST = std::min( DBL_MAX/10.0, 1e20 );

vector< pair< int, int > > cvec;/* conflict vector */
vector< int > neighs;/* used i fetch conflicts */
double *colLb;
double *colUb;
int nCols;
int nRows;

struct sort_sec_pair
{
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right)
    {
        if ( left.second != right.second )
            return left.second < right.second;

        return left.first < right.first;
    }
};

struct sort_sec_pair_reverse
{
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right)
    {
        if ( left.second != right.second )
            return left.second > right.second;

        return left.first < right.first;
    }
};

struct sort_columns
{
    bool operator()(const std::pair<int, double> &left, const std::pair<int,double> &right)
    {
        if ( fabs(left.second - right.second) > EPS )
            return ( left.second < right.second );

        return left.first < right.first;
    }
};

struct sort_columns_reverse
{
    bool operator()(const std::pair<int, double> &left, const std::pair<int,double> &right)
    {
        if ( fabs(left.second - right.second) > EPS )
            return ( left.second > right.second );

        return left.first < right.first;
    }
};

/* Searches for a clique in this constraint. */
void cliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs);

/* Searches for a clique involving the complement of variables in this constraint. */
void cliqueComplementDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs);

bool pairwiseAnalysisByGrouping(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs);

void pairwiseAnalysis(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs);

void processClique( const int n, const int *idx, CGraph *cgraph, const double *colLb, const double *colUb );

/* if size is large enough fech conflicts */
void fetchConflicts( const bool lastTime, CGraph *cgraph );

/* returns how much a variable can negatively contribute */
double mostNegativeContribution( const double coef, const char s, const double colLb, const double colUb );

/* returns how much a variable can positively contribute */
double unitaryContribution( const double coef, const char s, const double colLb, const double colUb );

/* greedy clique partitioning from "Conflict graphs in solving integer programming problems" */
vector<vector<int> > greedyCliquePartitioning(const CGraph *cgraph, const int nElements, const int *idxs, const double *coefs);

/* Returns the lower bound for LHS by fixing two variables (uses clique partitioning in partitions) */
double getLr(const CGraph *cgraph, const vector<vector<int> > &partitions, const vector<double> &coefs,
             const vector<int> &whichPartition, const double sumPosCoefs, const int x1, const int x2);

/* generates a maximal clique containing vertex */
/* candidates is the initial list of candidates nodeCoef */
vector<int> getMaximalClique(const CGraph *cgraph, const int vertex, const vector<int>& candidates);

/* uses clique partitioning to extend the conflict graph row-by-row */
void extendConflictGraphByRow(CGraph *cgraph, const int numElements, const int *idxs, const double *coefs, const double rhs);

/* Returns the first position of columns which the lower bound for LHS (considering activation of variables) is greater than rhs */
/* l=initial position for search in columns u=last position for search in columns */
int binary_search(const vector< pair<int, double> >& columns, double rhs, double coef, double sumNegCoefs, int l, int u);

/* Returns the first position of columns which the lower bound for LHS (considering deactivation of variables) is greater than rhs */
/* l=initial position for search in columns u=last position for search in columns */
int binary_search_complement(const vector< pair<int, double> >& columns, double rhs, double coef, double sumNegCoefs, int l, int u);

CGraph *osi_build_cgraph_pairwise( void *_lp )
{
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
        int nBools = 0, nPos = 0; // number of binary variables
        double sumNegCoefs = 0.0; //sum of all negative coefficients
        double minCoef = numeric_limits<double>::max();
      	double maxCoef = numeric_limits<double>::min();
        
        if ( (nElements<2) || (fabs(rhs[idxRow])>=LARGE_CONST) )
            continue;

        if ( sense[idxRow] == 'R' )  // lets not consider ranged constraints by now
        {
            printf("TODO: CHECK FOR RANGED CONSTRAINT (%s) rhs is %g\n", lp->getRowName(idxRow).c_str(), rhs[idxRow] );
            continue;
        }

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
            else nPos++;

            minCoef = min(minCoef, columns[i].second);
         	maxCoef = max(maxCoef, columns[i].second);

            /* inserting trivial conflicts: variable-complement */
            if(ctype[cidx] != 1) //consider only binary variables
                continue;
            cvec.push_back( pair<int, int>(cidx, cidx + nCols) );
        }

        /* considering just constraints which have only binary variables */
        if(nBools < nElements || (nPos == nElements && fabs(rhs[idxRow]) <= EPS))
            continue;

      	/* special case: GUB constraints */
		if ( DBL_EQUAL( minCoef, maxCoef ) &&  DBL_EQUAL( maxCoef, rhs[idxRow] * mult ) &&
		    DBL_EQUAL(minCoef, 1.0) && ((sense[idxRow]=='E') || (sense[idxRow]=='L'))
            && (row.getNumElements() > 3) ) 
		{
		    processClique( row.getNumElements(), (const int *)idx, cgraph, colLb, colUb );
		}

		else
		{
	        pairwiseAnalysis(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);

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
	            pairwiseAnalysis(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);            
	        }
	    }
	    
#undef FIXED_IN_ZERO
    }

    fetchConflicts(true, cgraph);
    cgraph_update_min_max_degree( cgraph );

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
                cvec.push_back( pair<int,int>(cidx1, cidx2) );

            if(coef1 + negDiscount > rhs + 0.001) /* cidx1 = 1 and cidx2 = 0 */
                cvec.push_back( pair<int,int>(cidx1, cidx2+nCols) );

            if(coef2 + negDiscount > rhs + 0.001) /* cidx1 = 0 and cidx2 = 1 */
                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2) );

            if(negDiscount > rhs + 0.001) /* cidx1 = 0 and cidx2 = 0 */
                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2+nCols) );
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
                    }
                }
            }
        }
        fetchConflicts(false, cgraph);
    }
    #undef FIXED_IN_ZERO
    return true;
}

vector<vector<int> > greedyCliquePartitioning(const CGraph *cgraph, const int nElements, const int *idxs, const double *coefs)
{
    const int numCols = cgraph_size(cgraph) / 2;
    int idxMap[numCols*2];
    vector<vector<int> > partition;
    vector<pair<int, double> > nodeCoef(nElements);
    vector<bool> marked(numCols*2, false);

    for(int i = 0; i < nElements; i++)
    {
        if(coefs[i] <= -EPS)
            nodeCoef[i] = pair<int, double>(idxs[i], fabs(coefs[i]));
        else
            nodeCoef[i] = pair<int, double>(numCols + idxs[i], coefs[i]);
    }

    sort(nodeCoef.begin(), nodeCoef.end(), sort_columns_reverse());
    fill(idxMap, idxMap + (numCols*2), -1);
    for(int i = 0; i < nElements; i++)
    	idxMap[nodeCoef[i].first] = i;

    for(int i = 0; i < nElements; i++)
    {
        if(!marked[i])
        {
        	vector<int> candidates;
        	for(int j = i+1; j < nElements; j++)
        		if(!marked[j])
        			candidates.push_back(nodeCoef[j].first);

            vector<int> clique = getMaximalClique(cgraph, nodeCoef[i].first, candidates);
            vector<int> tmp(clique.size());
            for(int j = 0; j < (int)clique.size(); j++)
            {
            	int mIdx = idxMap[clique[j]];
                tmp[j] = clique[j];
                marked[mIdx] = true;
            }
            partition.push_back(tmp);
        }
    }
    
    return partition;
}

vector<int> getMaximalClique(const CGraph *cgraph, const int vertex, const vector<int>& candidates)
{
    int numVertices = cgraph_size(cgraph);
    vector<int> clique;
    vector<pair<int, int> > cList;
    vector<bool> used(numVertices, false);

    clique.push_back(vertex);
    used[vertex] = true;

    for(int i = 0; i < (int)candidates.size(); i++)
        if(cgraph_conflicting_nodes(cgraph, vertex, candidates[i]))
            cList.push_back(pair<int, int>(candidates[i], cgraph_degree(cgraph, candidates[i])));

    sort(cList.begin(), cList.end(), sort_sec_pair_reverse());

    for(int i = 0; i < (int)cList.size(); i++)
    {
        if(!used[cList[i].first])
        {
            clique.push_back(cList[i].first);
            used[cList[i].first] = true;
            for(int j = i + 1; j < (int)cList.size(); j++)
                if(!cgraph_conflicting_nodes(cgraph, cList[i].first, cList[j].first))
                    used[cList[j].first] = true;
        }
    }

    return clique;
}

double getLr(const CGraph *cgraph, const vector<vector<int> > &partitions, const vector<double> &coefs,
             const vector<int> &whichPartition, const double sumPosCoefs, const int x1, const int x2)
{
    const int numCols = cgraph_size(cgraph) / 2;
    double Lr = sumPosCoefs;
    double coef1, coef2;
    int realX1, realX2, part1, part2;
    bool inPart1, inPart2; //checks if variables are in some partition
    /* variables with negative coefficients or complement of variables with positive coefficientes are in some partition. */

    if(cgraph_conflicting_nodes(cgraph, x1, x2)) //conflict already exists
        return DBL_MAX/2;

    realX1 = x1 % numCols;
    realX2 = x2 % numCols;
    coef1 = coefs[realX1];
    coef2 = coefs[realX2];
    inPart1 = (coef1 <= -EPS || x1 >= numCols); 
    inPart2 = (coef2 <= -EPS || x2 >= numCols); 
    part1 = whichPartition[realX1];
    part2 = whichPartition[realX2];

    assert(fabs(coef1) > EPS && fabs(coef2) > EPS);
    assert(part1 >= 0 && part2 >= 0);

    //x1 and x2 are in some partition
    if(part1 && part2)
        for(int i = 0; i < (int)partitions.size(); i++)
        {
            double coefSel = 0.0;

            if(part1 == i)
                coefSel = fabs(coef1);

            else if(part2 == i)
                coefSel = fabs(coef2);

            else
            { 
                const int var = partitions[i][0];
                const int realVar = var%numCols;
                coefSel = fabs(coefs[realVar]);
            }

            Lr -= coefSel;
        }

    //x1 and x2 are not in any partition
    else if(!inPart1 && !inPart2)
    {
        for(int i = 0; i < (int)partitions.size(); i++)
        {
            double coefSel = 0.0;
            int var = partitions[i][0];
            int realVar = var%numCols;

            if(!cgraph_conflicting_nodes(cgraph, x1, var) && !cgraph_conflicting_nodes(cgraph, x2, var))
                coefSel = fabs(coefs[realVar]);
            
            else
                for(int j = 1; j < (int)partitions[i].size(); j++)
                {
                    var = partitions[i][j];
                    realVar = var%numCols;
                    if(cgraph_conflicting_nodes(cgraph, x1, var))
                        continue;
                    if(cgraph_conflicting_nodes(cgraph, x2, var))
                        continue;
                    coefSel = max(coefSel, fabs(coefs[realVar]));
                }
            Lr -= coefSel;
        }
    }

    //x1 is in some partition and x2 not
    else if(inPart1 && !inPart2)
    {
        for(int i = 0; i < (int)partitions.size(); i++)
        {
            double coefSel = 0.0;
            int var = partitions[i][0];
            int realVar = var%numCols;

            if(part1 == i)
                coefSel = fabs(coef1);

            else
            {
                int var = partitions[i][0];
                int realVar = var%numCols;

                if(!cgraph_conflicting_nodes(cgraph, x2, var))
                    coefSel = fabs(coefs[realVar]);
                else
                    for(int j = 1; j < (int)partitions[i].size(); j++)
                    {
                        var = partitions[i][j];
                        realVar = var%numCols;
                        if(cgraph_conflicting_nodes(cgraph, x2, var))
                            continue;
                        coefSel = max(coefSel, fabs(coefs[realVar]));
                    }
            }
            Lr -= coefSel;
        }
    }

    //x2 is in some partition and x1 not
    else
        for(int i = 0; i < (int)partitions.size(); i++)
        {
            double coefSel = 0.0;

            if(part2 == i)
                coefSel = fabs(coef2);

            else
            {
                int var = partitions[i][0];
                int realVar = var%numCols;

                if(!cgraph_conflicting_nodes(cgraph, x1, var))
                    coefSel = fabs(coefs[realVar]);
                else
                    for(int j = 1; j < (int)partitions[i].size(); j++)
                    {
                        const int var = partitions[i][j];
                        const int realVar = var%numCols;
                        if(cgraph_conflicting_nodes(cgraph, x1, var))
                            continue;
                        coefSel = max(coefSel, fabs(coefs[realVar]));
                    }
            }
            Lr -= coefSel;
        }

    return Lr;
}

void extendConflictGraphByRow(CGraph *cgraph, const int nElements, const int *idxs, const double *coefs, const double rhs)
{
    const int numCols = cgraph_size(cgraph) / 2;
    double sumPosCoefs = 0.0;
    vector<double> rowCoefs(numCols, 0.0);
    vector<int> whichPartition(numCols, -1);
    vector<vector<int> >  partitions = greedyCliquePartitioning(cgraph, nElements, idxs, coefs);

    for(int i = 0; i < nElements; i++)
    {
        rowCoefs[idxs[i]%numCols] = coefs[i];
        if(coefs[i] >= EPS)
            sumPosCoefs += coefs[i];
    }        

    for(int i = 0; i < (int)partitions.size(); i++)
        for(int j = 0; j < (int)partitions[i].size(); j++)
        {
            int var = partitions[i][j];
            whichPartition[var%numCols] = i;
        }

    for(int j1 = 0; j1 < nElements; j1++)
    {
        const int cidx1 = idxs[j1];
        const double coef1 = coefs[j1];

        #define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )

        if (FIXED_IN_ZERO(cidx1))
            continue;

        for(int j2 = j1+1; j2 < nElements; j2++)
        {
            const int cidx2 = idxs[j2];
            const double coef2 = coefs[j2];

            if(FIXED_IN_ZERO(cidx2))
                continue;

            if(!cgraph_conflicting_nodes(cgraph, cidx1, cidx2))
            {
                double Lr = getLr(cgraph, partitions, rowCoefs, whichPartition, sumPosCoefs, cidx1, cidx2);
                if(Lr > rhs + 0.001)
                    cvec.push_back( pair<int,int>(cidx1, cidx2) );
            }

            if(!cgraph_conflicting_nodes(cgraph, cidx1+numCols, cidx2))
            {
                double Lr = getLr(cgraph, partitions, rowCoefs, whichPartition, sumPosCoefs, cidx1+numCols, cidx2);
                if(Lr > rhs + 0.001)
                    cvec.push_back( pair<int,int>(cidx1+numCols, cidx2) );
            }

            if(!cgraph_conflicting_nodes(cgraph, cidx1, cidx2+numCols))
            {
                double Lr = getLr(cgraph, partitions, rowCoefs, whichPartition, sumPosCoefs, cidx1, cidx2+numCols);
                if(Lr > rhs + 0.001)
                    cvec.push_back( pair<int,int>(cidx1, cidx2+numCols) );
            }

            if(!cgraph_conflicting_nodes(cgraph, cidx1+numCols, cidx2+numCols))
            {
                double Lr = getLr(cgraph, partitions, rowCoefs, whichPartition, sumPosCoefs, cidx1+numCols, cidx2+numCols);
                if(Lr > rhs + 0.001)
                    cvec.push_back( pair<int,int>(cidx1+numCols, cidx2+numCols) );
            }
        }
    }
     #undef FIXED_IN_ZERO
}

int binary_search(const vector< pair<int, double> >& columns, double rhs, double coef, double sumNegCoefs, int l, int u)
{
	int mid, position;
	if(l <= u)
	{
		mid = (l + u) / 2;
		double D = sumNegCoefs - min(0.0, coef) - min(0.0, columns[mid].second);
        double LHS = D + coef + columns[mid].second;
	  	if(rhs >= LHS)
	    	position = binary_search(columns, rhs, coef, sumNegCoefs, mid + 1, u);
	  	else
	    	position = binary_search(columns, rhs, coef, sumNegCoefs, l, mid-1);
	}
	else
	{
	 //Check boundaries
	 if(l < u || u < 0)
	    u++;//Fix boundaries
	 position =  u + 1;
	}
	return position;
}

int binary_search_complement(const vector< pair<int, double> >& columns, double rhs, double coef, double sumNegCoefs, int l, int u)
{
	int mid, position;
	if(l <= u)
	{
		mid = (l + u) / 2;
		double LHS = sumNegCoefs - min(0.0, coef) - min(0.0, columns[mid].second);
	  	if(rhs >= LHS)
	    	position = binary_search_complement(columns, rhs, coef, sumNegCoefs, mid + 1, u);
	  	else
	    	position = binary_search_complement(columns, rhs, coef, sumNegCoefs, l, mid-1);
	}
	else
	{
	 //Check boundaries
	 if(l < u || u < 0)
	    u++;//Fix boundaries
	 position =  u + 1;
	}
	return position;
}

CGraph *osi_build_cgraph( void *_lp )
{
    OsiSolverInterface *lp = (OsiSolverInterface *)_lp;

    if (lp->getNumIntegers()<2)
        return 0;

    int cgraphSize = lp->getNumCols() * 2;
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
        int nPos = 0; //number of positive coefficients
        double sumNegCoefs = 0.0; //sum of all negative coefficients
        double minCoef = numeric_limits<double>::max();
        double maxCoef = numeric_limits<double>::min();
        
        if ( (nElements<2) || (fabs(rhs[idxRow])>=LARGE_CONST) )
            continue;

        if ( sense[idxRow] == 'R' )  // lets not consider ranged constraints by now
        {
            printf("TODO: CHECK FOR RANGED CONSTRAINT (%s) rhs is %g\n", lp->getRowName(idxRow).c_str(), rhs[idxRow] );
            continue;
        }

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
            else nPos++;

            minCoef = min(minCoef, columns[i].second);
            maxCoef = max(maxCoef, columns[i].second);

            /* inserting trivial conflicts: variable-complement */
            if(ctype[cidx] != 1) //consider only binary variables
                continue;
            cvec.push_back( pair<int, int>(cidx, cidx + nCols) );
        }

        /* considering just constraints which have only binary variables */
        if(nBools < nElements || (nPos == nElements && fabs(rhs[idxRow]) <= EPS))
            continue;

        /* special case: GUB constraints */
        if ( DBL_EQUAL( minCoef, maxCoef ) &&  DBL_EQUAL( maxCoef, rhs[idxRow] * mult ) &&
            DBL_EQUAL(minCoef, 1.0) && ((sense[idxRow]=='E') || (sense[idxRow]=='L'))
            && (row.getNumElements() > 3) ) 
        {
            processClique( row.getNumElements(), (const int *)idx, cgraph, colLb, colUb );
        }

        else
        {
            sort(columns.begin(), columns.end(), sort_columns());
            cliqueDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);
            cliqueComplementDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);

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

                cliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);
                cliqueComplementDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);
            }
        }
    }

    fetchConflicts(true, cgraph);
    cgraph_update_min_max_degree( cgraph );

    return cgraph;
}

void cliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs)
{
    int nElements = (int)columns.size(), cliqueStart = -1;
    double maxLHS; //maxLHS = lower bound for LHS when the two variaveis with highest coefficients are activated.
    int cliqueSize = 0;

    maxLHS = sumNegCoefs - min(0.0, columns[nElements-2].second) - min(0.0, columns[nElements-1].second)
          + columns[nElements-2].second + columns[nElements-1].second;

    if(maxLHS <= rhs + EPS) return; //there is no clique involving activation of variables in this constraint.

    for(int i = 0; i < nElements - 1; i++)
    {
        double D = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[i+1].second);
        double LHS = D + columns[i].second + columns[i+1].second;

        if(LHS > rhs + EPS)
        {
            cliqueStart = i;
            break;
        }
    }

    assert(cliqueStart >= 0 && cliqueStart < nElements - 1);
    int n = nElements - cliqueStart, idxs[n];
    #define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )
    for(int i = cliqueStart, j = 0; i < nElements; i++)
    {
        if (FIXED_IN_ZERO(columns[i].first))
            continue;
        idxs[j++] = columns[i].first;
        cliqueSize++;
    }
    //process the first clique found
    processClique( cliqueSize, (const int *)idxs, cgraph, colLb, colUb );

    //now we have to check the variables that are outside of the clique found.
    for(int i = cliqueStart - 1; i >= 0; i--)
    {
    	int idx = columns[i].first;
    	double coef = columns[i].second;
    	int position = binary_search(columns, rhs, coef, sumNegCoefs, cliqueStart, nElements - 1);

		if(position < nElements) //clique was found
		{
			int n = nElements - position + 1, idxs[n];
			cliqueSize = 1;
			idxs[0] = idx;
		    for(int i = position, j = 1; i < nElements; i++)
		    {
		        if (FIXED_IN_ZERO(columns[i].first))
		            continue;
		        idxs[j++] = columns[i].first;
		        cliqueSize++;
		    }
		    processClique( cliqueSize, (const int *)idxs, cgraph, colLb, colUb );								
		}
		else break;
    }

    #undef FIXED_IN_ZERO
}

// void cliqueComplementDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs)
// {
//     int nElements = (int)columns.size(), cliqueCompStart = -1;
//     double maxLHS; //minLHS = lower bound for LHS when the two variaveis with smallest coefficients are deactivated.
//     int cliqueCompSize = 0;

//     maxLHS = sumNegCoefs - min(0.0, columns[1].second) - min(0.0, columns[0].second);

//     if(maxLHS <= rhs) return; //there is no clique involving the complement of variables in this constraint.

//     for(int i = nElements - 1; i > 0; i--)
//     {
//         double D = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[i-1].second);

//         if(D > rhs + 0.001)
//         {
//             cliqueCompStart = i;
//             break;
//         }
//     }

//     assert(cliqueCompStart > 0 && cliqueCompStart < nElements);
//     int n = cliqueCompStart + 1, idxs[n];
//     #define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )
//     for(int i = 0; i < n; i++)
//     {
//         if (FIXED_IN_ZERO(columns[i].first))
//             continue;
//         idxs[i] = columns[i].first + nCols; //binary complement
//         cliqueCompSize++;
//     }

//     //process the first clique found
//     processClique( cliqueCompSize, (const int *)idxs, cgraph, colLb, colUb );

//     //now we have to check the variables that are outside of the clique found.
//     for(int i = cliqueCompStart + 1; i < nElements; i++)
//     {
//     	printf("\n");
//     	int idx = columns[i].first;
//     	double coef = columns[i].second;
//     	int position = binary_search_complement(columns, rhs, coef, sumNegCoefs, 0, cliqueCompStart);

//     	printf("idx=%d coef=%.2lf position=%d\n", idx, coef, position);

// 		if(position <= cliqueCompStart)//clique was found
// 		{
// 			int n = position + 2, idxs[n];
// 			cliqueCompSize = 1;
// 			idxs[0] = idx + nCols;
// 		    for(int i = 0, j = 1; i <= position; i++)
// 		    {
// 		        if (FIXED_IN_ZERO(columns[i].first))
// 		            continue;
// 		        idxs[j++] = columns[i].first + nCols;
// 		        cliqueCompSize++;
// 		    }
// 		    processClique( cliqueCompSize, (const int *)idxs, cgraph, colLb, colUb );
// 		}
// 		else break;
//     }

//     #undef FIXED_IN_ZERO
// }

void cliqueComplementDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs)
{
    int nElements = (int)columns.size(), cliqueCompStart = -1;
    double maxLHS; //minLHS = lower bound for LHS when the two variaveis with smallest coefficients are deactivated.
    int cliqueCompSize = 0;
    vector<pair<int, double> > columnsInv(nElements);
    for(int i = 0; i < nElements; i++) columnsInv[i] = columns[nElements-i-1];

    maxLHS = sumNegCoefs - min(0.0, columnsInv[nElements-2].second) - min(0.0, columnsInv[nElements-1].second);

    if(maxLHS <= rhs + EPS) return; //there is no clique involving activation of variables in this constraint.

    for(int i = 0; i < nElements - 1; i++)
    {
        double LHS = sumNegCoefs - min(0.0, columnsInv[i].second) - min(0.0, columnsInv[i+1].second);

        if(LHS > rhs + EPS)
        {
            cliqueCompStart = i;
            break;
        }
    }

    assert(cliqueCompStart >= 0 && cliqueCompStart < nElements - 1);
    int n = nElements - cliqueCompStart, idxs[n];
    #define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )
    for(int i = cliqueCompStart, j = 0; i < nElements; i++)
    {
        if (FIXED_IN_ZERO(columnsInv[i].first))
            continue;
        idxs[j++] = columnsInv[i].first+nCols;
        cliqueCompSize++;
    }
    //process the first clique found
    processClique( cliqueCompSize, (const int *)idxs, cgraph, colLb, colUb );

    //now we have to check the variables that are outside of the clique found.
    for(int i = cliqueCompStart - 1; i >= 0; i--)
    {
    	int idx = columnsInv[i].first;
    	double coef = columnsInv[i].second;
    	int position = binary_search_complement(columnsInv, rhs, coef, sumNegCoefs, cliqueCompStart, nElements - 1);

		if(position < nElements) //clique was found
		{
			int n = nElements - position + 1, idxs[n];
			cliqueCompSize = 1;
			idxs[0] = idx+nCols;
		    for(int i = position, j = 1; i < nElements; i++)
		    {
		        if (FIXED_IN_ZERO(columnsInv[i].first))
		            continue;
		        idxs[j++] = columnsInv[i].first+nCols;
		        cliqueCompSize++;
		    }
		    processClique( cliqueCompSize, (const int *)idxs, cgraph, colLb, colUb );								
		}
		else break;
    }

    #undef FIXED_IN_ZERO
}