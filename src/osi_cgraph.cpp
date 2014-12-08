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

/* greedy clique partitioning from "Conflict graphs in solving integer programming problems" */
vector<vector<int> > greedyCliquePartitioning(const CGraph *cgraph, const int nElements, const int *idxs, const double *coefs);

/* Returns the lower bound for LHS by fixing two variables (uses clique partitioning in partitions) */
double getLr(const CGraph *cgraph, const vector<vector<int> > &partitions, const vector<double> &coefs,
             const vector<int> &whichPartition, const double sumPosCoefs, const int x1, const int x2);

/* generates a maximal clique containing vertex */
/* candidates is the initial list of candidates nodeCoef */
vector<int> getMaximalClique(const CGraph *cgraph, const int vertex, const vector<int>& candidates);

/* uses clique partitioning to extend the conflict graph row-by-row */
vector<pair<int, int> > extendConflictGraphByRow(CGraph *cgraph, const int numElements, const int *idxs, const double *coefs, const double rhs);

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
        
        if ( (nElements<2) || (fabs(rhs[idxRow])>=LARGE_CONST) )
            continue;

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
            else nPos++;

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
        if(nBools < nElements || (nPos == nElements && fabs(rhs[idxRow]) <= EPS))
            continue;

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

#ifdef DEBUG_CONF
        if (newConflicts)
        {
            printf("%s,%c,%g : ", lp->getRowName(idxRow).c_str(), sense[idxRow], rhs[idxRow]);
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
    }

    fetchConflicts(true, cgraph);
    cgraph_update_min_max_degree( cgraph );

    return cgraph;
}

CGraph *osi_build_cgraph( void *_lp )
{
    clock_t start = clock();
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
    double minCoefs[nRows], maxCoefs[nRows];
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

        minCoefs[idxRow] = DBL_MAX;
        maxCoefs[idxRow] = -DBL_MAX;
        
        if ( (nElements<2) || (fabs(rhs[idxRow])>=LARGE_CONST) )
            continue;

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

            minCoefs[idxRow] = min(minCoefs[idxRow], coefs[i]);
        	maxCoefs[idxRow] = max(maxCoefs[idxRow], coefs[i]);

            if(ctype[cidx] == 1)
                nBools++;

            if (FIXED_IN_ZERO(cidx))
                continue;

            if(columns[i].second <= -EPS)
                sumNegCoefs += columns[i].second;
            else nPos++;

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
        if(nBools < nElements || (nPos == nElements && fabs(rhs[idxRow]) <= EPS))
            continue;

        sort(columns.begin(), columns.end(), sort_columns());

        if(nPos < nElements || nElements <= MIN_PAIRWISE_ANALYSIS)
        {
            bool groupingWorked = pairwiseAnalysisByGrouping(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);
            if(!groupingWorked)
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
                groupingWorked = pairwiseAnalysisByGrouping(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);
                if(!groupingWorked)
                    pairwiseAnalysis(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);
            }   
        }

        else
        {
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

                bool foundClique = cliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);

                if(!foundClique)/* normal case, checking pairwise conflicts
                     * (these are computed increasing RHS in thisRhs) */
                {
                    bool groupingWorked = pairwiseAnalysisByGrouping(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);

                    if(!groupingWorked)
                        pairwiseAnalysis(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);
                }            
            }
        }

#ifdef DEBUG_CONF
        if (newConflicts)
        {
            printf("%s,%c,%g : ", lp->getRowName(idxRow).c_str(), sense[idxRow], rhs[idxRow]);
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
    }

    fetchConflicts(true, cgraph);

    double firstTime = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    unsigned long int firstGraph = 0;
    for(int i = 0; i < cgraph_size( cgraph ); i++)
        firstGraph += cgraph_degree(cgraph, i);
    firstGraph /= 2;

    start = clock();
    int bestRow = 0, bestNz = 0, bestConfs = 0;
    for(idxRow = 0; idxRow < nRows; idxRow++)
    {
        const CoinShallowPackedVector &row = M->getVector(idxRow);
        const int nElements = row.getNumElements();
        const int *idxs = row.getIndices();
        double *coefs = (double*) row.getElements();
        double thisRhs = rhs[idxRow];

        if(nElements > MAX_NONZEROS)
            continue;

        int maxDegree = 0;
        for(int i = 0; i < nElements; i++)
        {
            maxDegree = max(maxDegree, cgraph_degree(cgraph, idxs[i]));
            maxDegree = max(maxDegree, cgraph_degree(cgraph, idxs[i]+nCols));
        }

        if(maxDegree < 2)
        	continue;

        vector<pair<int, int> > test;

        if(sense[idxRow] == 'L')
        {
            test = extendConflictGraphByRow(cgraph, nElements, idxs, coefs, thisRhs);
            if((int)test.size() > bestConfs)
	        {
	        	bestConfs = (int)test.size();
	        	bestRow = idxRow;
	        	bestNz = nElements;
	        }
        }

        else if(sense[idxRow] == 'G')
        {
            thisRhs = -1.0 * thisRhs;
            for(int i = 0; i < nElements; i++)
                coefs[i] = -1.0 * coefs[i];
            test = extendConflictGraphByRow(cgraph, nElements, idxs, coefs, thisRhs);
            if((int)test.size() > bestConfs)
	        {
	        	bestConfs = (int)test.size();
	        	bestRow = idxRow;
	        	bestNz = nElements;
	        }
        }

        else if(sense[idxRow] == 'E')
        {
            test = extendConflictGraphByRow(cgraph, nElements, idxs, coefs, thisRhs);
            if((int)test.size() > bestConfs)
	        {
	        	bestConfs = (int)test.size();
	        	bestRow = idxRow;
	        	bestNz = nElements;
	        }
            thisRhs = -1.0 * thisRhs;
            for(int i = 0; i < nElements; i++)
                coefs[i] = -1.0 * coefs[i];
            test = extendConflictGraphByRow(cgraph, nElements, idxs, coefs, thisRhs);
            if((int)test.size() > bestConfs)
	        {
	        	bestConfs = (int)test.size();
	        	bestRow = idxRow;
	        	bestNz = nElements;
	        }
        }
        fetchConflicts(true, cgraph);
    }

    double secondTime = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    unsigned long int secondGraph = 0;
    for(int i = 0; i < cgraph_size( cgraph ); i++)
        secondGraph += cgraph_degree(cgraph, i);
    secondGraph /= 2;

    printf("%.2lf %lu %.2lf %lu\n", firstTime, firstGraph, firstTime+secondTime, secondGraph);
    printf("%s \t Nzs: %d \t Conflicts: %d\n", lp->getRowName(bestRow).c_str(), bestNz, bestConfs);

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
            {
                cvec.push_back( pair<int,int>(cidx1, cidx2) );
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1, cidx2) );
#endif
            }

            if(coef1 + negDiscount > rhs + 0.001) /* cidx1 = 1 and cidx2 = 0 */
            {
                cvec.push_back( pair<int,int>(cidx1, cidx2+nCols) );
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1, cidx2+nCols) );
#endif
            }

            if(coef2 + negDiscount > rhs + 0.001) /* cidx1 = 0 and cidx2 = 1 */
            {
                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2) );
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1+nCols, cidx2) );
#endif
            }

            if(negDiscount > rhs + 0.001) /* cidx1 = 0 and cidx2 = 0 */
            {
                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2+nCols) );
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
    int cliqueCompSize = 0, cliqueSize = 0;

    L00 = sumNegCoefs - min(0.0, columns[1].second) - min(0.0, columns[0].second);
    L11 = sumNegCoefs - min(0.0, columns[nElements-2].second) - min(0.0, columns[nElements-1].second)
          + columns[nElements-2].second + columns[nElements-1].second;

    if(L00 <= rhs && L11 <= rhs)
        return false;

    #define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )

    if(L00 > rhs + 0.001)
    {
        for(int i = nElements - 1; i > 0; i--)
        {
            double D = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[i-1].second);

            if(D > rhs + 0.001)
            {
                cliqueCompStart = i;
                break;
            }
        }
        assert(cliqueCompStart > 0 && cliqueCompStart < nElements);
        int n = cliqueCompStart + 1, idxs[n];
        for(int i = 0; i < n; i++)
        {
            if (FIXED_IN_ZERO(columns[i].first))
                continue;
            idxs[i] = columns[i].first + nCols; //binary complement
            cliqueCompSize++;
        }
        if(cliqueCompSize > 2)
        {
            processClique( cliqueCompSize, (const int *)idxs, cgraph, colLb, colUb );

            //looking for the same coefficient. For example: x + 2y + 2z + 3w >= 4 (- x - 2y - 2z - 3w <= -4)
            for(int j = cliqueCompStart + 1; j < nElements; j++)
            {
                if(fabs(columns[j].second - columns[cliqueCompStart].second) <= EPS)
                {
                	if(!FIXED_IN_ZERO(columns[j].first))
                	{
	                    idxs[n-1] = columns[j].first;
	                    processClique( cliqueCompSize, (const int *)idxs, cgraph, colLb, colUb );
	                    cliqueCompStart++;
	                }
                }
                else break;
            }
        }
    }

    if(L11 > rhs + 0.001)
    {
        for(int i = 0; i < nElements - 1; i++)
        {
            double D = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[i+1].second);
            double L = D + columns[i].second + columns[i+1].second;

            if(L > rhs + 0.001)
            {
                cliqueStart = i;
                break;
            }
        }
        
        assert(cliqueStart >= 0 && cliqueStart < nElements - 1);
        int n = nElements - cliqueStart, idxs[n];
        for(int i = cliqueStart, j = 0; i < nElements; i++)
        {
            if (FIXED_IN_ZERO(columns[i].first))
                continue;
            idxs[j++] = columns[i].first;
            cliqueSize++;
        }
        if(cliqueSize > 2)
        {
            processClique( cliqueSize, (const int *)idxs, cgraph, colLb, colUb );

            //looking for the same coefficient. For example: x + y + 2z + 2w <= 2
            for(int j = cliqueStart - 1; j >= 0; j--)
            {
                if(fabs(columns[j].second - columns[cliqueStart].second) <= EPS && !FIXED_IN_ZERO(columns[j].first))
                {
                    idxs[0] = columns[j].first;
                    processClique( cliqueSize, (const int *)idxs, cgraph, colLb, colUb );
                    cliqueStart--;
                }
                else break;
            }
        }
    }

    //checking for conflicts between variables that didnt appear in any clique
    int last = max(cliqueCompStart, cliqueStart);
    	last = max(last, 0); //just to handle the case when cliqueCompStart and cliqueStart are equal to -1
    for(int i = 0; i < last; i++)
    {
        if(FIXED_IN_ZERO(columns[i].first))
            continue;
        for(int j = i+1; j < nElements; j++)
        {
            if(FIXED_IN_ZERO(columns[j].first))
                continue;

            double negDiscount = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[j].second);
            int cidx1 = columns[i].first, cidx2 = columns[j].first;
            double coef1 = columns[i].second, coef2 = columns[j].second;

            if(coef1 + coef2 + negDiscount > rhs + 0.001)
            {
                cvec.push_back( pair<int,int>(cidx1, cidx2) );
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1, cidx2) );
#endif
            }

            if(coef1 + negDiscount > rhs + 0.001) /* cidx1 = 1 and cidx2 = 0 */
            {
                cvec.push_back( pair<int,int>(cidx1, cidx2+nCols) );
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1, cidx2+nCols) );
#endif
            }

            if(coef2 + negDiscount > rhs + 0.001) /* cidx1 = 0 and cidx2 = 1 */
            {
                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2) );
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1+nCols, cidx2) );
#endif
            }

            if(negDiscount > rhs + 0.001) /* cidx1 = 0 and cidx2 = 0 */
            {
                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2+nCols) );
#ifdef DEBUG_CONF
                newConflicts++;
                confS.push_back( pair<int,int>(cidx1+nCols, cidx2+nCols) );
#endif
            }
        }
        fetchConflicts( false, cgraph );
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

    //x1 e x2 estao em alguma particao
    if(part1 && part2)
        for(int i = 0; i < (int)partitions.size(); i++)
        {
            double coefSel = 0.0;

            if(part1 == i)
                coefSel = fabs(coef1);

            else if(part2 == i)
                coefSel = fabs(coef2);

            else 
                for(int j = 0; j < (int)partitions[i].size(); j++)
                {
                    const int var = partitions[i][j];
                    const int realVar = var%numCols;
                    coefSel = max(coefSel, fabs(coefs[realVar]));
                }

            Lr -= coefSel;
        }

    //x1 e x2 nao estao em nenhuma particao
    else if(!inPart1 && !inPart2)
        for(int i = 0; i < (int)partitions.size(); i++)
        {
            double coefSel = 0.0;
            for(int j = 0; j < (int)partitions[i].size(); j++)
            {
                const int var = partitions[i][j];
                const int realVar = var%numCols;
                if(cgraph_conflicting_nodes(cgraph, x1, var))
                    continue;
                if(cgraph_conflicting_nodes(cgraph, x2, var))
                    continue;
                coefSel = max(coefSel, fabs(coefs[realVar]));
            }
            Lr -= coefSel;
        }

    //x1 esta em uma particao e x2 nao
    else if(inPart1 && !inPart2)
    {
        for(int i = 0; i < (int)partitions.size(); i++)
        {
            double coefSel = 0.0;

            if(part1 == i)
                coefSel = fabs(coef1);

            else
            {
                for(int j = 0; j < (int)partitions[i].size(); j++)
                {
                    const int var = partitions[i][j];
                    const int realVar = var%numCols;
                    if(cgraph_conflicting_nodes(cgraph, x2, var))
                        continue;
                    coefSel = max(coefSel, fabs(coefs[realVar]));
                }
            }
            Lr -= coefSel;
        }
    }

    //x2 esta em uma particao e x1 nao
    else
        for(int i = 0; i < (int)partitions.size(); i++)
        {
            double coefSel = 0.0;

            if(part2 == i)
                coefSel = fabs(coef2);

            else
                for(int j = 0; j < (int)partitions[i].size(); j++)
                {
                    const int var = partitions[i][j];
                    const int realVar = var%numCols;
                    if(cgraph_conflicting_nodes(cgraph, x1, var))
                        continue;
                    coefSel = max(coefSel, fabs(coefs[realVar]));
                }
            Lr -= coefSel;
        }

    return Lr;
}

vector<pair<int, int> > extendConflictGraphByRow(CGraph *cgraph, const int nElements, const int *idxs, const double *coefs, const double rhs)
{
    const int numCols = cgraph_size(cgraph) / 2;
    double sumPosCoefs = 0.0;
    vector<double> rowCoefs(numCols, 0.0);
    vector<int> whichPartition(numCols, -1);
    vector<vector<int> >  partitions = greedyCliquePartitioning(cgraph, nElements, idxs, coefs);
    vector<pair<int, int> > newConf;

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
                {
                    cvec.push_back( pair<int,int>(cidx1, cidx2) );
                    newConf.push_back( pair<int,int>(cidx1, cidx2) );
#ifdef DEBUG_CONF
                    newConflicts++;
                    confS.push_back( pair<int,int>(cidx1, cidx2) );
#endif
                }
            }

            if(!cgraph_conflicting_nodes(cgraph, cidx1+numCols, cidx2))
            {
                double Lr = getLr(cgraph, partitions, rowCoefs, whichPartition, sumPosCoefs, cidx1+numCols, cidx2);
                if(Lr > rhs + 0.001)
                {
                    cvec.push_back( pair<int,int>(cidx1+numCols, cidx2) );
                    newConf.push_back( pair<int,int>(cidx1+numCols, cidx2) );
#ifdef DEBUG_CONF
                    newConflicts++;
                    confS.push_back( pair<int,int>(cidx1+numCols, cidx2) );
#endif
                }
            }

            if(!cgraph_conflicting_nodes(cgraph, cidx1, cidx2+numCols))
            {
                double Lr = getLr(cgraph, partitions, rowCoefs, whichPartition, sumPosCoefs, cidx1, cidx2+numCols);
                if(Lr > rhs + 0.001)
                {
                    cvec.push_back( pair<int,int>(cidx1, cidx2+numCols) );
                    newConf.push_back( pair<int,int>(cidx1, cidx2+numCols) );
#ifdef DEBUG_CONF
                    newConflicts++;
                    confS.push_back( pair<int,int>(cidx1, cidx2+numCols) );
#endif
                }
            }

            if(!cgraph_conflicting_nodes(cgraph, cidx1+numCols, cidx2+numCols))
            {
                double Lr = getLr(cgraph, partitions, rowCoefs, whichPartition, sumPosCoefs, cidx1+numCols, cidx2+numCols);
                if(Lr > rhs + 0.001)
                {
                    cvec.push_back( pair<int,int>(cidx1+numCols, cidx2+numCols) );
                    newConf.push_back( pair<int,int>(cidx1+numCols, cidx2+numCols) );
#ifdef DEBUG_CONF
                    newConflicts++;
                    confS.push_back( pair<int,int>(cidx1+numCols, cidx2+numCols) );
#endif
                }
            }
        }
    }
     #undef FIXED_IN_ZERO

    return newConf;
}