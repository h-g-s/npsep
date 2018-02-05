#include <cstring>
#include <cmath>
#include <limits>
#include <cfloat>
#include <algorithm>
#include <climits>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include "OsiSolverInterface.hpp"
#include "CoinTime.hpp"

extern "C"
{
#include "memory.h"
#include "node_heap.h"
#include "build_cgraph.h"
#include "lp.h"
}

using namespace std;

#define EPS 1e-8
#define MIN_CLIQUE_ROW 300 /* mininum size for a row to be considered a clique row */
#define DBL_EQUAL( v1, v2 ) ( fabs(v1-v2) < EPS )

// conflict vector dimensions
#define CVEC_CAP   1800000
#define CVEC_FLUSH 1700000

#define MAX_TIME_CG 60

const double LARGE_CONST = std::min( DBL_MAX/10.0, 1e20 );
double startCG;
char recomputeDegree = 0;

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
            return ( left.second < right.second );

        return left.first < right.first;
    }
};

//void pairwiseAnalysis(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs);

void processClique( const int n, const int *idx, CGraph *cgraph, vector< pair< int, int > > &cvec);

/* if size is large enough fetch conflicts */
void fetchConflicts( const bool lastTime, CGraph *cgraph, vector< pair< int, int > > &cvec);

/* Returns the first position of columns which the lower bound for LHS (considering activation of variables) is greater than rhs */
/* colStart=initial position for search in columns, colEnd=last position for search in columns */
/* partialLHS = LHS calculated with only one variable */
int binary_search(const vector< pair<int, double> >& columns, double partialLHS, double rhs, int colStart, int colEnd);

/* Returns the first position of columns which the lower bound for LHS (considering deactivation of variables) is greater than rhs */
/* colStart=initial position for search in columns, colEnd=last position for search in columns */
/* partialLHS = LHS calculated with only one variable */
int binary_search_complement(const vector< pair<int, double> >& columns, double partialLHS, double rhs, int colStart, int colEnd);

/* Searches for cliques involving the activation of variables in this constraint. */
void cliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs, vector< pair< int, int > > &cvec);

/* Searches for cliques involving the complement of variables in this constraint. */
void cliqueComplementDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs, vector< pair< int, int > > &cvec);

/* Searches for cliques involving variables and complements of variables in this constraint. */
void mixedCliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs, vector< pair< int, int > > &cvec);

//void pairwiseAnalysis(CGraph* cgraph, const vector<pair<int, double> >& columns, const double sumNegCoefs, const double rhs)
//{
//    int nElements = (int)columns.size();
//    int nCols = cgraph_size(cgraph) / 2;
//    for(int j1 = 0; j1 < nElements; j1++)
//    {
//        const int cidx1 = columns[j1].first;
//        const double coef1 = columns[j1].second;
//
//        for(int j2 = j1+1; j2 < nElements; j2++)
//        {
//            const int cidx2 = columns[j2].first;
//            const double coef2 = columns[j2].second;
//            const double negDiscount = sumNegCoefs - min(0.0, coef1) - min(0.0, coef2);
//
//            if(coef1 + coef2 + negDiscount > rhs + EPS)
//                cvec.push_back( pair<int,int>(cidx1, cidx2) );
//
//            if(coef1 + negDiscount > rhs + EPS) /* cidx1 = 1 and cidx2 = 0 */
//                cvec.push_back( pair<int,int>(cidx1, cidx2+nCols) );
//
//            if(coef2 + negDiscount > rhs + EPS) /* cidx1 = 0 and cidx2 = 1 */
//                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2) );
//
//            if(negDiscount > rhs + EPS) /* cidx1 = 0 and cidx2 = 0 */
//                cvec.push_back( pair<int,int>(cidx1+nCols, cidx2+nCols) );
//        }
//        fetchConflicts( false, cgraph );
//    }
//}

void processClique( const int n, const int *idx, CGraph *cgraph, vector< pair< int, int > > &cvec)
{
    if ( n >= MIN_CLIQUE_ROW )
    {
        vector< int > vidx;
        vidx.resize(n);

        memcpy( &(vidx[0]), idx, sizeof(int)*n );
        cgraph_add_clique( cgraph, &(vidx[0]), n );
        recomputeDegree = 1;
    }
    else
    {
        int i1, i2, nm1 = n-1;

        for ( i1=0 ; (i1<nm1) ; ++i1 )
            for ( i2=i1+1 ; (i2<n) ; ++i2 )
                cvec.push_back( pair<int,int>(idx[i1],idx[i2]) );
        fetchConflicts(false, cgraph, cvec);
    }
}

void fetchConflicts( const bool lastTime, CGraph *cgraph, vector< pair< int, int > > &cvec)
{
    if ( (!lastTime) && (cvec.size()<CVEC_FLUSH) )
        return;

    if ( cvec.empty() )
        return;
    pair< int, int > last;
    int currNode;
    vector< pair<int,int> >::const_iterator vIt;
    vector< int > neighs;
    neighs.reserve( 8192 );

    cvec.push_back( pair<int,int>( INT_MAX, INT_MAX ) );

    sort( cvec.begin(), cvec.end() );
    currNode = cvec.begin()->first;
    last = pair< int,int >( -1, -1 );

    for ( vIt=cvec.begin() ; (vIt!=cvec.end()) ; last=*vIt,++vIt )
    {
        if (last == *vIt)  // skipping repetitions
            continue;

        if ( vIt->first!=currNode )
        {
            if (neighs.size() > 0)
            {
                cgraph_add_node_conflicts_no_sim( cgraph, currNode, &(neighs[0]), (int)neighs.size() );
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
            if (neighs.size() > 0)
            {
                cgraph_add_node_conflicts_no_sim( cgraph, currNode, &(neighs[0]), (int)neighs.size() );
                neighs.clear();
            }
            currNode = vIt->second;
        }
        neighs.push_back( vIt->first );
    }

    cvec.clear();
}

int binary_search(const vector< pair<int, double> >& columns, double partialLHS, double rhs, int colStart, int colEnd)
{
    int mid;
    while(colStart <= colEnd)
    {
        mid = (colStart + colEnd) / 2;
        double LHS = partialLHS - min(0.0, columns[mid].second) + columns[mid].second;

        if(rhs + EPS >= LHS)
            colStart = mid + 1;
        else
            colEnd = mid - 1;
    }

    return colEnd + 1;
}

int binary_search_complement(const vector< pair<int, double> >& columns, double partialLHS, double rhs, int colStart, int colEnd)
{
    int mid;
    while(colStart <= colEnd)
    {
        mid = (colStart + colEnd) / 2;
        double LHS = partialLHS - min(0.0, columns[mid].second);

        if(rhs + EPS >= LHS)
            colEnd = mid - 1;
        else
            colStart = mid + 1;
    }

    return colStart - 1;
}

void cliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs, vector< pair< int, int > > &cvec)
{
    int nElements = (int)columns.size(), cliqueStart = -1;
    double maxLHS; //maxLHS = lower bound for LHS when the two variables with highest coefficients are activated.
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
    for(int i = cliqueStart, j = 0; i < nElements; i++)
    {
        idxs[j++] = columns[i].first;
        cliqueSize++;
    }
    //process the first clique found
    processClique( cliqueSize, (const int *)idxs, cgraph, cvec );

    //now we have to check the variables that are outside of the clique found.
    for(int i = cliqueStart - 1; i >= 0; i--)
    {
        int idx = columns[i].first;
        double coef = columns[i].second;
        double partialLHS = sumNegCoefs - min(0.0, coef) + coef;

        maxLHS = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[nElements-1].second)
                 + columns[i].second + columns[nElements-1].second;

        if(maxLHS <= rhs + EPS) return;

        int position = binary_search(columns, partialLHS, rhs, cliqueStart, nElements - 1);

        assert(position >= 0 && position < nElements);

        int n = nElements - position + 1, idxs[n];
        cliqueSize = 1;
        idxs[0] = idx;
        for(int i = position, j = 1; i < nElements; i++)
        {
            idxs[j++] = columns[i].first;
            cliqueSize++;
        }
        processClique( cliqueSize, (const int *)idxs, cgraph, cvec );
    }
}

void cliqueComplementDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs, vector< pair< int, int > > &cvec)
{
    int nElements = (int)columns.size(), cliqueCompStart = -1;
    int nCols = cgraph_size(cgraph) / 2;
    double maxLHS; //maxLHS = lower bound for LHS when the two variaveis with smallest coefficients are deactivated.
    int cliqueCompSize = 0;

    maxLHS = sumNegCoefs - min(0.0, columns[1].second) - min(0.0, columns[0].second);

    if(maxLHS <= rhs + EPS) return; //there is no clique involving the complement of variables in this constraint.

    for(int i = nElements - 1; i > 0; i--)
    {
        double D = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[i-1].second);

        if(D > rhs + EPS)
        {
            cliqueCompStart = i;
            break;
        }
    }

    assert(cliqueCompStart > 0 && cliqueCompStart < nElements);
    int n = cliqueCompStart + 1, idxs[n];
    for(int i = 0; i < n; i++)
    {
        idxs[i] = columns[i].first + nCols; //binary complement
        cliqueCompSize++;
    }

    //process the first clique found
    processClique( cliqueCompSize, (const int *)idxs, cgraph, cvec );

    //now we have to check the variables that are outside of the clique found.
    for(int i = cliqueCompStart + 1; i < nElements; i++)
    {
        int idx = columns[i].first;
        double coef = columns[i].second;
        double partialLHS = sumNegCoefs - min(0.0, coef);

        maxLHS = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[0].second);

        if(maxLHS <= rhs + EPS) return;

        int position = binary_search_complement(columns, partialLHS, rhs, 0, cliqueCompStart);

        assert(position >=0 && position < nElements);

        int n = position + 2, idxs[n];
        cliqueCompSize = 1;
        idxs[0] = idx + nCols;
        for(int i = 0, j = 1; i <= position; i++)
        {
            idxs[j++] = columns[i].first + nCols;
            cliqueCompSize++;
        }
        processClique( cliqueCompSize, (const int *)idxs, cgraph, cvec );
    }
}

void mixedCliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs, vector< pair< int, int > > &cvec)
{
    int nElements = (int)columns.size();
    int nCols = cgraph_size(cgraph) / 2;

    //Looking for conflicts like (x = 0, y = 1)
    for(int i = 0; i < nElements - 1; i++)
    {
        int idx = columns[i].first;
        double coef = columns[i].second;
        double maxLHS = sumNegCoefs - min(0.0, coef) - min(0.0, columns[nElements-1].second) + columns[nElements-1].second;

        if(maxLHS <= rhs + EPS) break; //there are no conflicts here

        double partialLHS = sumNegCoefs - min(0.0, coef);
        int position = binary_search(columns, partialLHS, rhs, i+1, nElements - 1);

        assert(position >= i+1 && position < nElements);

        //Conflicts between the complement of variable with index idx and
        //all the variables with indexes in the range [position, nElements-1]
        for(int j = position; j < nElements; j++)
            cvec.push_back( pair<int,int>(idx + nCols, columns[j].first) );
        fetchConflicts(false, cgraph, cvec);
    }

    //Looking for conflicts like (x = 1, y = 0)
    for(int i = nElements - 2; i >= 0; i--)
    {
        int idx = columns[i].first;
        double coef = columns[i].second;

        double maxLHS = sumNegCoefs - min(0.0, coef) - min(0.0, columns[i+1].second) + coef;

        if(maxLHS <= rhs + EPS) break; //there are no conflicts here

        double partialLHS = sumNegCoefs - min(0.0, coef) + coef;
        int position = binary_search_complement(columns, partialLHS, rhs, i+1, nElements - 1);

        assert(position >= i+1 && position < nElements);

        //Conflicts between the variable with index idx and
        //all the complements of variables with indexes in the range [i+1, position]
        for(int j = i+1; j <= position; j++)
            cvec.push_back( pair<int,int>(idx, columns[j].first + nCols) );
        fetchConflicts(false, cgraph, cvec);
    }
}

CGraph *build_cgraph_osi( const void *_solver )
{
    recomputeDegree = 0;
    startCG = CoinCpuTime();
    vector< pair< int, int > > cvec;
    cvec.reserve( CVEC_CAP );

    const OsiSolverInterface *solver = (const OsiSolverInterface *) _solver;

    if (solver->getNumIntegers()<2)
        return NULL;

    int cgraphSize = solver->getNumCols() * 2;
    CGraph *cgraph = cgraph_create( cgraphSize );
    const char *ctype = solver->getColType();
    const CoinPackedMatrix *M = solver->getMatrixByRow();
    const double *rhs = solver->getRightHandSide();
    const char *sense = solver->getRowSense();
    int nCols = solver->getNumCols();
    int nRows = solver->getNumRows();
    int idxRow;

    for(int i = 0; i < nCols; i++)
        /* inserting trivial conflicts: variable-complement */
        if(ctype[i] == 1) //consider only binary variables
            cvec.push_back( pair<int, int>(i, i + nCols) );

    for(idxRow = 0; idxRow < nRows; idxRow++)
    {

        if(CoinCpuTime() - startCG >= MAX_TIME_CG)
            break;

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
            printf("TODO: CHECK FOR RANGED CONSTRAINT (%s) rhs is %g\n", solver->getRowName(idxRow).c_str(), rhs[idxRow] );
            continue;
        }

        double mult = (sense[idxRow] == 'G') ? -1.0 : 1.0;
        for(int i = 0; i < nElements; i++)
        {
            const int cidx = idx[i];

            columns[i].first = idx[i];
            columns[i].second = coefs[i] * mult;

            if(ctype[cidx] == 1)
                nBools++;

            if(columns[i].second <= -EPS)
                sumNegCoefs += columns[i].second;
            else nPos++;

            minCoef = min(minCoef, columns[i].second);
            maxCoef = max(maxCoef, columns[i].second);
        }

        if(nBools < nElements)
            continue;

        /* special case: GUB constraints */
        if ( DBL_EQUAL( minCoef, maxCoef ) &&  DBL_EQUAL( maxCoef, rhs[idxRow] * mult ) &&
                DBL_EQUAL(minCoef, 1.0) && ((sense[idxRow]=='E') || (sense[idxRow]=='L'))
                && (row.getNumElements() > 3) )
        {
            processClique( row.getNumElements(), (const int *)idx, cgraph, cvec );
        }

        else
        {
            sort(columns.begin(), columns.end(), sort_columns());
            cliqueDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult, cvec);
            cliqueComplementDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult, cvec);
            mixedCliqueDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult, cvec);

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

                cliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow], cvec);
                cliqueComplementDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow], cvec);
                mixedCliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow], cvec);
            }
        }
    }

    fetchConflicts(true, cgraph, cvec);

    if(recomputeDegree)
        cgraph_recompute_degree( cgraph );
    else
        cgraph_update_min_max_degree( cgraph );

    return cgraph;
}

CGraph *build_cgraph_lp( const void *_mip )
{
    recomputeDegree = 0;
    startCG = CoinCpuTime();
    vector< pair< int, int > > cvec;
    cvec.reserve( CVEC_CAP );

    const LinearProgram *mip = (const LinearProgram *) _mip;

    if(lp_num_binary_cols(mip) < 2)
        return NULL;

    const int nCols = lp_cols(mip);
    const int nRows = lp_rows(mip);
    const int cgraphSize = nCols * 2;
    CGraph *cgraph = cgraph_create( cgraphSize );
    int *idxs = new int[nCols];
    double *coefs = new double[nCols];

    for(int i = 0; i < nCols; i++)
        /* inserting trivial conflicts: variable-complement */
        if(lp_is_binary(mip, i))//consider only binary variables
            cvec.push_back( pair<int, int>(i, i + nCols) );

    for(int idxRow = 0; idxRow < nRows; idxRow++)
    {
        const int nElements = lp_row(mip, idxRow, idxs, coefs);
        const double rhs = lp_rhs(mip, idxRow);
        const char sense = lp_sense(mip, idxRow);
        vector< pair<int, double> > columns(nElements);
        int nBools = 0; // number of binary variables
        int nPos = 0; //number of positive coefficients
        double sumNegCoefs = 0.0; //sum of all negative coefficients
        double minCoef = numeric_limits<double>::max();
        double maxCoef = numeric_limits<double>::min();

        if ( (nElements<2) || (fabs(rhs)>=LARGE_CONST) )
            continue;

        if ( sense == 'R' )  // lets not consider ranged constraints by now
        {
            char name[256];
            printf("TODO: CHECK FOR RANGED CONSTRAINT (%s) rhs is %g\n", lp_row_name(mip, idxRow, name), rhs );
            continue;
        }

        double mult = (sense == 'G') ? -1.0 : 1.0;
        for(int i = 0; i < nElements; i++)
        {
            columns[i].first = idxs[i];
            columns[i].second = coefs[i] * mult;

            if(lp_is_binary(mip, columns[i].first))
                nBools++;

            if(columns[i].second <= -EPS)
                sumNegCoefs += columns[i].second;
            else nPos++;

            minCoef = min(minCoef, columns[i].second);
            maxCoef = max(maxCoef, columns[i].second);
        }

        if(nBools < nElements)
            continue;

        /* special case: GUB constraints */
        if ( DBL_EQUAL( minCoef, maxCoef ) &&  DBL_EQUAL( maxCoef, rhs * mult ) &&
                DBL_EQUAL(minCoef, 1.0) && ((sense=='E') || (sense=='L'))
                && (nElements > 3) )
        {
            processClique( nElements, idxs, cgraph, cvec );
        }

        else
        {
            sort(columns.begin(), columns.end(), sort_columns());
            cliqueDetection(cgraph, columns, sumNegCoefs, rhs * mult, cvec);
            cliqueComplementDetection(cgraph, columns, sumNegCoefs, rhs * mult, cvec);
            mixedCliqueDetection(cgraph, columns, sumNegCoefs, rhs * mult, cvec);

            /*equality constraints are converted into two inequality constraints (<=).
            the first one is analyzed above and the second (multiplying the constraint by -1) is analyzed below.
            Example: x + y + z = 2 ==>  (x + y + z <= 2) and (- x - y - z <= -2)*/
            if(sense == 'E')
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

                cliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs, cvec);
                cliqueComplementDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs, cvec);
                mixedCliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs, cvec);
            }
        }
    }

    fetchConflicts(true, cgraph, cvec);

    if(recomputeDegree)
        cgraph_recompute_degree( cgraph );
    else
        cgraph_update_min_max_degree( cgraph );

    delete[] idxs;
    delete[] coefs;

    return cgraph;
}
