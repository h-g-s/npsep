#include "cut.h"
#include "memory.h"
#include "macros.h"
#include <string.h>
#include <assert.h>
#include <float.h>

#define INITIAL_CAPACITY 500
#define INCREASE_CAPACITY 500

struct _Cut {
    int n, numActiveCols;
    double rhs, violation;
    int *idx;
    double *coef;

    double fitness;
};

struct _CutPool {
    Cut **cuts;
    int n, capacity;

    int numCols;
    int *bestCutByCol;
    int *cutFrequency;
};

Cut* cut_create( const int *idxs, const double *coefs, int nz, double rhs, const double *x)
{
    Cut *cut;

    ALLOCATE(cut, Cut);
    ALLOCATE_VECTOR(cut->idx, int, nz);
    ALLOCATE_VECTOR(cut->coef, double, nz);

    cut->n = nz;
    cut->rhs = rhs;
    cut->numActiveCols = 0;

    double lhs = 0.0;
    double minCoef = DBL_MAX, maxCoef = -DBL_MAX;
    for(int i = 0; i < cut->n; i++)
    {
        cut->idx[i] = idxs[i];
        cut->coef[i] = coefs[i];

        if(fabs(x[cut->idx[i]]) > 1e-6)
        {
            lhs += (cut->coef[i] * x[cut->idx[i]]);
            cut->numActiveCols++;
            minCoef = MIN(minCoef, cut->coef[i]);
            maxCoef = MAX(maxCoef, cut->coef[i]);
        }
    }

    assert(cut->numActiveCols > 0);
    const double diffOfCoefs = fabs(maxCoef - minCoef) + fabs(maxCoef - cut->rhs) + fabs(minCoef - cut->rhs);
    cut->violation = lhs - rhs;
    cut->fitness = ( (cut->violation / ((double)cut->numActiveCols)) * 100000.0 ) +
                     ( (1.0 / (diffOfCoefs + 1.0)) * 100.0 );

    assert(cut->violation > 1e-6);

    return cut;
}

void cut_free( Cut **_cut )
{
    Cut *cut = *_cut;
    free(cut->idx);
    free(cut->coef);
    free(cut);
    *_cut = NULL;
}

int cut_size( const Cut *cut ) { return cut->n; }
const int* cut_get_idxs( const Cut *cut ) { return cut->idx; }
const double* cut_get_coefs( const Cut *cut ) { return cut->coef; }
double cut_get_rhs( const Cut *cut ) { return cut->rhs; }
double cut_get_violation( const Cut *cut ) { return cut->violation; }
int cut_get_num_active_cols( const Cut *cut ) { return cut->numActiveCols; }
double cut_get_fitness( const Cut *cut ) { return cut->fitness; }

int bin_search( const int *v, const int n, const int x )
{
    int mid, left = 0, right = n-1;

    while(left <= right)
    {
        mid = (left + right) / 2;

        if(v[mid] == x) return mid;
        else if(v[mid] < x) left = mid + 1;
        else right = mid - 1;
    }
    return -1;
}

int cut_check_domination(  const Cut *cutA, const Cut *cutB )
{
    if(cutA->violation <= cutB->violation - EPS)
        return 0;

    /* rhsA == 0 && rhsB < 0 */
    if(fabsl(cutA->rhs) < EPS && cutB->rhs <= -EPS)
        return 0;

    /* rhsA > 0 && rhsB == 0 */
    if(cutA->rhs > EPS && fabsl(cutB->rhs) < EPS)
        return 0;

    int sizeA = cutA->n, sizeB = cutB->n;
    const int *idxsA = cutA->idx, *idxsB = cutB->idx;
    const double *coefsA = cutA->coef, *coefsB = cutB->coef;
    double normConstA, normConstB;
    char analyzed[sizeB];

    if(fabsl(cutA->rhs) < EPS || fabsl(cutB->rhs) < EPS)
        normConstA = normConstB = 1.0L;
    else
    {
        normConstA = cutA->rhs;
        normConstB = cutB->rhs;
    }

    memset(analyzed, 0, sizeof(char) * sizeB);

    int i;
    for(i = 0; i < sizeA; i++)
    {
        int idxA = idxsA[i];
        double normCoefA = (coefsA[i] / normConstA);
        int posB = bin_search(idxsB, sizeB, idxA);

        if(posB == -1 && normCoefA <= -EPS)
            return 0;
        else if(posB != -1)
        {
            double normCoefB = (coefsB[posB] / normConstB);
            analyzed[posB] = 1;
            if(normCoefA < normCoefB - EPS)
                return 0;
        }
    }

    for(i = 0; i < sizeB; i++)
    {
        if(analyzed[i]) continue;

        int idxB = idxsB[i];
        double normCoefB = (coefsB[i] / normConstB);
        int posA = bin_search(idxsA, sizeA, idxB);

        if(posA == -1 && normCoefB > EPS)
            return 0;
    }

    return 1;
}

int cut_is_equal( const Cut *cutA, const Cut *cutB )
{
    if(cutA->violation != cutB->violation || cutA->n != cutB->n)
        return 0;

    /* rhsA == 0 && rhsB != 0 */
    if(fabsl(cutA->rhs) < EPS && fabsl(cutB->rhs) > EPS)
        return 0;

    /* rhsA != 0 && rhsB == 0 */
    if(fabsl(cutA->rhs) > EPS && fabsl(cutB->rhs) < EPS)
        return 0;

    const int *idxsA = cutA->idx, *idxsB = cutB->idx;
    const double *coefsA = cutA->coef, *coefsB = cutB->coef;
    double normConstA, normConstB;

    if(fabsl(cutA->rhs) < EPS && fabsl(cutB->rhs) < EPS)
        normConstA = normConstB = 1.0L;
    else
    {
        normConstA = cutA->rhs;
        normConstB = cutB->rhs;
    }

    int i;
    for(i = 0; i < cutA->n; i++)
    {
        int idxA = idxsA[i], idxB = idxsB[i];
        double normCoefA = (coefsA[i] / normConstA), normCoefB = (coefsB[i] / normConstB);

        if(idxA != idxB) return 0;
        if(fabsl(normCoefA - normCoefB) > EPS) return 0;
    }
    return 1;
}

int cut_domination( const Cut *cutA, const Cut *cutB )
{
    /* checks if cutA and cutB are equivalent */
    if(cut_is_equal(cutA, cutB))
        return 0;

    /* checks if cutA dominates cutB */
    if(cut_check_domination(cutA, cutB))
        return 1;

    /* checks if cutB dominates cutA */
    if(cut_check_domination(cutB, cutA))
        return 2;

    /* cutA and cutB are not dominated */
    return 3;
}

CutPool* cut_pool_create(const int numCols)
{
    CutPool *cutpool;
    ALLOCATE(cutpool, CutPool);

    cutpool->capacity = INITIAL_CAPACITY;
    cutpool->n = 0;
    cutpool->numCols = numCols;

    ALLOCATE_VECTOR(cutpool->cuts, Cut*, cutpool->capacity);
    ALLOCATE_VECTOR(cutpool->bestCutByCol, int, cutpool->numCols);
    FILL(cutpool->bestCutByCol, 0, cutpool->numCols, -1);
    ALLOCATE_VECTOR(cutpool->cutFrequency, int, cutpool->capacity);

    return cutpool;
}

void cut_pool_free( CutPool **_cutpool )
{
    CutPool *cutpool = *_cutpool;

    int i;
    for(i = 0; i < cutpool->n; i++)
        cut_free(&cutpool->cuts[i]);

    free(cutpool->cuts);
    free(cutpool->bestCutByCol);
    free(cutpool->cutFrequency);
    free(cutpool);
    *_cutpool = NULL;
}

void update_best_cut_by_col(CutPool *cutpool, const int idxCut)
{
    const Cut *cut = cutpool->cuts[idxCut];
    const int nz = cut->n;
    const int *idxs = cut->idx;
    const double fitness = cut->fitness;

    for(int i = 0; i < nz; i++)
    {
        const int idx = idxs[i];

        if(cutpool->bestCutByCol[idx] == -1)
        {
            cutpool->bestCutByCol[idx] = idxCut;
            cutpool->cutFrequency[idxCut]++;
        }
        else
        {
            const int currBestCut = cutpool->bestCutByCol[idx];
            const double currFitness = cutpool->cuts[currBestCut]->fitness;

            if(fitness > currFitness + EPS)
            {
                cutpool->cutFrequency[idxCut]++;
                cutpool->cutFrequency[currBestCut]--;
                cutpool->bestCutByCol[idx] = idxCut;
                assert(cutpool->cutFrequency[currBestCut] >= 0);
                assert(cutpool->cutFrequency[currBestCut] < cutpool->numCols);
            }
        }
    }
}

int cut_pool_insert( CutPool *cutpool, Cut *newCut )
{
    if(cutpool->n == cutpool->capacity)
    {
        cutpool->capacity = cutpool->capacity + INCREASE_CAPACITY;
        cutpool->cuts = realloc(cutpool->cuts, sizeof(Cut*) * cutpool->capacity);
        cutpool->cutFrequency = realloc(cutpool->cutFrequency, sizeof(int) * cutpool->capacity);
    }

    int position = cutpool->n;
    cutpool->cuts[position] = newCut;
    cutpool->cutFrequency[position] = 0;
    update_best_cut_by_col(cutpool, position);

    if(!cutpool->cutFrequency[position])
        return 0;

    cutpool->n++;

    return 1;
}

void cut_pool_update(CutPool *cutPool)
{
    if(cutPool->n < 2)
        return;

    char *removed;
    int nRemoved = 0;

    ALLOCATE_VECTOR_INI(removed, char, cutPool->n);

    for(int i = 0; i < cutPool->n; i++)
    {
        if(removed[i])
            continue;

        const Cut *cutA = cutPool->cuts[i];

        if(!cutPool->cutFrequency[i])
        {
            removed[i] = 1;
            nRemoved++;
            cut_free(&cutPool->cuts[i]);
            cutPool->cuts[i] = NULL;
            continue;
        }

        for(int j = i + 1; j < cutPool->n; j++)
        {
            if(removed[j] || !cutPool->cutFrequency[j])
                continue;

            const Cut *cutB = cutPool->cuts[j];
            int chkDm = cut_domination(cutA, cutB);

            if(chkDm == 0 || chkDm == 2)//cutA is equivalent to cutB or cutB dominates cutA
            {
                removed[i] = 1;
                nRemoved++;
                cut_free(&cutPool->cuts[i]);
                cutPool->cuts[i] = NULL;
                break;
            }

            else if(chkDm == 1)//cutA dominates cutB
            {
                removed[j] = 1;
                nRemoved++;
                cut_free(&cutPool->cuts[j]);
                cutPool->cuts[j] = NULL;
                continue;
            }
        }
    }

    assert(nRemoved >= 0 && nRemoved < cutPool->n);

    int last = 1;
    for(int i = 0; i < cutPool->n; i++)
    {
        if(!removed[i])
            continue;

        last = MAX(last, i+1);

        while(last < cutPool->n)
        {
            if(!removed[last])
                break;
            last++;
        }

        if(last >= cutPool->n)
            break;

        cutPool->cuts[i] = cutPool->cuts[last];
        removed[last] = 1;
        last++;
    }

    cutPool->n = cutPool->n - nRemoved;

    free(removed);
}

int cut_pool_size( const CutPool *cutpool )
{
    if(!cutpool)
        return 0;
    return cutpool->n;
}

Cut* cut_pool_get_cut( const CutPool *cutpool, int idx )
{
    assert(idx >= 0 && idx < cutpool->n);
    return cutpool->cuts[idx];
}

int cut_pool_cut_frequency( const CutPool *cutpool, int idx )
{
    assert(idx >= 0 && idx < cutpool->n);
    return cutpool->cutFrequency[idx];
}
