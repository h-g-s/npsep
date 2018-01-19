#include "knapsack_separation.hpp"
#include <float.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <map>

extern "C"
{
    #include "macros.h"
    #include "clique.h"
    #include "bron_kerbosch.h"
    #include "strutils.h"
}

#define NOT_INSERTED                    0
#define INSERTED_AT_SEPARATION          1
#define INSERTED_AT_SIMPLE_EXTENSION    2
#define INSERTED_AT_CLIQUE_EXTENSION    3

#define KSEP_DEF_MIN_VIOL               0.02
#define KSEP_DEF_MAX_RC                 DBL_MAX/10.0
#define KSEP_DEF_MAX_IT_BK              INT_MAX/10

#define KSEP_STR_MAX_RC                 "maxRC"
#define KSEP_STR_MAX_IT_BK              "maxItBK"

#define CAPACITY_LIMIT  1000000
#define BOTTOM_UP_LIMIT 1000000

extern double startTime, maxTime;

struct _KnapsackSeparation
{
    int nRows, nCols;
    int *nz;
    int **knpRowIdxs;
    double **knpRowCoefs;
    char **isComplement;
    char *rowHasCompVars;
    double *rhs;
    double maxRC;

    /*lifting*/
    int maxItBK;

    /*dynamic programming*/
    double *maxCoef, *totalWeight;

    const CGraph *cg;
};

double fractionalPart( const double x )
{
    double nextInteger = ceil( x );
    double previousInteger = floor( x );

    return MIN( nextInteger-x, x-previousInteger );
}

KnapsackSeparation *knapsack_sep_create(const OsiSolverInterface &si)
{
    KnapsackSeparation *ksep = new KnapsackSeparation;

    const int nRows = si.getNumRows(), nCols = si.getNumCols();
    const char *rowSense = si.getRowSense();
    const double *rhs = si.getRightHandSide();
    const CoinPackedMatrix *M = si.getMatrixByRow();

    int i, j, nz = 0;

    ksep->nCols = nCols;
    ksep->nRows = 0;
    ksep->cg = NULL;
    ksep->maxRC = KSEP_DEF_MAX_RC;
    ksep->maxItBK = KSEP_DEF_MAX_IT_BK;

    ksep->knpRowIdxs = new int*[nRows];
    ksep->knpRowIdxs[0] = new int[si.getNumElements()];
    ksep->knpRowCoefs = new double*[nRows];
    ksep->knpRowCoefs[0] = new double[si.getNumElements()];
    ksep->isComplement = new char*[nRows];
    ksep->isComplement[0] = new char[si.getNumElements()];
    ksep->rowHasCompVars = new char[nRows];
    ksep->nz = new int[nRows]();
    ksep->rhs = new double[nRows];
    ksep->maxCoef = new double[nRows];
    ksep->totalWeight = new double[nRows];

    for(i = 0; i < nRows; i++)
    {
        const CoinShallowPackedVector &row = M->getVector(i);
        const int nzRow = row.getNumElements();
        const int *idxRow = row.getIndices();
        const double *coefRow = row.getElements();
        const char sense = rowSense[i];
        char onlyBin = True;

        if(nzRow < 2 || sense != 'L')
            continue;

        ksep->rhs[ksep->nRows] = rhs[i];
        ksep->nz[ksep->nRows] = nzRow;
        ksep->rowHasCompVars[ksep->nRows] = 0;
        ksep->maxCoef[ksep->nRows] = fabs(coefRow[0]);
        ksep->totalWeight[ksep->nRows] = 0.0;

        double minCoef = DBL_MAX, sumCoef = 0.0, maxFrac = 0.0;

        for(j = 0; j < nzRow; j++)
        {
            if(!si.isBinary(idxRow[j]))
            {
                onlyBin = False;
                break;
            }

            ksep->knpRowIdxs[0][nz] = idxRow[j];

            if(coefRow[j] <= -EPS)
            {
                ksep->knpRowCoefs[0][nz] = fabs(coefRow[j]);
                ksep->isComplement[0][nz] = 1;
                ksep->rhs[ksep->nRows] += ksep->knpRowCoefs[0][nz];
                ksep->rowHasCompVars[ksep->nRows] = 1;
            }
            else
            {
                ksep->knpRowCoefs[0][nz] = coefRow[j];
                ksep->isComplement[0][nz] = 0;
            }

            minCoef = MIN(minCoef, ksep->knpRowCoefs[0][nz]);
            ksep->maxCoef[ksep->nRows] = MAX(ksep->maxCoef[ksep->nRows], ksep->knpRowCoefs[0][nz]);
            ksep->totalWeight[ksep->nRows] += ksep->knpRowCoefs[0][nz];
            sumCoef += ksep->knpRowCoefs[0][nz];
            maxFrac = MAX(maxFrac, fractionalPart(ksep->knpRowCoefs[0][nz]));
            nz++;
        }

        maxFrac = MAX(maxFrac, fractionalPart(ksep->rhs[ksep->nRows]));

        int nEntries, capacity;
        if(maxFrac > 1e-5)
        	capacity = ((int)(floor(ksep->maxCoef[ksep->nRows] * 1000.0)))
        		     + ((int)(ceil(ksep->rhs[ksep->nRows] * 1000.0)));
        else
        	capacity = ((int)(ksep->maxCoef[ksep->nRows])) + ((int)(ksep->rhs[ksep->nRows]));
        nEntries = (ksep->nz[ksep->nRows] + 1) * (capacity + 1);

        if( !onlyBin || nEntries > CAPACITY_LIMIT || nEntries <= 0 || ksep->rhs[ksep->nRows] <= 1.0||
        	sumCoef <= ksep->rhs[ksep->nRows] || minCoef > ksep->rhs[ksep->nRows] )
        {
            nz = nz - j;
            continue;
        }

        /* se existem coeficientes fracionários, transforma para inteiros */
        if(maxFrac > 1e-5)
        {
            ksep->totalWeight[ksep->nRows] = 0.0;
            for(j = 1; j <= ksep->nz[ksep->nRows]; j++)
            {
                ksep->knpRowCoefs[0][nz-j] = floor(ksep->knpRowCoefs[0][nz-j] * 1000.0);
                assert(ksep->knpRowCoefs[0][nz-j] > 0);
                ksep->maxCoef[ksep->nRows] = MAX(ksep->maxCoef[ksep->nRows], ksep->knpRowCoefs[0][nz-j]);
                ksep->totalWeight[ksep->nRows] += ksep->knpRowCoefs[0][nz-j];
            }
            ksep->rhs[ksep->nRows] = ceil(ksep->rhs[ksep->nRows] * 1000.0);
        }

        ksep->nRows++;
    }

    for(i = 1; i < ksep->nRows; i++)
    {
        ksep->knpRowIdxs[i] = ksep->knpRowIdxs[i-1] + ksep->nz[i-1];
        ksep->knpRowCoefs[i] = ksep->knpRowCoefs[i-1] + ksep->nz[i-1];
        ksep->isComplement[i] = ksep->isComplement[i-1] + ksep->nz[i-1];
    }

    return ksep;
}

KnapsackSeparation* knapsack_sep_create(const LinearProgram *lp)
{
    KnapsackSeparation *ksep = new KnapsackSeparation;

    const int nRows = lp_rows(lp), nCols = lp_cols(lp);
    int i, j, nz = 0;
    int *rowIdxs = new int[nCols];
    double *rowCoefs = new double[nCols];

    ksep->nCols = nCols;
    ksep->nRows = 0;
    ksep->cg = NULL;
    ksep->maxRC = KSEP_DEF_MAX_RC;
    ksep->maxItBK = KSEP_DEF_MAX_IT_BK;

    ksep->knpRowIdxs = new int*[nRows];
    ksep->knpRowIdxs[0] = new int[lp_nz(lp)];
    ksep->knpRowCoefs = new double*[nRows];
    ksep->knpRowCoefs[0] = new double[lp_nz(lp)];
    ksep->isComplement = new char*[nRows];
    ksep->isComplement[0] = new char[lp_nz(lp)];
    ksep->rowHasCompVars = new char[nRows];
    ksep->nz = new int[nRows]();
    ksep->rhs = new double[nRows];
    ksep->maxCoef = new double[nRows];
    ksep->totalWeight = new double[nRows];

    for(i = 0; i < nRows; i++)
    {
        const int nzRow = lp_row(lp, i, rowIdxs, rowCoefs);
        const char sense = lp_sense(lp, i);
        char onlyBin = True;

        if(nzRow < 2 || sense != 'L')
            continue;

        ksep->rhs[ksep->nRows] = lp_rhs(lp, i);
        ksep->nz[ksep->nRows] = nzRow;
        ksep->rowHasCompVars[ksep->nRows] = 0;
        ksep->maxCoef[ksep->nRows] = fabs(rowCoefs[0]);
        ksep->totalWeight[ksep->nRows] = 0.0;

        double minCoef = DBL_MAX, sumCoef = 0.0, maxFrac = 0.0;

        for(j = 0; j < nzRow; j++)
        {
            if(!lp_is_binary(lp, rowIdxs[j]))
            {
                onlyBin = False;
                break;
            }

            ksep->knpRowIdxs[0][nz] = rowIdxs[j];

            if(rowCoefs[j] <= -EPS)
            {
                ksep->knpRowCoefs[0][nz] = fabs(rowCoefs[j]);
                ksep->isComplement[0][nz] = 1;
                ksep->rhs[ksep->nRows] += ksep->knpRowCoefs[0][nz];
                ksep->rowHasCompVars[ksep->nRows] = 1;
            }
            else
            {
                ksep->knpRowCoefs[0][nz] = rowCoefs[j];
                ksep->isComplement[0][nz] = 0;
            }

            minCoef = MIN(minCoef, ksep->knpRowCoefs[0][nz]);
            ksep->maxCoef[ksep->nRows] = MAX(ksep->maxCoef[ksep->nRows], ksep->knpRowCoefs[0][nz]);
            ksep->totalWeight[ksep->nRows] += ksep->knpRowCoefs[0][nz];
            sumCoef += ksep->knpRowCoefs[0][nz];
            maxFrac = MAX(maxFrac, fractionalPart(ksep->knpRowCoefs[0][nz]));
            nz++;
        }

        maxFrac = MAX(maxFrac, fractionalPart(ksep->rhs[ksep->nRows]));

        int nEntries, capacity;
        if(maxFrac > 1e-5)
        	capacity = ((int)(floor(ksep->maxCoef[ksep->nRows] * 1000.0)))
        		     + ((int)(ceil(ksep->rhs[ksep->nRows] * 1000.0)));
        else
        	capacity = ((int)(ksep->maxCoef[ksep->nRows])) + ((int)(ksep->rhs[ksep->nRows]));
        nEntries = (ksep->nz[ksep->nRows] + 1) * (capacity + 1);

        if( !onlyBin || nEntries > CAPACITY_LIMIT || nEntries <= 0 || ksep->rhs[ksep->nRows] <= 1.0||
        	sumCoef <= ksep->rhs[ksep->nRows] || minCoef > ksep->rhs[ksep->nRows] )
        {
            nz = nz - j;
            continue;
        }

        /* se existem coeficientes fracionários, transforma para inteiros */
        if(maxFrac > 1e-5)
        {
            ksep->totalWeight[ksep->nRows] = 0.0;
            for(j = 1; j <= ksep->nz[ksep->nRows]; j++)
            {
                ksep->knpRowCoefs[0][nz-j] = floor(ksep->knpRowCoefs[0][nz-j] * 1000.0);
                assert(ksep->knpRowCoefs[0][nz-j] > 0);
                ksep->maxCoef[ksep->nRows] = MAX(ksep->maxCoef[ksep->nRows], ksep->knpRowCoefs[0][nz-j]);
                ksep->totalWeight[ksep->nRows] += ksep->knpRowCoefs[0][nz-j];
            }
            ksep->rhs[ksep->nRows] = ceil(ksep->rhs[ksep->nRows] * 1000.0);
        }

        ksep->nRows++;
    }

    for(i = 1; i < ksep->nRows; i++)
    {
        ksep->knpRowIdxs[i] = ksep->knpRowIdxs[i-1] + ksep->nz[i-1];
        ksep->knpRowCoefs[i] = ksep->knpRowCoefs[i-1] + ksep->nz[i-1];
        ksep->isComplement[i] = ksep->isComplement[i-1] + ksep->nz[i-1];
    }

    delete[] rowIdxs;
    delete[] rowCoefs;

    return ksep;
}

void knapsack_sep_free(KnapsackSeparation **_ksep)
{
    KnapsackSeparation *ksep = *_ksep;

    delete[] ksep->knpRowIdxs[0];
    delete[] ksep->knpRowIdxs;
    delete[] ksep->knpRowCoefs[0];
    delete[] ksep->knpRowCoefs;
    delete[] ksep->isComplement[0];
    delete[] ksep->isComplement;
    delete[] ksep->rowHasCompVars;
    delete[] ksep->nz;
    delete[] ksep->rhs;
    delete[] ksep->maxCoef;
    delete[] ksep->totalWeight;
    delete ksep;
    _ksep = NULL;
}

//LinearProgram* generate_ip_program(KnapsackSeparation *ksep, const double *x, const int rowIdx)
//{
//    assert(rowIdx >= 0 && rowIdx < ksep->nRows);
//
//    const int numCols = ksep->nz[rowIdx];
//    LinearProgram *knp = lp_create();
//    int *idxs = new int[numCols];
//    double *coefs = new double[numCols];
//
//    for(int i = 0; i < numCols; i++)
//    {
//        char name[256];
//        const int idxCol = ksep->knpRowIdxs[rowIdx][i];
//
//        sprintf(name, "item_%d", i);
//
//        if(idxCol > EPS)
//        {
//            const double objValue = (ksep->isComplement[rowIdx][i]) ? x[idxCol] : (1.0 - x[idxCol]);
//            lp_add_col(knp, objValue, 0.0, 1.0, True, name, 0, NULL, NULL);
//        }
//        else
//            lp_add_col(knp, 0.0, 0.0, 0.0, True, name, 0, NULL, NULL);
//
//        idxs[i] = i;
//        coefs[i] = ksep->knpRowCoefs[rowIdx][i];
//    }
//
//    lp_add_row(knp, numCols, idxs, coefs, "capacity", 'G', ksep->rhs[rowIdx] + 1.0);
//
//    delete[] idxs;
//    delete[] coefs;
//
//    lp_set_print_messages(knp, 0);
//    lp_set_direction(knp, LP_MIN);
//    lp_set_max_seconds(knp, maxTime-ELAPSED_TIME(startTime));
//
//    return knp;
//}

int simple_lifting(KnapsackSeparation *ksep, const double *x, const double *rc, const int rowIdx, char *selected)
{
    int newVars = 0;
    double maxWeight = 0;

    for(int i = 0; i < ksep->nz[rowIdx]; i++)
        if(selected[i] == INSERTED_AT_SEPARATION || selected[i] == INSERTED_AT_CLIQUE_EXTENSION)
            maxWeight = MAX(maxWeight, ksep->knpRowCoefs[rowIdx][i]);

    for(int i = 0; i < ksep->nz[rowIdx]; i++)
    {
        const int idxCol = ksep->knpRowIdxs[rowIdx][i];
        if( (selected[i] == NOT_INSERTED) && (ksep->knpRowCoefs[rowIdx][i] >= maxWeight)
            && ((x[idxCol] > EPS) || (rc[idxCol] <= ksep->maxRC)) )
        {
            selected[i] = INSERTED_AT_SIMPLE_EXTENSION;
            newVars++;
        }
    }

    return newVars;
}

int clique_lifting(KnapsackSeparation *ksep, const double *x, const double *rc, const int rowIdx, char *selected)
{
    if(!ksep->cg)
    {
        fprintf(stderr, "Cannot perform lifting using cliques. CGraph is not filled\n");
        exit(EXIT_FAILURE);
    }

    const int cgSize = cgraph_size(ksep->cg);
    const int nzRow = ksep->nz[rowIdx];
    int *base = new int[nzRow];
    int *candidates = new int[nzRow];
    int *nindexes = new int[cgSize];
    int *nodeWeight = new int[cgSize];
    int nBase = 0, nCandidates = 0, newVars = 0;
    double *minWeight = new double[nzRow];
    double *maxWeight = new double[nzRow];
    double sumMinWeight = 0.0, sumMaxWeight = 0.0, minMax = DBL_MAX;
    double minRC = DBL_MAX, maxRC = -DBL_MAX;

    for(int i = 0; i < nzRow; i++)
    {
        if(selected[i] == INSERTED_AT_SEPARATION)
        {
            const double coef = ksep->knpRowCoefs[rowIdx][i];
            base[nBase] = i;
            minWeight[nBase] = maxWeight[nBase] = coef;
            sumMinWeight += coef;
            sumMaxWeight += coef;
            minMax = MIN(minMax, coef);
            nBase++;
        }
        else if(selected[i] == NOT_INSERTED)
        {
            const int idxCol = ksep->knpRowIdxs[rowIdx][i];

            if( (x[idxCol] > EPS) || (rc[idxCol] <= ksep->maxRC) )
            {
                candidates[nCandidates++] = i;

                if(x[idxCol] < EPS)
                {
                    minRC = MIN(minRC, rc[idxCol]);
                    maxRC = MAX(maxRC, rc[idxCol]);
                }
            }
        }
        else
        {
            fprintf(stderr, "Variable has already been selected!\n");
            exit(EXIT_FAILURE);
        }
    }

    if(nCandidates == 0)
    {
        delete[] base;
        delete[] candidates;
        delete[] minWeight;
        delete[] maxWeight;
        delete[] nindexes;
        delete[] nodeWeight;
        return 0;
    }

    for(int i = 0; i < nBase; i++)
    {
        const int baseIdx = ksep->knpRowIdxs[rowIdx][base[i]];
        int numVertices = 0;
        int *posIdx = new int[cgSize];
        FILL(posIdx, 0, cgSize, -1);
        FILL(nindexes, 0, cgSize, -1);

        for(int j = 0; j < nCandidates; j++)
        {
            const int candIdx = ksep->knpRowIdxs[rowIdx][candidates[j]];
            const double candCoef = ksep->knpRowCoefs[rowIdx][candidates[j]];

            //se a var já foi adicionada, desconsiderar nova inclusão
            if(selected[candidates[j]] != NOT_INSERTED)
                continue;

            //se a var candidata nao possui conflito com a base
            if(!cgraph_conflicting_nodes(ksep->cg, baseIdx, candIdx))
                continue;

            const double newLocalMinWeight = MIN(minWeight[i], candCoef);
            const double newSumMinWeight = sumMinWeight - minWeight[i] + newLocalMinWeight;

            double newMinMax = minMax;
            double newSumMaxWeight = sumMaxWeight;

            if(candCoef > maxWeight[i] + EPS)
            {
                newSumMaxWeight = newSumMaxWeight - maxWeight[i] + candCoef;
                newMinMax = candCoef;
                for(int k = 0; k < nBase; k++)
                    if(k != i)
                        newMinMax = MIN(newMinMax, maxWeight[k]);
            }

            //verificando se a inclusao da nova var elimina a propriedade de cobertura mínima
            if(newSumMinWeight <= ksep->rhs[rowIdx] || newSumMaxWeight - newMinMax > ksep->rhs[rowIdx] + EPS)
                continue;

            nindexes[candIdx] = numVertices;
            posIdx[candIdx] = candidates[j];

            if(x[candIdx] > EPS)
                nodeWeight[numVertices] = ((int)ceil(x[candIdx] * 100000.0));
            else
            {
                if(fabs(maxRC-minRC) < EPS)
                    nodeWeight[numVertices] = 10;
                else
                {
                    const double norm = ((rc[candIdx] - minRC) / (maxRC - minRC)) + 1.0; /*range: [1,2]*/
                    nodeWeight[numVertices] = ((int)nearbyint(norm*10));
                }
            }
            numVertices++;
        }

        if(numVertices == 0)
        {
            delete[] posIdx;
            continue;
        }

        //construindo grafo e rodando bk
        {
            CGraph *cg = cgraph_create_induced_subgraph(ksep->cg, nindexes);
            for(int j = 0; j < numVertices; j++)
                cgraph_set_node_weight(cg, j, nodeWeight[j]);

            BronKerbosch *bk = bk_create(cg);
            bk_set_max_it(bk, ksep->maxItBK);
            bk_set_min_weight(bk, 0);
            bk_run(bk);
            const CliqueSet *clqSet = bk_get_clq_set(bk);
            if(clq_set_number_of_cliques(clqSet) == 0)
            {
                delete[] posIdx;
                continue;
            }

            int idxMaxW = 0, maxW = clq_set_weight(clqSet, 0);
            for(int j = 1; j < clq_set_number_of_cliques(clqSet); j++)
                if(clq_set_weight(clqSet, j) > maxW)
                {
                    maxW = clq_set_weight(clqSet, j);
                    idxMaxW = j;
                }

            const int *clqEl = clq_set_clique_elements(clqSet, idxMaxW);
            const int clqSize = clq_set_clique_size(clqSet, idxMaxW);
            for(int j = 0; j < clqSize; j++)
            {
                const int origIdx = cgraph_get_original_node_index(cg, clqEl[j]);
                const int pos = posIdx[origIdx];
                assert(pos >= 0 && pos < ksep->nz[rowIdx]);
                const double coef = ksep->knpRowCoefs[rowIdx][pos];
                selected[pos] = INSERTED_AT_CLIQUE_EXTENSION;
                newVars++;

                /*atualizando o minimo*/
                if(coef + EPS < minWeight[i])
                {
                    sumMinWeight = sumMinWeight - minWeight[i] + coef;
                    minWeight[i] = coef;
                }
                /*atualizando o maximo*/
                if(coef > maxWeight[i] + EPS)
                {
                    sumMaxWeight = sumMaxWeight - maxWeight[i] + coef;
                    maxWeight[i] = coef;
                    minMax = DBL_MAX;
                    for(int k = 0; k < nBase; k++)
                        minMax = MIN(minMax, maxWeight[k]);
                }
            }
            bk_free(bk);
            cgraph_free(&cg);
            delete[] posIdx;
        }
    }

    delete[] base;
    delete[] candidates;
    delete[] minWeight;
    delete[] maxWeight;
    delete[] nindexes;
    delete[] nodeWeight;

    return newVars;
}


double knapsack(KnapsackSeparation *ksep, const int idxRow, double **table, const double *x)
{
    int i, j;
    const int capacity = ((int)ksep->rhs[idxRow]) + ((int)ksep->maxCoef[idxRow]);
    const int nz = ksep->nz[idxRow];

    for(i = 1; i <= nz; i++)
    {
        const int idxCol = ksep->knpRowIdxs[idxRow][i-1];
        assert(idxCol >= 0 && idxCol < ksep->nCols);
        const int weight = ((int)ksep->knpRowCoefs[idxRow][i-1]);
        double profit;

        if(ksep->isComplement[idxRow][i-1])
            profit = ((1.0 - x[idxCol]) * 1000.0);
        else
            profit = (x[idxCol] * 1000.0);

        for(j = 1; j <= capacity; j++)
        {
            if(weight <= j)
                table[i][j] = MAX(profit + table[i-1][j-weight], table[i-1][j]);
            else
                table[i][j] = table[i-1][j];
        }
    }

    return table[nz][capacity];
}

double recursiveKnapsack(KnapsackSeparation *ksep, int idxRow, std::map<std::pair<int, int>, double>&table,
                         const double *x, int item, int currentCap)
{
    std::map<std::pair<int, int>, double>::iterator it = table.find(std::pair<int, int>(item, currentCap));

    if(it != table.end())
        return it->second;

    if(item == 0 || currentCap == 0)
        return 0.0;

    double value = 0.0;
    const int itemIdx = ksep->knpRowIdxs[idxRow][item-1];
    const int itemWeight = ((int)ksep->knpRowCoefs[idxRow][item-1]);
    double profit = (ksep->isComplement[idxRow][item-1]) ? ((1.0 - x[itemIdx]) * 1000.0)
                                                         : (x[itemIdx] * 1000.0);

    if(itemWeight <= currentCap)
        value = MAX(profit + recursiveKnapsack(ksep, idxRow, table, x, item-1, currentCap-itemWeight),
                    recursiveKnapsack(ksep, idxRow, table, x, item-1, currentCap));
    else
        value = recursiveKnapsack(ksep, idxRow, table, x, item-1, currentCap);

    const std::pair<int, int> newEntry(item, currentCap);
    std::pair<std::map<std::pair<int, int>, double>::iterator, bool> ret;
    ret = table.insert(std::pair<std::pair<int, int>, double>(newEntry, value));
    assert(ret.second);

    return value;
}

void separate_by_row(KnapsackSeparation *ksep, const int idxRow, const double *x, const double *rc, CutPool *cutPool)
{
    const int nz = ksep->nz[idxRow];
    const int rowRHS = ((int)ksep->rhs[idxRow]);
    const int capacity = rowRHS + ((int)ksep->maxCoef[idxRow]);
    const char bottomUp = ((nz+1) * (capacity+1)) <= BOTTOM_UP_LIMIT;

    if(ksep->totalWeight[idxRow] <= ksep->rhs[idxRow])
        return;

    int minWeight, sumWeight;
    double lhs, rhs, prevObj = 0.0;
    int *idx = new int[nz];
    double *coef = new double[nz];
    double **table = NULL;
    std::map<std::pair<int, int>, double> sparseTable;
    std::map<std::pair<int, int>, double>::iterator it;

    if(bottomUp)
    {
    	table = new double*[nz+1];
	    for(int i = 0; i <= nz; i++)
	    {
	        table[i] = new double[capacity+1];
	        for(int j = 0; j <= capacity; j++)
	            table[i][j] = 0.0;
	    }
	    knapsack(ksep, idxRow, table, x);
    }
    else
    	recursiveKnapsack(ksep, idxRow, sparseTable, x, nz, capacity);

    for(int cap = capacity; cap > rowRHS; cap--)
    {
        int currCap = cap;
        lhs = rhs = 0.0;
        minWeight = INT_MAX;
        sumWeight = 0;
        double obj;

        if(bottomUp)
        	obj = table[nz][currCap];
        else
        {
        	it = sparseTable.find(std::pair<int, int>(nz, currCap));
        	if(it == sparseTable.end())
        		obj = prevObj;
        	else
        		obj = it->second;
        }

        if(fabs(obj - prevObj) <= 1e-5)
            continue;

        prevObj = obj;
        char *selected = new char[nz]();

        for(int i = nz; i > 0; i--)
        {
        	double obj1, obj2;

        	if(bottomUp)
        	{
        		obj1 = table[i][currCap];
        		obj2 = table[i-1][currCap];
        	}
        	else
        	{
        		it = sparseTable.find(std::pair<int, int>(i, currCap));
        		assert(it != sparseTable.end());
        		obj1 = it->second;
        		it = sparseTable.find(std::pair<int, int>(i-1, currCap));
        		assert(it != sparseTable.end());
        		obj2 = it->second;
        	}

            if(fabs(obj1 - obj2) > 1e-5)
            {
                const int idxCol = ksep->knpRowIdxs[idxRow][i-1];
                const int weight = ((int)ksep->knpRowCoefs[idxRow][i-1]);

                currCap -= weight;
                assert(currCap >= 0);
                minWeight = MIN(minWeight, weight);
                sumWeight += weight;

                selected[i-1] = INSERTED_AT_SEPARATION;

                if(ksep->isComplement[idxRow][i-1])
                    lhs -= x[idxCol];
                else
                {
                    lhs += x[idxCol];
                    rhs += 1.0;
                }
            }
        }

        rhs -= 1.0;
        const double viol = lhs - rhs;

        /* cut is not a minimal cover or all items fit in the knapsack */
        if( (sumWeight - minWeight > rowRHS) || (sumWeight <= rowRHS) )
        {
            delete[] selected;
            continue;
        }

        /* invalid cut */
        if(viol < KSEP_DEF_MIN_VIOL)
        {
            delete[] selected;
            break;
        }

        /* extending cover cuts */
        if(ksep->rowHasCompVars[idxRow])
            simple_lifting(ksep, x, rc, idxRow, selected);
        else
        {
            clique_lifting(ksep, x, rc, idxRow, selected);
            simple_lifting(ksep, x, rc, idxRow, selected);
        }

        int cutSize = 0;
        for(int k = 0; k < nz; k++)
        {
            if(selected[k] != NOT_INSERTED)
            {
                const int idxCol = ksep->knpRowIdxs[idxRow][k];
                idx[cutSize] = idxCol;

                if(ksep->isComplement[idxRow][k])
                {
                    coef[cutSize] = -1.0;

                    if(selected[k] == INSERTED_AT_SIMPLE_EXTENSION)
                        rhs -= 1.0;

                    assert(selected[k] != INSERTED_AT_CLIQUE_EXTENSION);
                }
                else
                    coef[cutSize] = 1.0;

                cutSize++;
            }
        }

        Cut *cut = cut_create(idx, coef, cutSize, rhs, x);
        //#pragma omp critical(cut_pool_insertion)
        {
            if(!cut_pool_insert(cutPool, cut))
                cut_free(&cut);
        }

        delete[] selected;
    }

    delete[] idx;
    delete[] coef;

    if(bottomUp)
    {
    	for(int i = 0; i <= nz; i++)
        	delete[] table[i];
    	delete[] table;
    }
}

int knapsack_sep_separate(KnapsackSeparation *ksep, const double *x, const double *rc, CutPool *cutPool)
{
    //#pragma omp parallel for shared(ksep, x, cutPool)
    for(int i = 0; i < ksep->nRows; i++)
        separate_by_row(ksep, i, x, rc, cutPool);

    return 0;
}

void knapsack_sep_set_graph(KnapsackSeparation *ksep, const CGraph *cg)
{
    ksep->cg = cg;
}

void knapsack_sep_set_params_parse_cmd_line( KnapsackSeparation *ksep, const int argc, const char **argv )
{
    int i;

    char param[256] = "";
    char paramName[256] = "";
    char paramValue[256] = "";

    for ( i=1 ; ( i<argc ) ; ++i )
    {
        strncpy( param, argv[i], 256 );
        if ( strstr( param, "=" ) == NULL )
            continue;

        getParamName( paramName, param );
        getParamValue( paramValue, param );

        if( strcasecmp( KSEP_STR_MAX_RC, paramName ) == 0 )
        {
            ksep->maxRC = atof( paramValue );
            continue;
        }

        if( strcasecmp( KSEP_STR_MAX_IT_BK, paramName ) == 0 )
        {
            ksep->maxItBK = atof( paramValue );
            continue;
        }
    }
}

int knapsack_sep_get_num_cols(const KnapsackSeparation *ksep)
{
    return ksep->nCols;
}

int knapsack_sep_get_num_rows(const KnapsackSeparation *ksep)
{
    return ksep->nRows;
}
