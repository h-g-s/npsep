#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "cgraph.h"
#include "memory.h"
#include "clique.h"
#include "vint_set.h"
#include "clique_extender.h"
#include "vectormgm.h"
#include "macros.h"
#include "bron_kerbosch.h"

#define CLQE_DEF_MAX_CANDIDATES   256
#define CLQE_DEF_MAX_COST         500.0
#define CLQE_DEF_MAX_GEN          5
#define CLQE_DEF_RC_PERCENTAGE    0.6

#define CLQE_DEF_BK_MAX_IT        5000

/* additional space for candidates */
#define CANDIDATES_SLACK 100

struct _CliqueExtender
{
    const CGraph *cgraph;
    CliqueSet *clqSet;
    int *candidates;
    int *newClique;
    int newCliqueSize;
    int candidatesCap;

    int maxCandidates;
    int maxCost;
    int maxClqGen;
    double rcPercentage;

    NeighIterator *nit;

    int costsCap;
    int *costs;
};

void clqe_check_nodes_cap( CliqueExtender *clqe )
{
    const CGraph *cgraph = clqe->cgraph;

    const int nodes = cgraph_size( cgraph );

    if ( clqe->candidatesCap < nodes )
    {
        clqe->candidatesCap = nodes;
        if ( !clqe->candidates )
        {
            clqe->candidates = xmalloc( sizeof(int)*nodes*CANDIDATES_SLACK );
            clqe->newClique = xmalloc( sizeof(int)*nodes );
        }
        else
        {
            clqe->candidates = xrealloc( clqe->candidates, sizeof(int)*clqe->candidatesCap*CANDIDATES_SLACK );
            clqe->newClique = xrealloc( clqe->newClique, sizeof(int)*clqe->candidatesCap );
        }
    }
}

CliqueExtender *clqe_create()
{
    CliqueExtender *clqe = xmalloc( sizeof(CliqueExtender) );

    clqe->cgraph  = NULL;
    clqe->clqSet  = clq_set_create();
    clqe->candidatesCap = 0;
    clqe->candidates = NULL;
    clqe->newClique = NULL;

    clqe->costs = NULL;
    clqe->costsCap = 0;

    clqe->nit = nit_create();

    clqe->maxCandidates = CLQE_DEF_MAX_CANDIDATES;
    clqe->maxCost = CLQE_DEF_MAX_COST;
    clqe->maxClqGen = CLQE_DEF_MAX_GEN;
    clqe->rcPercentage = CLQE_DEF_RC_PERCENTAGE;

    return clqe;
}

int clqe_get_best_candidates_clique_insertion( CliqueExtender *clqe, const IntSet *clique,
        const CliqueExtendingMethod clqem )
{
    /* node with the smallest degree */
    int nodeSD = -1, degree = INT_MAX, i;

    const int cliqueSize = vint_set_size( clique ), *cliqueEl = vint_set_get_elements( clique );

    const CGraph *cgraph = clqe->cgraph;
    const int cgSize = cgraph_size(cgraph);
    const int nCols = cgSize/2;
    assert( (cgSize % 2) == 0 );
    char *iv = xcalloc(cgSize, sizeof(char));

    /* picking node with the smallest degree */
    for ( i=0 ; (i<cliqueSize) ; ++i )
    {
        if ( cgraph_degree( cgraph, cliqueEl[i] ) < degree )
        {
            degree = cgraph_degree( cgraph, cliqueEl[i] );
            nodeSD = cliqueEl[i];
        }
    }

    int nCandidates = 0;
    int *candidates = clqe->candidates;
    const int *costs = clqe->costs;

    assert(clqem == CLQEM_PRIORITY_GREEDY || clqem == CLQEM_RANDOM || clqem == CLQEM_MAX_DEGREE);

    if(clqem == CLQEM_PRIORITY_GREEDY) //clique extender method uses greedy selection (reduced cost)
    {
        NeighIterator *nit = clqe->nit;
        int selected = -1;

        nit_start( nit, cgraph, nodeSD, costs  );

        while ( ( (selected=nit_next(nit))!=INT_MAX ) && (nCandidates<clqe->maxCandidates) )
        {
            if( (costs != NULL && costs[selected] > clqe->maxCost) || iv[selected] == 1)
                continue;

            /* need to have conflict with all nodes in clique and all others inserted */
            for(i = 0; i < cliqueSize; i++)
            {
                int complement = (cliqueEl[i] < nCols) ? (cliqueEl[i] + nCols) : (cliqueEl[i] - nCols);
                if ( (!cgraph_conflicting_nodes( cgraph, cliqueEl[i], selected )) || (selected==cliqueEl[i]) || (selected==complement) )
                    break;
            }
            if (i<cliqueSize)
                continue;
            for ( i=0 ; (i<nCandidates) ; ++i )
                if (!cgraph_conflicting_nodes( cgraph, candidates[i], selected ))
                    break;
            if (i<nCandidates)
                continue;
            candidates[nCandidates++] = selected;

            int complement = (selected < nCols) ? (selected + nCols) : (selected - nCols);
            iv[selected] = 1;
            iv[complement] = 1;
        }
    }
    else if(clqem == CLQEM_MAX_DEGREE)
    {
        int *degree = xmalloc(sizeof(int) * cgraph_size(cgraph));
        NeighIterator *nit = clqe->nit;
        int selected = -1;

        for(i = 0; i < cgraph_size(cgraph); i++)
            degree[i] = cgraph_max_degree(cgraph) - cgraph_degree(cgraph, i);

        nit_start( nit, cgraph, nodeSD, degree  );

        while ( ( (selected=nit_next(nit))!=INT_MAX ) && (nCandidates<clqe->maxCandidates) )
        {
            if( (costs != NULL && costs[selected] > clqe->maxCost) || iv[selected] == 1)
                continue;

            for(i = 0; i < cliqueSize; i++)
            {
                int complement = (cliqueEl[i] < nCols) ? (cliqueEl[i] + nCols) : (cliqueEl[i] - nCols);
                if ( (!cgraph_conflicting_nodes( cgraph, cliqueEl[i], selected )) || (selected==cliqueEl[i]) || (selected==complement) )
                    break;
            }
            if (i < cliqueSize)
                continue;
            for(i = 0; i < nCandidates; i++)
                if (!cgraph_conflicting_nodes( cgraph, candidates[i], selected ))
                    break;
            if (i < nCandidates)
                continue;
            candidates[nCandidates++] = selected;

            int complement = (selected < nCols) ? (selected + nCols) : (selected - nCols);
            iv[selected] = 1;
            iv[complement] = 1;
        }
        free(degree);
    }
    else //clique extender method uses random selection
    {
        int *neighs = xmalloc(sizeof(int) * cgraph_size(cgraph) * 2);
        int nConflicts = cgraph_get_all_conflicting( cgraph, nodeSD, neighs, cgraph_size(cgraph) * 2);
        int j, selected;

        if(nConflicts < clqe->maxCandidates)
        {
            for(i = 0; i < nConflicts; i++)
            {
                selected = neighs[i];

                if( (costs != NULL && costs[selected] > clqe->maxCost) || iv[selected] == 1)
                    continue;

                for(j = 0; j < cliqueSize; j++)
                {
                    int complement = (cliqueEl[j] < nCols) ? (cliqueEl[j] + nCols) : (cliqueEl[j] - nCols);
                    if ( (!cgraph_conflicting_nodes( cgraph, cliqueEl[j], selected )) || (selected==cliqueEl[j]) || (selected==complement) )
                        break;
                }
                if(j < cliqueSize)
                    continue;
                for(j = 0; j < nCandidates; j++)
                    if(!cgraph_conflicting_nodes(cgraph, candidates[j], selected))
                        break;
                if(j < nCandidates)
                    continue;
                candidates[nCandidates++] = selected;

                int complement = (selected < nCols) ? (selected + nCols) : (selected - nCols);
                iv[selected] = 1;
                iv[complement] = 1;
            }
        }
        else
        {
            int r, remaining = nConflicts;
            char *isSelected = xmalloc(sizeof(char) * nConflicts);
            memset(isSelected, 0, sizeof(char) * nConflicts);

            while(nCandidates < clqe->maxCandidates && remaining > 0)
            {
                do {
                    r = rand() % nConflicts;
                    selected = neighs[r];
                } while(isSelected[r]);

                isSelected[r] = 1;
                remaining--;

                if( (costs != NULL && costs[selected] > clqe->maxCost) || iv[selected] == 1)
                    continue;

                for(j = 0; j < cliqueSize; j++)
                {
                    int complement = (cliqueEl[j] < nCols) ? (cliqueEl[j] + nCols) : (cliqueEl[j] - nCols);
                    if ( (!cgraph_conflicting_nodes( cgraph, cliqueEl[j], selected )) || (selected==cliqueEl[j]) || (selected==complement) )
                        break;
                }
                if(j < cliqueSize)
                    continue;
                for(j = 0; j < nCandidates; j++)
                    if(!cgraph_conflicting_nodes(cgraph, candidates[j], selected))
                        break;
                if(j < nCandidates)
                    continue;
                candidates[nCandidates++] = selected;

                int complement = (selected < nCols) ? (selected + nCols) : (selected - nCols);
                iv[selected] = 1;
                iv[complement] = 1;
            }
            free(isSelected);
        }

        free(neighs);
    }

    free(iv);

    return nCandidates;
}

typedef struct {
    int idx;
    int weight;
} CliqueWeight;

int cmp_clq_weight( const void *e1, const void *e2 )
{
    CliqueWeight i1 = (*((const CliqueWeight*) e1));
    CliqueWeight i2 = (*((const CliqueWeight*) e2));
    if(i1.weight != i2.weight)
        return i2.weight - i1.weight;
    return i1.idx - i2.idx;
}

int exact_clique_extension( CliqueExtender *clqe, const IntSet *clique, const int weight )
{
    int i, j, nCandidates = 0, newCliques = 0, cgSize = cgraph_size(clqe->cgraph);
    int *nindexes = xmalloc(sizeof(int) * cgSize);
    int *candidates = xmalloc(sizeof(int) * cgSize);
    const int clqSize = vint_set_size(clique);
    const int *clqEl = vint_set_get_elements(clique);
    int nodeSD = -1, degree = INT_MAX;
    int minRC =  INT_MAX, maxRC = -INT_MAX;
    int *newClique = xmalloc(sizeof(int) * cgSize);
    const int nCols = cgSize/2;
    assert( (cgSize % 2) == 0 );
    char *iv = xcalloc(cgSize, sizeof(char));

    for(i = 0; i < cgSize; i++)
        nindexes[i] = -1;

    for(i = 0; i < clqSize; i++)
    {
        newClique[i] = clqEl[i];
        if(cgraph_degree(clqe->cgraph, clqEl[i]) < degree)
        {
            nodeSD = clqEl[i];
            degree = cgraph_degree(clqe->cgraph, clqEl[i]);
        }
    }

    int *neighs = xmalloc(sizeof(int) * cgSize * 2);
    int n = cgraph_get_all_conflicting(clqe->cgraph, nodeSD, neighs, cgSize*2);

    for(i = 0; i < n; i++)
    {
        int neigh = neighs[i];

        if(clqe->costs[neigh] > clqe->maxCost || iv[neigh] == 1)
            continue;

        for(j = 0; j < clqSize; j++)
        {
            int complement = (clqEl[j] < nCols) ? (clqEl[j] + nCols) : (clqEl[j] - nCols);
            if( (neigh == clqEl[j]) || (neigh == complement) || (!cgraph_conflicting_nodes(clqe->cgraph, clqEl[j], neigh)) )
                break;
        }
        if(j >= clqSize)
        {
            candidates[nCandidates++] = neigh;
            minRC = MIN(minRC, clqe->costs[neigh]);
            maxRC = MAX(maxRC, clqe->costs[neigh]);

            int complement = (neigh < nCols) ? (neigh + nCols) : (neigh - nCols);
            iv[neigh] = 1;
            iv[complement] = 1;
        }
    }

    int maxValue = minRC + floor(clqe->rcPercentage * (maxRC-minRC));
    j = nCandidates;
    nCandidates = 0;
    for(i = 0; i < j; i++)
        if(clqe->costs[candidates[i]] > maxValue)
            nindexes[candidates[i]] = -1;
        else nindexes[candidates[i]] = nCandidates++;

    free(neighs);
    free(candidates);

    if(nCandidates == 0)
    {
        free(newClique);
        free(nindexes);
        free(iv);
        return 0;
    }

    CGraph *cg = cgraph_create_induced_subgraph(clqe->cgraph, nindexes);
    free(nindexes);

    if(minRC == maxRC)
        for(i = 0; i < nCandidates; i++)
            cgraph_set_node_weight(cg, i, 1);
    else
        for(i = 0; i < nCandidates; i++)
        {
            int origIdx = cgraph_get_original_node_index(cg, i);
            double normRC = 1.0 - ( ((double)(clqe->costs[origIdx] - minRC)) / ((double)(maxRC - minRC)) );
            cgraph_set_node_weight(cg, i, cgraph_weight(normRC));
        }

    BronKerbosch *bk = bk_create(cg);
    bk_set_max_it(bk, CLQE_DEF_BK_MAX_IT);
    bk_set_min_weight(bk, 0);
    bk_run(bk);
    const CliqueSet *clqSet = bk_get_clq_set(bk);
    const int numClqs = clq_set_number_of_cliques(clqSet);
    if(clqSet && numClqs)
    {
        if(numClqs <= clqe->maxClqGen)
        {
            for(i = 0; i < numClqs; i++)
            {
                const int *extClqEl = clq_set_clique_elements(clqSet, i);
                const int extClqSize = clq_set_clique_size(clqSet, i);
                int newCliqueSize = clqSize;

                for(j = 0; j < extClqSize; j++)
                {
                    int origIdx = cgraph_get_original_node_index(cg, extClqEl[j]);
                    newClique[clqSize + j] = origIdx;
                    newCliqueSize++;
                }
                newCliques += clq_set_add(clqe->clqSet, newCliqueSize, newClique, weight);
            }
        }
        else
        {
            CliqueWeight *clqw = xmalloc(sizeof(CliqueWeight) * numClqs);
            for(i = 0; i < numClqs; i++)
            {
                clqw[i].idx = i;
                clqw[i].weight = clq_set_clique_size(clqSet, i) * clq_set_weight(clqSet, i);
            }
            qsort(clqw, numClqs, sizeof(CliqueWeight), cmp_clq_weight);
            for(i = 0; i < clqe->maxClqGen; i++)
            {
                const int *extClqEl = clq_set_clique_elements(clqSet, clqw[i].idx);
                const int extClqSize = clq_set_clique_size(clqSet, clqw[i].idx);
                int newCliqueSize = clqSize;

                for(j = 0; j < extClqSize; j++)
                {
                    int origIdx = cgraph_get_original_node_index(cg, extClqEl[j]);
                    newClique[clqSize + j] = origIdx;
                    newCliqueSize++;
                }
                newCliques += clq_set_add(clqe->clqSet, newCliqueSize, newClique, weight);
            }
            free(clqw);
        }
    }

    bk_free(bk);
    cgraph_free(&cg);
    free(newClique);
    free(iv);

    return newCliques;
}

int clqe_extend( CliqueExtender *clqe, const CGraph *cgraph, const IntSet *clique,
                 const int weight, const CliqueExtendingMethod clqem )
{
    int result = 0;

#ifdef DEBUG
    assert( (clqe) && (cgraph) && (clique) && ((vint_set_size(clique))) );
    int en1, en2;
    if (!clq_validate( cgraph, vint_set_size(clique), vint_set_get_elements(clique), &en1, &en2 ))
    {
        fprintf( stderr, "ERROR clqe_extend : Nodes %d and %d are not in conflict.\n", en1, en2 );
        exit( EXIT_FAILURE );
    }
#endif

    clqe->cgraph = cgraph;

    clqe_check_nodes_cap( clqe );

    if(clqem == CLQEM_EXACT)
    {
        if( (clqe->costs == NULL) || (clqe->costsCap<cgraph_size(cgraph)) )
        {
            fprintf(stderr, "Cannot start exact clique extension. No costs were informed.\n");
            exit(EXIT_FAILURE);
        }
        result += exact_clique_extension(clqe, clique, weight);
        return result;
    }

    if( (clqem == CLQEM_PRIORITY_GREEDY) && ((clqe->costs == NULL) || (clqe->costsCap < cgraph_size(cgraph))) )
        fprintf( stderr, "Warning: using random selection for extension since no costs were informed.\n");

    int nCandidates = clqe_get_best_candidates_clique_insertion( clqe, clique, clqem );
    assert(nCandidates >= 0 && nCandidates <= clqe->maxCandidates);

    if(!nCandidates)
        return 0;

    /* clique can be extended, starting to fill new clique */
    memcpy( clqe->newClique, vint_set_get_elements( clique ), sizeof(int)*vint_set_size(clique) );
    clqe->newCliqueSize = vint_set_size( clique );

    int i;
    for(i = 0; i < nCandidates; i++)
    {
        const int selectedNode = clqe->candidates[i];
        assert( selectedNode >= 0 && selectedNode < cgraph_size(cgraph) );
        clqe->newClique[clqe->newCliqueSize++] = selectedNode;
    }

    result += clq_set_add( clqe->clqSet, clqe->newCliqueSize, clqe->newClique, weight );

    return result;
}

void clqe_set_clear( CliqueExtender *clqe )
{
    clq_set_clear( clqe->clqSet );
}

const CliqueSet *clqe_get_cliques( CliqueExtender *clqe )
{
    return clqe->clqSet;
}

void clqe_set_costs( CliqueExtender *clqe, const int costs[], const int n )
{
    vmg_adjust_vector_capacity( (void**)&(clqe->costs), &(clqe->costsCap), n, sizeof(int) );

    memcpy( clqe->costs, costs, sizeof(int)*n );

    if ( clqe->costsCap > n )
        memset( clqe->costs+n, 0, (clqe->costsCap-n)*sizeof(int) );
}

const int *clqe_get_costs( CliqueExtender *clqe )
{
    return clqe->costs;
}

void clqe_free( CliqueExtender  **clqe )
{
    if ( (*clqe)->newClique )
        free( (*clqe)->newClique );
    if ( (*clqe)->candidates )
        free( (*clqe)->candidates );

    clq_set_free( &(((*clqe)->clqSet) ) );

    if ( (*clqe)->costs )
        free( (*clqe)->costs );

    nit_free( &((*clqe)->nit) );

    free( *clqe );
    (*clqe) = NULL;
}

void clqe_set_max_candidates( CliqueExtender *clqe, const int max_size )
{
    assert(max_size > 0);
    clqe->maxCandidates = max_size;
}

int clqe_get_max_candidates( CliqueExtender *clqe )
{
    return clqe->maxCandidates;
}

void clqe_set_max_cost( CliqueExtender *clqe, const int maxCost )
{
    clqe->maxCost = maxCost;
}

int clqe_get_max_cost( CliqueExtender *clqe )
{
    return clqe->maxCost;
}

void clqe_set_max_clq_gen( CliqueExtender *clqe, const int maxClqGen )
{
    assert(maxClqGen > 0);
    clqe->maxClqGen = maxClqGen;
}

int clqe_get_max_clq_gen( CliqueExtender *clqe )
{
    return clqe->maxClqGen;
}

void clqe_set_rc_percentage( CliqueExtender *clqe, const double rcPercentage )
{
    clqe->rcPercentage = rcPercentage;
}

double clqe_get_rc_percentage( CliqueExtender *clqe )
{
    return clqe->rcPercentage;
}
