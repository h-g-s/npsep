#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <float.h>
#include "clique_separation.h"
#include "memory.h"
#include "vectormgm.h"
#include "macros.h"
#include "bron_kerbosch.h"
#include "vint_queue.h"
#include "strutils.h"


/* default values */
#define CLQ_SEP_DEF_VERBOSE             0
#define CLQ_SEP_DEF_MIN_VIOL            0.02
#define CLQ_SEP_DEF_MIN_FRAC            0.001
#define CLQ_SEP_DEF_MAX_IT_BK           INT_MAX/10
#define CLQ_SEP_DEF_CLQE_EXTEND         CLQEM_PRIORITY_GREEDY
#define CLQ_SEP_DEF_CLQE_MAX_RC         DBL_MAX

/* command line param names */
#define CLQ_SEP_STR_VERBOSE              "verbose"
#define CLQ_SEP_STR_MIN_VIOL             "minViol"
#define CLQ_SEP_STR_MIN_FRAC             "minFrac"
#define CLQ_SEP_STR_MAX_IT_BK            "maxItBK"
#define CLQ_SEP_STR_CLQE_EXTEND          "extendC"
#define CLQ_SEP_STR_CLQE_MAX_CANDIDATES  "maxClqECand"
#define CLQ_SEP_STR_CLQE_MAX_RC          "maxRC"
#define CLQ_SEP_STR_CLQE_MAX_GEN         "maxClqEGen"
#define CLQ_SEP_STR_CLQE_RC_PERCENTAGE   "rcPercentage"

double fracPart( const double x );

struct _CliqueSeparation
{
    /* original conflict graph */
    const CGraph *cgraph;

    /* indicates if a variable will be considered in the pre-processed separation graph */
    int *iv;

    /* since cgraph may change of size, storing how many nodes we allocated space */
    int nodeCap;

    /* clique extender */
    CliqueExtender *clqe;

    /* clique set with original nodes, only translated cliques */
    CliqueSet *clqSetOrig;

    /* for each clique in clqSetOrig, stores if it was extended */
    char *extended;
    int extendedCap;

    /* final clique set, possibly with additional extended cliques */
    CliqueSet *clqSet;

    /* minViol: default 0.02 */
    double minViol;

    /* extendCliques: 0: off - 1: random (default) - 2: max degree selection - 3: greedy extension - 4: exact extension */
    int extendCliques;

    /* costs based in reduced cost for extending cliques */
    int *costs;
    char hasCosts;
    double maxRC;

    /* minimum fractional value that a variable must have to be considered for separation */
    double minFrac;

    void *bk;

    /* max iterations for bron-kerbosch */
    int maxItBK;

    int verbose;
};

/* private function */
void clq_sep_check_node_cap( CliqueSeparation *clq_sep );

CliqueSeparation *clq_sep_create( const CGraph *origGraph )
{
    CliqueSeparation *clqSep = xmalloc( sizeof(CliqueSeparation) );

    clqSep->minViol       = CLQ_SEP_DEF_MIN_VIOL;
    clqSep->minFrac       = CLQ_SEP_DEF_MIN_FRAC;
    clqSep->extendCliques = CLQ_SEP_DEF_CLQE_EXTEND;
    clqSep->verbose       = CLQ_SEP_DEF_VERBOSE;
    clqSep->maxItBK       = CLQ_SEP_DEF_MAX_IT_BK;
    clqSep->maxRC         = CLQ_SEP_DEF_CLQE_MAX_RC;
    clqSep->hasCosts = 0;
    clqSep->nodeCap = cgraph_size(origGraph);

    clqSep->iv = xmalloc( sizeof(int)*clqSep->nodeCap );
    clqSep->costs = xmalloc( sizeof(int)*clqSep->nodeCap );

    clqSep->cgraph = origGraph;

    clqSep->clqe = clqe_create();

    clqSep->clqSetOrig = clq_set_create();
    clqSep->clqSet = clq_set_create();

    clqSep->extendedCap = 1024;
    clqSep->extended = xmalloc( sizeof(char)*clqSep->extendedCap );

    return clqSep;
}

void clq_sep_set_rc( CliqueSeparation *sep, const double rc[] )
{
    clq_sep_check_node_cap( sep );
    int i;
    for ( i=0 ; (i<cgraph_size(sep->cgraph)) ; ++i )
        sep->costs[i] = (int)((rc[i]*1000.0)+0.5);

    sep->hasCosts = 1;
}

void clq_sep_update_ppgraph_weights( CGraph *ppcg, const int cols, const double x[] )
{
    /* weights for fractional variables */
    int i;
    for ( i=0 ; (i<cgraph_size(ppcg)) ; ++i )
    {
        const int origIdx = cgraph_get_original_node_index(ppcg,i);
        const double f = x[origIdx];
        cgraph_set_node_weight( ppcg, i, cgraph_weight( f ) );
    }
}

void clq_sep_separate( CliqueSeparation *sep, const double x[] )
{
    const CGraph *cgraph = sep->cgraph;

    clq_set_clear( sep->clqSet );
    clq_set_clear( sep->clqSetOrig );

    clq_sep_check_node_cap( sep );

    CliqueSet *clqSetOrig = sep->clqSetOrig;
    clq_set_clear( clqSetOrig ); /* before extension, orig indexes */
    CliqueSet *clqSet = sep->clqSet;
    clq_set_clear( clqSet ); /* final clique set */

    int *iv = sep->iv;
    const double minFrac = sep->minFrac;

    {
        IntQueue queue;
        int i, *neighs, csize = cgraph_size(cgraph);
        char inQueue[csize]; //'1' if a vertices is already in queue - '0' otherwise
        neighs = xmalloc(sizeof(int) * csize * 10);
        vint_queue_init(&queue, csize);

        for(i = 0; i < csize; i++)
        {
            if(cgraph_degree(cgraph, i) == 0) //deleting variables that are not binary
            {
                iv[i] = -1;
                inQueue[i] = '0'; //just to avoid memory problems
            }

            else if(cgraph_degree(cgraph, i) == 1 || (fracPart(x[i]) < minFrac && x[i] < 0.99))
            { //integer variables and variables that have conflict only with their complements are deleted
                iv[i] = -1;
                vint_queue_push(&queue, i);
                inQueue[i] = '1';
            }

            else
            {
                iv[i] = cgraph_degree(cgraph, i);
                inQueue[i] = '0';
            }
        }

        while(!vint_queue_is_empty(&queue))
        {
            int v;
            vint_queue_pop(&queue, &v);
            int nsize = cgraph_get_all_conflicting(cgraph, v, neighs, csize * 10);

            for(i = 0; i < nsize; i++)
            {
                int u = neighs[i];

                if(iv[u] == -1) continue;

                assert(iv[u] > 0);
                iv[u]--; //considering v was deleted

                if(iv[u] == 0)
                {
                    iv[u] = -1;

                    if(inQueue[u] == '0')
                    {
                        vint_queue_push(&queue, u);
                        inQueue[u] = '1';
                    }
                }
            }
        }

        int idx = 0;
        for(i = 0; i < csize; i++)
            if(iv[i] > 0)
                iv[i] = idx++;

        free(neighs);
        vint_queue_clean(&queue);
    }


    CGraph *ppcg = cgraph_create_induced_subgraph( cgraph, iv );
    clq_sep_update_ppgraph_weights( ppcg, cgraph_size(cgraph), x );

    if (sep->verbose)
        cgraph_print_summary( ppcg, "pre-processed graph - part 1" );

    /* separation works with integer weights */
    const int minW = (int)(1000.0 + (sep->minViol*1000.0));

    if(cgraph_size(ppcg)>=2)
    {
        sep->bk = bk_create( ppcg );
        clock_t startBK = clock();
        bk_set_max_it(sep->bk, sep->maxItBK);
        bk_set_min_weight(sep->bk, minW);
        bk_run( sep->bk );
        clock_t endBK = clock();
        if (sep->verbose)
        {
            printf("bk took %.3g seconds\n", ((double)endBK-startBK)/((double)CLOCKS_PER_SEC) );
        }

        const CliqueSet *bkClqSet = bk_get_clq_set(sep->bk);

        if (bkClqSet)
        {
            if (clq_set_number_of_cliques( bkClqSet ))
            {
#ifdef DEBUG
                int nc = clq_set_number_of_cliques( bkClqSet );
                int ic;
                for ( ic = 0 ; (ic<nc) ; ++ic )
                {
                    const IntSet *is = clq_set_get_clique( bkClqSet, ic );
                    int n1, n2;
                    if (!clq_validate( ppcg, vint_set_size(is), vint_set_get_elements(is), &n1, &n2 ))
                    {
                        fprintf( stderr, "Nodes %d and %d are not in conflict in ppcg.\n", n1, n2 );
                        exit( EXIT_FAILURE );
                    }
                    int j;
                    for ( j=0 ; (j<vint_set_size(is)) ; ++j )
                    {
                        const int vidx = vint_set_get_elements(is)[j];
                        assert( vidx >=0 );
                        assert( vidx < cgraph_size(ppcg) );
                    }
                }
#endif
                clq_set_add_using_original_indexes( clqSetOrig, bkClqSet , cgraph_get_original_node_indexes( ppcg ) );
            }
        }

        bk_free( (sep->bk) );
        sep->bk = NULL;
    }

    /* extending cliques */
    vmg_adjust_vector_capacity( (void**)&(sep->extended), &(sep->extendedCap), clq_set_number_of_cliques(clqSetOrig), sizeof(char) );
    char *extended = sep->extended;
    memset( extended, 0, sizeof(char)*clq_set_number_of_cliques( clqSetOrig ) );

    if (sep->extendCliques)
    {
        clock_t startExtend = clock();

        CliqueExtender *clqe = sep->clqe;

        if(sep->hasCosts)
            clqe_set_costs( clqe, sep->costs, cgraph_size(cgraph) );

        int i;
        for ( i=0 ; (i<clq_set_number_of_cliques( clqSetOrig )) ; ++i )
            extended[i] = clqe_extend( clqe, cgraph, clq_set_get_clique(clqSetOrig,i), clq_set_weight(clqSetOrig,i), sep->extendCliques );

        /* adding all extended cliques */
        clq_set_add_cliques( clqSet, clqe_get_cliques( clqe ) );
        clock_t endExtend = clock();
        const double timeExtend = ((double)endExtend-startExtend) /  ((double)CLOCKS_PER_SEC);
        if (sep->verbose)
        {
            printf("clique extension took %.3f seconds.\n", timeExtend);
        }
    }

    /* adding cliques which were not extended */
    {
        int i;
        for ( i=0 ; (i<clq_set_number_of_cliques(clqSetOrig)) ; ++i )
            if ( !extended[i] )
                clq_set_add( clqSet, clq_set_clique_size(clqSetOrig,i), clq_set_clique_elements(clqSetOrig,i), clq_set_weight(clqSetOrig,i) );
    }

    /* need to be informed again next call */
    sep->hasCosts = 0;

    cgraph_free( &ppcg );
}

const CliqueSet *clq_sep_get_cliques( CliqueSeparation *sep )
{
    return sep->clqSet;
}

void clq_sep_free( CliqueSeparation **clqSep )
{
    clqe_free( &((*clqSep)->clqe) );

    clq_set_free( &((*clqSep)->clqSetOrig) );
    clq_set_free( &((*clqSep)->clqSet) );

    free( (*clqSep)->iv );
    free( (*clqSep)->costs );
    free( (*clqSep)->extended );

    free(*clqSep);

    *clqSep = NULL;
}

int clq_sep_get_verbose( CliqueSeparation *sep )
{
    return sep->verbose;
}

double clq_sep_get_min_viol( CliqueSeparation *sep )
{
    return sep->minViol;
}

void clq_sep_set_params_parse_cmd_line( CliqueSeparation *clqsp, const int argc, const char **argv )
{
    int i;

#define STR_SIZE 256
    char param[STR_SIZE] = "";
    char paramName[STR_SIZE] = "";
    char paramValue[STR_SIZE] = "";

    for ( i=1 ; ( i<argc ) ; ++i )
    {
        strncpy( param, argv[i], STR_SIZE );
        if ( strstr( param, "=" ) == NULL )
            continue;

        getParamName( paramName, param );
        getParamValue( paramValue, param );

        if ( strcasecmp( CLQ_SEP_STR_VERBOSE, paramName ) == 0 )
        {
            clqsp->verbose = atoi( paramValue );
            continue;
        }

        if ( strcasecmp( CLQ_SEP_STR_MIN_VIOL, paramName ) == 0 )
        {
            clqsp->minViol = atof( paramValue );
            continue;
        }

        if( strcasecmp( CLQ_SEP_STR_MIN_FRAC, paramName ) == 0 )
        {
            clqsp->minFrac = atof( paramValue );
        }

        if( strcasecmp( CLQ_SEP_STR_MAX_IT_BK, paramName ) == 0 )
        {
            clqsp->maxItBK = atoi( paramValue );
        }

        if ( strcasecmp( CLQ_SEP_STR_CLQE_EXTEND, paramName ) == 0 )
        {
            clqsp->extendCliques = atoi( paramValue );
            continue;
        }

        if ( strcasecmp( CLQ_SEP_STR_CLQE_MAX_CANDIDATES, paramName ) == 0 )
        {
            clqe_set_max_candidates(clqsp->clqe, atoi( paramValue ));
            continue;
        }

        if( strcasecmp( CLQ_SEP_STR_CLQE_MAX_RC, paramName ) == 0 )
        {
            clqsp->maxRC = atof( paramValue );
            clqe_set_max_cost(clqsp->clqe, (int)((clqsp->maxRC*1000.0)+0.5));
        }

        if( strcasecmp( CLQ_SEP_STR_CLQE_MAX_GEN, paramName ) == 0 )
        {
            clqe_set_max_clq_gen(clqsp->clqe, atoi( paramValue ));
        }

        if( strcasecmp( CLQ_SEP_STR_CLQE_RC_PERCENTAGE, paramName ) == 0 )
        {
            clqe_set_rc_percentage(clqsp->clqe, atof( paramValue ));
        }
    }
#undef STR_SIZE
}

void clq_sep_params_print( const CliqueSeparation *clqsp )
{
    printf( "Clique Separation Parameters:\n" );
    printf( "\tverbose                   : %d\n", clqsp->verbose );
    printf( "\tminimum violation         : %.4lf\n", clqsp->minViol );
    printf( "\tminimum frac. val         : %.4lf\n", clqsp->minFrac );
    printf( "\tmaximum iterations bk     : %d\n", clqsp->maxItBK );
    printf( "\textend cliques            : %d\n", clqsp->extendCliques );
    printf( "\tmaximum rc                : %.4lf\n", clqsp->maxRC);
    printf( "\tmaximum clqs in exact sep : %d\n", clqe_get_max_clq_gen(clqsp->clqe) );
    printf( "\tmaximum range of rc       : %.4lf\n", clqe_get_rc_percentage(clqsp->clqe));
    fflush(stdout);
}

void clq_sep_params_help_cmd_line()
{
    printf("  -%s=int           : print detailed messages on-off\n", CLQ_SEP_STR_VERBOSE );
    printf("  -%s=float         : minimum violation for a cut to be selected\n", CLQ_SEP_STR_MIN_VIOL );
    printf("  -%s=float         : minimum violation for a cut to be selected\n", CLQ_SEP_STR_MIN_FRAC );
    printf("  -%s=int           : max iterations for BK algorithm\n", CLQ_SEP_STR_MAX_IT_BK );
    printf("  -%s={0,1,2,3,4}   : clique extension: off, random, max degree, greedy and exact, respectively.\n", CLQ_SEP_STR_CLQE_EXTEND );
    printf("  -%s=float           : maximum reduced cost value that a variable must have to be considered for extension.\n", CLQ_SEP_STR_CLQE_MAX_RC );
    printf("  -%s=int        : maximum number of cliques generated by the exact extension of an original one.\n", CLQ_SEP_STR_CLQE_MAX_GEN );
    printf("  -%s=float    : maximum range of reduced cost value considered in the exact extension.\n", CLQ_SEP_STR_CLQE_RC_PERCENTAGE );
    fflush(stdout);
}

void clq_sep_check_node_cap( CliqueSeparation *clq_sep )
{
    const int cGraphSize = cgraph_size( clq_sep->cgraph );

    if ( cGraphSize > clq_sep->nodeCap )
    {
        clq_sep->nodeCap = MAX( cGraphSize, clq_sep->nodeCap*2 );
        clq_sep->iv = xrealloc( clq_sep->iv, sizeof(int)*clq_sep->nodeCap );
        clq_sep->costs = xrealloc( clq_sep->costs, sizeof(int)*clq_sep->nodeCap );
    }
}

double fracPart( const double x )
{
    double nextInteger = ceil( x );
    double previousInteger = floor( x );

    return MIN( nextInteger-x, x-previousInteger );
}

double clq_sep_get_max_it_bk( const CliqueSeparation *clqSep )
{
    return clqSep->maxItBK;
}

void clq_sep_set_max_it_bk( CliqueSeparation *clqSep, size_t maxItBk )
{
    clqSep->maxItBK = maxItBk;
}

void clq_sep_set_min_viol( CliqueSeparation *sep, const double viol )
{
    sep->minViol = viol;
}

void clq_sep_set_verbose( CliqueSeparation *sep, const char verbose )
{
    sep->verbose = verbose;
}

void clq_sep_set_extend_method( CliqueSeparation *sep, const int extendC )
{
    assert(extendC >= 0 && extendC <= 2);
    sep->extendCliques = extendC;
}

int clq_sep_get_extend_method( CliqueSeparation *sep )
{
    return sep->extendCliques;
}
