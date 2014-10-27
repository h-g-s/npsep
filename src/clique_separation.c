#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "clique_separation.h"
#include "clique_extender.h"
#include "strUtils.h"
#include "memory.h"
#include "grasp.h"
#include "vectormgm.h"
#include "macros.h"
#include "clique_enum.h"
#include "bron_kerbosch.h"


/* default values */
#define CLQ_SEP_DEF_MAX_DEPTH      3
#define CLQ_SEP_DEF_MAX_PASSES     3
#define CLQ_SEP_DEF_MIN_VIOL    0.02
#define CLQ_SEP_DEF_MIN_FRAC   0.001
#define CLQ_SEP_DEF_EXTEND         2
#define CLQ_SEP_DEF_ENUM          10
#define CLQ_SEP_DEF_VERBOSE        1
#define CLQ_SEP_DEF_TLBK          20   /* time limit bron-kerbosch */

/* command line param names */
#define CLQ_SEP_STR_MAX_DEPTH   "maxDepth"
#define CLQ_SEP_STR_MAX_PASSES  "maxPasses"
#define CLQ_SEP_STR_MIN_VIOL    "minViol"
#define CLQ_SEP_STR_MIN_FRAC    "minFrac"
#define CLQ_SEP_STR_EXTEND      "extendC"
#define CLQ_SEP_STR_ENUM        "enumM"
#define CLQ_SEP_STR_VERBOSE     "verbose"


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

   CliqueEnumerator *clqEnum;
   IntSet clqEnumNeighs;  /* used in enumeration */

   /* clique set with original nodes, only translated cliques */
   CliqueSet *clqSetOrig;

   /* for each clique in clqSetOrig, stores if it was extended */
   char *extended;
   int extendedCap;

   /* final clique set, possibly with additional extended cliques */
   CliqueSet *clqSet;

   /* maxDepth: default 3 */
   int maxDepth;

   /* maxPasses: default 3 */
   int maxPasses;

   /* minViol: default 0.05 */
   double minViol;

   /* extendCliques: 0 - off  1 - extend randomly  2 - greedy extend (default) */
   int extendCliques;

   /* 0 - no enumeration    > 0 : enumerate all cliques containing node i such that d(i) < = enumUsage */
   int enumUsage;

   /* costs based in reduced cost for extending cliques */
   int *costs;
   char hasCosts;

   /* minimum fractional value that a variable
      must have to be considered for separation */
   double minFrac;

   void *bk;

   /* time limit execution for bron-kerbosch */
   double maxTimeBK;

   int verbose;
};

/* private function */
void clq_sep_check_node_cap( CliqueSeparation *clq_sep );

int clq_sep_enum_all_cliques_low_degree_nodes( CliqueSeparation *sep, int iv[], const CGraph *ppGraph );

CliqueSeparation *clq_sep_create( const CGraph *origGraph )
{
   CliqueSeparation *clqSep = xmalloc( sizeof(CliqueSeparation) );

   clqSep->maxDepth      = CLQ_SEP_DEF_MAX_DEPTH;
   clqSep->maxPasses     = CLQ_SEP_DEF_MAX_PASSES;
   clqSep->minViol       = CLQ_SEP_DEF_MIN_VIOL;
   clqSep->minFrac       = CLQ_SEP_DEF_MIN_FRAC;
   clqSep->extendCliques = CLQ_SEP_DEF_EXTEND;
   clqSep->enumUsage     = CLQ_SEP_DEF_ENUM;
   clqSep->verbose       = CLQ_SEP_DEF_VERBOSE;
   clqSep->maxTimeBK     = CLQ_SEP_DEF_TLBK;
   clqSep->hasCosts = 0;
   clqSep->nodeCap = cgraph_size(origGraph);
   clqSep->clqEnum = clq_enum_create( 0 );
   vint_set_init( &(clqSep->clqEnumNeighs) );

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
   printf("verbose mode now %d.\n", ((int)sep->verbose) );
   const CGraph *cgraph = sep->cgraph;

   clq_set_clear( sep->clqSet );
   clq_set_clear( sep->clqSetOrig );

   clq_sep_check_node_cap( sep );

   CliqueSet *clqSetOrig = sep->clqSetOrig; clq_set_clear( clqSetOrig ); /* before extension, orig indexes */
   CliqueSet *clqSet = sep->clqSet; clq_set_clear( clqSet ); /* final clique set */

   int *iv = sep->iv;
   const double minFrac = sep->minFrac;

   {
      int i, idx = 0;
      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
         iv[i] = -1;

      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
         if ((cgraph_degree(cgraph, i)<2)||(fracPart(x[i])<minFrac))
            continue;
         else
            iv[i] = idx++;
   }

   CGraph *ppcg = cgraph_create_induced_subgraph( cgraph, iv );
   cgraph_set_low_degree( ppcg, sep->enumUsage );
   clq_sep_update_ppgraph_weights( ppcg, cgraph_size(cgraph), x );

   /* if enumeration is used, iv will be update*/

   /*cgraph_save( ppcg, "ppgraph.clq" );*/

   if (sep->verbose)
      cgraph_print_summary( ppcg, "pre-processed graph - part 1" );

   /* separation works with integer weights */
   const int minW = (int)(1000.0 + (sep->minViol*1000.0));

   char enumerationComplete = 0;
   if ( sep->enumUsage > 0 )
   {
      clq_enum_set_min_weight( sep->clqEnum, minW );
      int enumNodes = clq_sep_enum_all_cliques_low_degree_nodes( sep, iv, ppcg );
      if (enumNodes < cgraph_size(ppcg))
      {
         /* more processing will be needed, creating a another preprocesses graph
           without already enumerated nodes */
         cgraph_free( &ppcg );
         ppcg = cgraph_create_induced_subgraph( cgraph, iv );
         cgraph_set_low_degree( ppcg, sep->enumUsage );
         clq_sep_update_ppgraph_weights( ppcg, cgraph_size(cgraph), x );
         if (sep->verbose)
            cgraph_print_summary( ppcg, "pre-processed graph - part 2" );
      }
      else
      {
         if (sep->verbose)
            printf("pre-processed graph - part 2 - no nodes left. all cliques have been enumerated.\n");
         enumerationComplete = 1;
      }
   }

   int firstGraspClique = clq_set_number_of_cliques( clqSetOrig );
   if ( (!enumerationComplete) && (cgraph_size(ppcg)>=2) )
   {
      sep->bk = bk_create( ppcg );
      clock_t startBK = clock();
      int stillWorkToDo = bk_run( sep->bk, minW, sep->maxTimeBK );
      clock_t endBK = clock();
      if (sep->verbose)
      {
         printf("bk took %.3g seconds\n", ((double)endBK-startBK)/((double)CLOCKS_PER_SEC) );
      }

      CliqueSet *bkClqSet = bk_get_clq_set(sep->bk);

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

      if (stillWorkToDo)
      {
         Grasp *grasp = grasp_create( ppcg, minW );
         if (sep->verbose)
         {
            printf("running grasp\n");
         }
         clock_t graspStart = clock();
         grasp_run( grasp );
         const CliqueSet *graspCliques = grasp_solution_set(grasp); /* before extension, pp indexes */
         clq_set_add_using_original_indexes( clqSetOrig, graspCliques, cgraph_get_original_node_indexes( ppcg ) );
         grasp_free( &grasp );
         clock_t graspEnd = clock();
         if (sep->verbose)
         {
            printf("grasp took %.3f seconds.\n", ((double)graspEnd-graspStart)/((double)CLOCKS_PER_SEC) );
         }
      }

      bk_free( (sep->bk) );
      sep->bk = NULL;
   }

   /* extending cliques */
   vmg_adjust_vector_capacity( (void**)&(sep->extended), &(sep->extendedCap), clq_set_number_of_cliques(clqSetOrig), sizeof(char) );
   char *extended = sep->extended;
   memset( extended, 0, sizeof(char)*clq_set_number_of_cliques( clqSetOrig ) );

   /* since grasp ran in a restricted subgraph, cliques found may be dominated by the ones
      found in the enumeration phase, checking ... */
   {
      int i;
      for ( i=firstGraspClique ; (i<clq_set_number_of_cliques(clqSetOrig)) ; ++i )
      {
         const IntSet *graspClique = clq_set_get_clique( clqSetOrig, i );
         int j;
         for ( j=0 ; (j<firstGraspClique) ; ++j )
            if ( clq_dominates( clq_set_get_clique( clqSetOrig, j ), graspClique ) )
               extended[i] = 1;
      }
   }

   if (sep->extendCliques)
   {
      clock_t startExtend = clock();

      CliqueExtender *clqe = sep->clqe;

      CliqueExtendingMethod clqem = CLQEM_RANDOM;
      if (sep->hasCosts)
      {
         clqe_set_costs( clqe, sep->costs, cgraph_size(cgraph) );
         clqem = CLQEM_PRIORITY_GREEDY;
      }


      int i;
      for ( i=0 ; (i<clq_set_number_of_cliques( clqSetOrig )) ; ++i )
         if (!extended[i])
            extended[i] = clqe_extend( clqe, cgraph, clq_set_get_clique(clqSetOrig,i), clq_set_weight(clqSetOrig,i), clqem );

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

   clq_enum_free( &((*clqSep)->clqEnum) );
   vint_set_clean( &((*clqSep)->clqEnumNeighs) );

   free(*clqSep);

   *clqSep = NULL;
}

void clq_sep_params_print( const CliqueSeparation *clqsp )
{
   printf( "Clique Separation Parameters:\n" );
   printf( "\tmax depth         : %d\n", clqsp->maxDepth );
   printf( "\tmax passes        : %d\n", clqsp->maxPasses );
   printf( "\tminimum violation : %.4f\n", clqsp->minViol );
   printf( "\tminimum frac. val : %.4f\n", clqsp->minFrac );
   printf( "\textend cliques    : %d\n", clqsp->extendCliques );
   printf( "\tenum max. degree  : %d\n", clqsp->enumUsage );
   printf( "\tverbose           : %d\n", clqsp->verbose );
   printf( "\tmaximum time bk   : %d\n", CLQ_SEP_DEF_TLBK );
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

      if ( strcasecmp( CLQ_SEP_STR_MAX_DEPTH, paramName ) == 0 )
      {
         clqsp->maxDepth = atoi( paramValue );
         continue;
      }

      if ( strcasecmp( CLQ_SEP_STR_MAX_PASSES, paramName ) == 0 )
      {
         clqsp->maxPasses = atoi( paramValue );
         continue;
      }

      if ( strcasecmp( CLQ_SEP_STR_MIN_VIOL, paramName ) == 0 )
      {
         clqsp->minViol = atof( paramValue );
         continue;
      }

      if ( strcasecmp( CLQ_SEP_STR_EXTEND, paramName ) == 0 )
      {
         clqsp->extendCliques = atoi( paramValue );
         continue;
      }

      if ( strcasecmp( CLQ_SEP_STR_ENUM, paramName ) == 0 )
      {
         clqsp->enumUsage = atoi( paramValue );
         continue;
      }

      if ( strcasecmp( CLQ_SEP_STR_VERBOSE, paramName ) == 0 )
      {
         clqsp->verbose = atoi( paramValue );
         continue;
      }
   }
#undef STR_SIZE
}

void clq_sep_params_help_cmd_line()
{
   printf("\t-%s=int     : sets maximum depth\n", CLQ_SEP_STR_MAX_DEPTH );
   printf("\t-%s=int     : sets maximum passes\n", CLQ_SEP_STR_MAX_PASSES );
   printf("\t-%s=float   : minimum violation for a cut to be selected\n", CLQ_SEP_STR_MIN_VIOL );
   printf("\t-%s=float   : minimum violation for a cut to be selected\n", CLQ_SEP_STR_MIN_FRAC );
   printf("\t-%s={0,1,2} : clique extension: off, random and greedy using reduced costs, respectively\n", CLQ_SEP_STR_EXTEND );
   printf("\t-%s=int     : clique enumeration: 0 off, >= 1 enumerate cliques containing nodes with degree at most \"int\"\n", CLQ_SEP_STR_ENUM );
   printf("\t-%s=int     : print detailed messages on-off\n", CLQ_SEP_STR_VERBOSE );
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

/* enumerated all cliques containing nodes with a maximum degree.
   after the enumeration, these nodes can be safely removed from the graph
   returns how many node were enumerated
   */
int clq_sep_enum_all_cliques_low_degree_nodes( CliqueSeparation *sep, int iv[], const CGraph *ppGraph )
{
   int i;
   const int nodes = cgraph_size( ppGraph );
   if (!nodes)
      return 0;
   CliqueEnumerator *clqEnum = sep->clqEnum;
   IntSet *neighs = &(sep->clqEnumNeighs);
   const int neighsCap = nodes*100;
   vint_set_check_capacity( neighs, neighsCap );
   int *vneighs = vint_set_force_elements_access( neighs );
   const int nCliquesBefore = clq_set_number_of_cliques( sep->clqSetOrig );
   int removed = 0;

   clock_t start = clock();

   for ( i=0 ; (i<nodes) ; i++ )
   {
      if ( cgraph_degree( ppGraph, i ) <= sep->enumUsage )
      {
         const int nodeNeighs = cgraph_get_all_conflicting( ppGraph, i, vneighs, neighsCap );
         const int nodeWeight = cgraph_get_node_weight( ppGraph, i );
         vint_set_force_size( neighs, nodeNeighs );
#ifdef DEBUG
         vint_set_force_check( neighs );
#endif
         if (neighs->size)
         {
            clq_enum_run( clqEnum, ppGraph, 1, &i, nodeWeight, neighs );
            clq_set_add_using_original_indexes( sep->clqSetOrig, clq_enum_get_cliques( clqEnum ), cgraph_get_original_node_indexes( ppGraph ) );
         }
         iv[ cgraph_get_original_node_index(ppGraph,i) ] = -1;
         ++removed;
      }
   }

   const int nCliquesAfter = clq_set_number_of_cliques( sep->clqSetOrig );
   clock_t end = clock();
   const double secs = (((double)end-start)/((double)CLOCKS_PER_SEC));

   if (sep->verbose)
   {
      printf("Enumeration finished in %.3f seconds. Found %d violated cliques. Removed %d nodes from graph.\n", secs, nCliquesAfter-nCliquesBefore, removed );
      fflush(stdout); fflush(stderr);
   }

   /* updating incidence vetor */
   int idx = 0;
   const int nodesOrig = cgraph_size( sep->cgraph );
   for ( i=0 ; (i<nodesOrig) ; ++i )
      if ( iv[i] != -1 )
         iv[i] = idx++;

#ifdef DEBUG
   const CliqueSet *clqSet = sep->clqSetOrig;
   const int nCliques = clq_set_number_of_cliques( clqSet );
   for ( i=0 ; (i<nCliques) ; ++i )
   {
      const IntSet *clique = clq_set_get_clique( clqSet, i );
      const int *el = vint_set_get_elements( clique );
      const int cSize = vint_set_size( clique );
      int n1,n2;
      if (!clq_validate( sep->cgraph, cSize, el, &n1, &n2 ))
      {
         fprintf( stderr, "Invalid clique found in enumeration. Nodes %d and %d are not neighbors.\n", n1, n2 );
         exit( EXIT_FAILURE );
      }
   }
#endif

   return removed;
}

int clq_sep_get_max_passes( const CliqueSeparation *clqSep )
{
   return clqSep->maxPasses;
}

double clq_sep_get_max_time_bk( const CliqueSeparation *clqSep )
{
   return clqSep->maxTimeBK;
}

void clq_sep_set_max_time_bk( CliqueSeparation *clqSep, const double maxTime )
{
   clqSep->maxTimeBK = maxTime;
}

void clq_sep_set_min_viol( CliqueSeparation *sep, const double viol )
{
   sep->minViol = viol;
}

void clq_sep_set_verbose( CliqueSeparation *sep, const char verbose )
{
   sep->verbose = verbose;
}

