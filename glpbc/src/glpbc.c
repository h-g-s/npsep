#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
#include <glpk.h>
#include <time.h>
#include "cgraph.h"
#include "grasp.h"
#include "conflict_discover.h"
#include "memory.h"
#include "glp_hooks.h"
#include "strUtils.h"
#include "clique_extender.h"
#include "vectormgm.h"
#include "lp_callbacks.h"

#define MAX( a, b ) (( (a)>(b) ) ? (a) : (b))
#define MIN( a, b ) (( (a)<(b) ) ? (a) : (b))
#define MIN_FRAC 1e-4
#define MAX_EXT 128

void cut_callback( glp_tree *tree, void *info );

/*static int callsCGRAPH = 0;*/
static double timeCGRAPH = 0.0;

typedef struct
{
   CGraph *cgraph;
   ConflictDiscover *cd;
   int *idxExtension;
   int *newNodes;
   int newNodesCap;

   int *idx;
   int idxCap;

   int *callsAtDepth;
   int lastDepth;
   int lastNode;
   int pass;
   int *ivExtended;
   int ivExtendedCap;
   CliqueExtender *clqe;
   CliqueSet *clqSet;
   CliqueSet *clqSetEx;
} CBData;

static char probName[256];
static FILE *flogdb = NULL;  /* log dual bound */
static time_t startLogdb;
static int useCliqueCuts = 1;
static int useOtherCuts = 0;
static int timeLimit = 300;
static double bestBound = 0.0;
static int nViolatedCliquesFound = 0;
static int maxDepth;
static int minViol;
static int maxPasses;
static int extendCliques = 1;

/* parameters:
   instance useClique useOther timeLimit extendCliques */

int main ( int argc, char **argv )
{
   if ( argc < 9 )
   {
      fprintf( stderr, "Parameters: instance useClique(0,1) useOtherCuts(0,1) timeLimit maxDepth minViol maxPasses extendCliques\n" );
      exit(1);
   }

   useCliqueCuts = atoi( argv[2] );
   useOtherCuts = atoi( argv[3] );
   timeLimit = atoi( argv[4] );
   maxDepth = atoi( argv[5] );
   minViol = atoi( argv[6] );
   maxPasses = atoi( argv[7] );
   extendCliques = atoi( argv[8] );

   printf("Parameters:\n");
   printf("\tclique %d\n", useCliqueCuts );
   printf("\tother %d\n", useOtherCuts );
   printf("\ttimeLimit %d\n", timeLimit );
   printf("\tmaxDepth %d\n", maxDepth );
   printf("\tminViol %d\n", minViol );
   printf("\tmaxPasses %d\n", maxPasses );
   printf("\textendCliques %d\n", extendCliques );

   getFileName( probName, argv[1] );

   glp_prob *lp = glp_create_prob();
   if (glp_read_lp( lp, NULL, argv[1] ))
   {
      fprintf( stderr, "Error reading file: %s\n", argv[1] );
      exit( EXIT_FAILURE );
   }

   glp_iocp ioParams;
   glp_init_iocp( &ioParams );

   CGraph *cgraph = cgraph_create( glp_get_num_cols( lp ) );

   LPCallbacks lpc;
   cd_ini_lp_callbacks( &lpc );
   ConflictDiscover *cd = cd_create( cgraph, lp, &lpc );

   CBData cbData;

   cbData.cgraph = cgraph;
   cbData.cd = cd;
   cbData.idxExtension = xmalloc(sizeof(int)*(glp_get_num_cols( lp )+1));
   cbData.newNodesCap = ((glp_get_num_cols( lp )+1)*2);
   cbData.newNodes = xmalloc(sizeof(int)*(cbData.newNodesCap));

   cbData.idx = NULL;
   cbData.idxCap = 0;

   cbData.callsAtDepth = xmalloc(sizeof(int)*glp_get_num_cols( lp ));
   cbData.lastDepth = -1;
   cbData.lastNode = -1;
   cbData.clqe = clqe_create();
   cbData.clqSet = clq_set_create();
   cbData.clqSetEx = clq_set_create();
   cbData.ivExtended = NULL;
   cbData.ivExtendedCap = 0;

   memset( cbData.callsAtDepth, 0, sizeof(int)*glp_get_num_cols( lp ) );

   ioParams.cb_info  = &cbData;
   ioParams.cb_func  = cut_callback;
   ioParams.fp_heur  = GLP_OFF;
   ioParams.bt_tech  = GLP_BT_BLB;

   if (useOtherCuts)
   {
      ioParams.mir_cuts = GLP_ON;
      ioParams.gmi_cuts = GLP_ON;
      ioParams.cov_cuts = GLP_ON;
   }
   else
   {
      ioParams.mir_cuts = GLP_OFF;
      ioParams.gmi_cuts = GLP_OFF;
      ioParams.cov_cuts = GLP_OFF;
   }
   ioParams.tm_lim   = timeLimit*1000;

   int result = 0;

   int status = glp_simplex( lp, NULL );
   if (status)
   {
      fprintf( stderr, "Could not solve linear relaxation.\n");
      result = 1;
      goto TERMINATE;
   }
   else
   {
      switch ( glp_get_status( lp ) )
      {
         case GLP_INFEAS:
         case GLP_NOFEAS:
         case GLP_UNBND:
            fprintf( stderr, "Linear Program Inseasible or Unbounded.\n" );
            result = 1;
            goto TERMINATE;
            break;
         default:
            break;
      }
   }

   /* log generation */
   char fName[256];
   char cCuts[32] = "-";
   if (useCliqueCuts)
      sprintf( cCuts, "cCuts" );
   char oCuts[32] = "-";
   if (useOtherCuts)
      sprintf( oCuts, "oCuts" );
   sprintf( fName, "%s-log-dbound-%s-%s-%d-%d-%d.txt", probName, cCuts, oCuts, maxDepth, minViol, maxPasses );
   flogdb = fopen( fName, "w" );
   if (!flogdb)
   {
      fprintf( stderr, "Could not open file %s.\n", fName );
      exit( EXIT_FAILURE );
   }
   startLogdb = time( NULL );
   fprintf( flogdb, "%d %.6f 0\n", 0, glp_get_obj_val(lp) ); fflush( flogdb );


   status = glp_intopt( lp, &(ioParams) );
   switch (status)
   {
      case GLP_EFAIL :
      case GLP_EROOT :
      case GLP_EBOUND :
         fprintf( stderr, "Error Solving MIP.");
         result = 1;
         goto TERMINATE;
   }

   fprintf( flogdb, "%d %.6f %d\n", timeLimit, bestBound, nViolatedCliquesFound ); fflush( flogdb );


TERMINATE:
   glp_delete_prob( lp );

   if ( cbData.cgraph )
      cgraph_free( &cbData.cgraph );

   if ( cbData.cd )
      cd_free( &(cbData.cd) );

   if ( cbData.idx )
      free( (cbData.idx) );

   if ( cbData.idxExtension )
      free( (cbData.idxExtension) );

   if ( cbData.newNodes )
      free( (cbData.newNodes) );

   if ( cbData.callsAtDepth )
      free( (cbData.callsAtDepth) );

   if (cbData.clqe)
      clqe_free( &(cbData.clqe) );

   if (cbData.clqSet)
      clq_set_free( &(cbData.clqSet) );

   if (cbData.clqSetEx)
      clq_set_free( &(cbData.clqSetEx) );

   if (cbData.ivExtended)
      free(cbData.ivExtended);

   if ( flogdb )
      fclose( flogdb );

   return result;
}

void cut_callback( glp_tree *tree, void *info )
{
   glp_prob *lp = NULL;

   CBData *cbd = info;

   CGraph *cgraph = cbd->cgraph;
   ConflictDiscover *cd = cbd->cd;

   int depth = 0;
   int node  = glp_ios_curr_node( tree );
   if (node != 0 )
      depth = glp_ios_node_level( tree, node );

   switch ( glp_ios_reason(tree) )
   {
      case GLP_ICUTGEN :
         {

            time_t currTime = time( NULL );
            int seconds = (int) difftime( currTime, startLogdb );
            static int lastSec = -1;
            static double lastBound = -9999999.0;
            if ( ( seconds != lastSec ) && ( seconds % 3 == 0)  )
            {
               lastSec = seconds;
               const int bestNode = glp_ios_best_node( tree );
               if (bestNode>0)
               {
                  bestBound = glp_ios_node_bound( tree, bestNode );
                  if ( fabs( bestBound-lastBound ) > 1e-7 )
                  {
                     fprintf( flogdb, "%d %.6f %d\n", seconds, bestBound, nViolatedCliquesFound ); fflush( flogdb );
                     lastBound = bestBound;
                  }
               }
            }

            if (!useCliqueCuts)
               break;

            if ((depth != cbd->lastDepth) || (node != cbd->lastNode))
            {
               cbd->pass = 1;
               cbd->lastDepth = depth;
               cbd->lastNode = node;
            }
            else
               cbd->pass++;

            if ( ( cbd->pass > maxPasses ) || ( depth > maxDepth ) )
               break;

            lp = glp_ios_get_prob( tree );

            assert( (lp) && (cgraph) );
            clock_t start = clock();

            cd_update_graph( cd, lp );

            vmg_adjust_vector_capacity( (void**) &(cbd->idx), &(cbd->idxCap), (glp_get_num_cols(lp)+1), sizeof(int) );

            int *nodeIdx = cbd->idx;
            /* filling fractional solution */

            int i;
            for ( i=0 ; (i<glp_get_num_cols(lp)) ; ++i )
               cgraph_set_node_weight( cgraph, i, cgraph_weight( glp_get_col_prim( lp, i+1 ) ) );

            const int *w = cgraph_get_node_weights( cgraph );

            /* selecting nodes to be considered for separation */
            int ni = 0;
            for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
               if (((cgraph_degree(cgraph,i)<=2)||(w[i]>=980)||(w[i]<20)))  /* only fractional variables v with d(v)>2 in cgraph */
                  nodeIdx[i] = -1;
               else
                  nodeIdx[i] = ni++;

            if ( ni < 3 )
               break;

            CGraph *ppcg = cgraph_create_induced_subgraph( cgraph, nodeIdx );

            Grasp *grasp = grasp_create( ppcg, minViol );

            CliqueExtender *clqe = cbd->clqe;
            CliqueSet *clqSet = cbd->clqSet;
            CliqueSet *clqSetEx = cbd->clqSetEx;

            grasp_set_alpha_choice( grasp, GRASP_ALPHA_REACTIVE );
            grasp_set_max_no_improvement( grasp, 8192 );
            grasp_run( grasp );

            const CliqueSet *clqSetPP = grasp_solution_set( grasp );

            nViolatedCliquesFound += clq_set_number_of_cliques(clqSetPP);

#ifdef DEBUG
            /*printf( "PP Nodes\n" ); fflush( stdout );
            clq_set_print( clqSetPP );*/
#endif

            /* going to original indexes */
            clq_set_clear( clqSet );
            clq_set_clear( clqSetEx );
            clq_set_add_using_original_indexes( clqSet, clqSetPP, cgraph_get_original_node_indexes(ppcg) );
#ifdef DEBUG
            /*printf( "Orig Nodes\n" ); fflush( stdout );
            clq_set_print( clqSet );*/
#endif

            vmg_adjust_vector_capacity( (void**) &(cbd->ivExtended), &(cbd->ivExtendedCap) , clq_set_number_of_cliques( clqSet ), sizeof(int) );
            int *extended = cbd->ivExtended;
            memset( extended, 0, sizeof(int)*clq_set_number_of_cliques( clqSet ) ); /* no one extended by now */
            if ( extendCliques )
            {
               /* filling priorities based in the reduced cost, those with
                  the smallest reduced cost will have higher priorities */
               {
                  int max = -INT_MAX, i;
                  for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
                  {
                     cbd->idx[i] = (glp_get_col_dual( lp, i+1 ) * 1000.0);
                     if ( cbd->idx[i] > max )
                        max = cbd->idx[i];
                  }
                  for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
                     cbd->idx[i] = max - cbd->idx[i];
               }
               clqe_set_priorities( clqe, cgraph_size(cgraph), cbd->idx );
               /* extending cliques */
               int i;
               for ( i=0 ; (i<clq_set_number_of_cliques(clqSet)) ; ++i )
                  if (clqe_extend( clqe, cgraph, clq_set_get_clique( clqSet, i ), CLQEM_RANDOM ))
                     extended[i] = 1;

               clq_set_clear( clqSetEx );
               clq_set_add_cliques( clqSetEx, clqe_get_cliques(clqe) );
            }
#ifdef DEBUG
/*            printf( "Ext Nodes\n" ); fflush( stdout );
            clq_set_print( clqSetEx ); */
#endif

            /* those which were not extended will be added as they are */
            {
               int i;
               int *extended = cbd->ivExtended;
               for ( i=0 ; (i<clq_set_number_of_cliques(clqSet)); ++i )
               {
                  if (!extended[i])
                     clq_set_add( clqSetEx, clq_set_clique_size(clqSet,i), clq_set_clique_elements(clqSet,i), 0 );
               }
            }

            clock_t end = clock();
            timeCGRAPH += ((double)(end-start))/((double)CLOCKS_PER_SEC);
            /*printf( "Conflict graph now has %d columns (ppgraph %d). CGRAPH called %d times using %.3f seconds.\n", cgraph_size( cgraph ), cgraph_size(ppcg) , (++callsCGRAPH), timeCGRAPH );*/

            int minCSize = INT_MAX;
            int maxCSize = -1;
            { /* adding cuts */

               int i;
               for ( i=0 ; (i<clq_set_number_of_cliques(clqSetEx)) ; ++i )
               {
                  const int size = clq_set_clique_size( clqSetEx, i );
                  if ( size > maxCSize )
                     maxCSize = size;
                  if ( size < minCSize )
                     minCSize = size;
                  int *idx = cbd->idx;
                  memcpy( idx, clq_set_clique_elements(clqSetEx,i), sizeof(int)*clq_set_clique_size(clqSetEx,i) );
                  lp_update_idxs_to_glpk( idx, size );

                  glp_ios_add_row( tree, NULL, 0, 0, size, idx-1, lp_get_ones_vector(size)-1, GLP_UP, 1.0 );
               }
            }

            printf("n: %d d: %d p: %d - separation found %d violated clique inequalities violation: [%.3f,%.3f]. nzs: {%d...%d}\n", node, depth, cbd->pass, clq_set_number_of_cliques(clqSetPP), (((double)grasp_get_worst_weight(grasp))/1000.0), (((double)grasp_get_best_weight(grasp))/1000.0), minCSize, maxCSize );

            grasp_free( &grasp );
            cgraph_free( &ppcg );

            /*if (callsCGRAPH==1)
              glp_ios_terminate( tree );*/
         }
         break;
   }

   return;
}

