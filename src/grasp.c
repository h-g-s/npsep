#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <assert.h>
#include "grasp.h"
#include "memory.h"
#include "vint_set.h"
#include "clique.h"
#include "macros.h"


/* default parameters */
#define MAX_NO_IMPROVEMENT 4096
#define MAX_SECONDS 60
/* default alpha value */
#define ALPHA 0.9
/* update probabilities interval */
#define REACTIVE_INTERVAL 1024

#define RA_ROULLETE_SIZE 1024

const double alpha_choices[] =
   { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

/* reactive grasp
   to update chances of selecting a given alpha value, after some iterations:
   k is an alpha value
   av( k ) = score( k ) / count( k )
   Qk = ( av( k ) / Zbest )^theta
   Pa( k ) = Qk / sum_k( av(k) )  */

struct _Grasp
{
   const CGraph *cgraph;
   const int *w;

   int alpha_choice;
   double alpha;

   int iteration;

   /* reactive grasp related info */
   int *ra_score;    /* score: sumation of weights in iterations where cliqueSet is updated
                        this prevents the algorithm of using a value of alpha which generated
                        repeated cliques with high weights */
   int *ra_count;
   double *ra_q;
   double *ra_p;
   int *ra_roullete;
   int ra_roullete_size;
   int ra_alphaIndex;

   int max_ni_it;
   int max_seconds;
   int minW;

   int *candidates;
   int nCandidates;

   int *nodesLeft;
   int nNodesLeft;
   float *evalNodesLeft;

   int *conflicts;
   int conflictsCap;

   int *clique;
   int cliqueSize;
   int cliqueWeight;

   int bestWeight;
   int worstWeight;

   CliqueSet *clqSet;

   /* information to build candidate list */
   float minNodeW;
   float maxNodeW;
   float diffNodeW;
   float minDegree;
   float maxDegree;
   float difDegree;
};

int frac_part_weight( int w )
{
   return MAX( w - 500, 500 - w );
}

Grasp *grasp_create( const CGraph *_cgraph, const int _minW )
{
   Grasp *grasp = xmalloc( sizeof(Grasp) );

   grasp->cgraph = _cgraph;
   grasp->w = cgraph_get_node_weights( _cgraph );

   grasp->alpha_choice = GRASP_ALPHA_REACTIVE;
   grasp->alpha = ALPHA;
   grasp->max_ni_it = MAX_NO_IMPROVEMENT;
   grasp->max_seconds = MAX_SECONDS;

   grasp->ra_score = NULL;
   grasp->ra_count = NULL;
   grasp->ra_q = NULL;
   grasp->ra_p = NULL;
   grasp->ra_roullete = NULL;

   grasp->candidates = xmalloc( sizeof(int)*cgraph_size( grasp->cgraph ) );
   grasp->clique = xmalloc( sizeof(int)*cgraph_size( grasp->cgraph ) );
   grasp->nodesLeft = xmalloc( sizeof(int)*cgraph_size( grasp->cgraph ) );
   grasp->conflictsCap = cgraph_size( grasp->cgraph )*10;   /* when processing temporary space is needed */
   grasp->conflicts = xmalloc( sizeof(int)*grasp->conflictsCap );
   grasp->evalNodesLeft = xmalloc( sizeof(float)*cgraph_size( grasp->cgraph ) );

   grasp->minNodeW = INT_MAX;
   grasp->maxNodeW = -1;

   {
      int i;
      for ( i=0 ; (i<cgraph_size(_cgraph)) ; ++i )
      {
         grasp->minNodeW = MIN( grasp->minNodeW, grasp->w[i] );
         grasp->maxNodeW = MAX( grasp->maxNodeW, grasp->w[i] );
      }
   }
   grasp->diffNodeW = grasp->maxNodeW - grasp->minNodeW;

   grasp->minW = _minW;

   grasp->minDegree = cgraph_min_degree( grasp->cgraph );
   grasp->maxDegree = cgraph_max_degree( grasp->cgraph );

   fflush( stdout ); fflush( stderr );
   assert( grasp->maxDegree >= grasp->minDegree );

   grasp->difDegree = grasp->maxDegree - grasp->minDegree;

   grasp->clqSet = clq_set_create();

   return grasp;
}

void grasp_set_alpha( Grasp *grasp, const double _alpha )
{
   assert( _alpha >= 0.0 ); assert( _alpha <= 1.0 );

   grasp->alpha = _alpha;
}

void grasp_set_max_no_improvement( Grasp *grasp, const int _max_ni_it )
{
   grasp->max_ni_it = _max_ni_it;
}

void grasp_run( Grasp *grasp )
{
   grasp->ra_alphaIndex = -1;
   grasp->iteration = 1;
   grasp->bestWeight = -1;
   grasp->worstWeight = INT_MAX;
   clock_t start = clock();
   int seconds;

   if ( grasp->alpha_choice==GRASP_ALPHA_REACTIVE )
      grasp_ini_reactive_info( grasp );

START_ITERATION:
   {

      grasp_iteration( grasp );

      if ( grasp->cliqueWeight >= grasp->minW )
      {
         int updated = clq_set_add( grasp->clqSet, grasp->cliqueSize, grasp->clique, grasp->cliqueWeight );

         if ( grasp->alpha_choice == GRASP_ALPHA_REACTIVE )
         {
            if (updated)
            {
               grasp->ra_score[grasp->ra_alphaIndex] += grasp->cliqueWeight;
            }
            grasp->ra_count[grasp->ra_alphaIndex]++;
         }
      }

      clock_t currentClock = clock();
      seconds = ((double)(currentClock-start)) / ((double)CLOCKS_PER_SEC);

         /* printf( "iteration %d - weight %d alpha %g alpha index %d\n", grasp->iteration, grasp->cliqueWeight, grasp->alpha, grasp->ra_alphaIndex ); */


      if ( grasp->bestWeight < grasp->cliqueWeight )
         grasp->bestWeight = grasp->cliqueWeight;

      if ( grasp->worstWeight > grasp->cliqueWeight )
         grasp->worstWeight = grasp->cliqueWeight;


      /* finished iteration, if it is a reactive grasp, updating score and count */
      if ( ( grasp->alpha_choice==GRASP_ALPHA_REACTIVE ) && ( grasp->iteration%REACTIVE_INTERVAL==0 ) )
         grasp_recompute_reactive_probabilities( grasp);
   }
   if ( ( grasp->iteration<grasp->max_ni_it ) && (seconds<grasp->max_seconds) )
   {
      grasp->iteration++;
      goto START_ITERATION;
   }

   printf("grasp ran in %d seconds\n", seconds ); fflush( stdout );

   /*printf("best weight: %d\n", grasp->bestWeight);
   printf("%d maximal cliques were found with weight greater or equal to %d\n", clq_set_number_of_cliques(grasp->clqSet), grasp->minW );
   printf("total weight: %d\n", clq_set_weight_sum( grasp->clqSet ) );*/
}

void grasp_set_alpha_choice( Grasp *grasp, const int choice )
{
   grasp->alpha_choice = choice;
}

void grasp_iteration( Grasp *grasp )
{
   const CGraph *cgraph = grasp->cgraph;

   grasp_select_alpha( grasp );

   int nConflicts;
   int *conflicts = grasp->conflicts;

   grasp->cliqueSize = 0;
   grasp->cliqueWeight = 0;

   int nodeToEnter;

   grasp_fill_nodes_left( grasp );

   if ( grasp->nNodesLeft == 0 )
      return;

   assert( grasp->nNodesLeft );

INSERT_NODE_INTO_CLIQUE:

   grasp_build_candidate_list( grasp );

   assert( grasp->nCandidates );

   nodeToEnter = grasp_select_from_candidate_list( grasp );

   assert( (nodeToEnter >= 0) && (nodeToEnter < cgraph_size(cgraph)) );

   grasp->clique[ grasp->cliqueSize++ ] = nodeToEnter;
   grasp->cliqueWeight += grasp->w[nodeToEnter];

   nConflicts = cgraph_get_all_conflicting( cgraph, nodeToEnter, conflicts, grasp->conflictsCap );
   grasp_update_nodes_left( grasp, nConflicts, conflicts, nodeToEnter );

   if ( grasp->nNodesLeft )
      goto INSERT_NODE_INTO_CLIQUE;

/* here we must have a maximal clique */
#ifdef DEBUG
   int *newNodes = xmalloc( cgraph_size(cgraph)*100 );

   qsort( grasp->clique, grasp->cliqueSize, sizeof(int), vint_set_cmp_int );
   int nodesToAdd = cgraph_get_candidates_clique_insertion( (CGraph*)cgraph, grasp->cliqueSize, grasp->clique, newNodes, cgraph_size(cgraph)*100 );
   if (nodesToAdd)
   {
      fprintf( stderr, "ERROR: GRASP is generating non-maximal cliques.\n" );
      fprintf( stderr, "  clique being built: " );
      {
         int i;
         for ( i=0 ; (i<grasp->cliqueSize) ; ++i )
            fprintf( stderr, "%d ", grasp->clique[i] );
         fprintf( stderr, "\n" );
      }
      fprintf( stderr, "  possible nodes for insertion:  " );
      {
         int i;
         for ( i=0 ; (i<nodesToAdd) ; ++i )
            fprintf( stderr, "%d ", newNodes[i] );
         fprintf( stderr, "\n" );
      }
      exit( EXIT_FAILURE );
   }

   free( newNodes );
#endif
}

void grasp_build_candidate_list( Grasp *grasp )
{
   float minEval;
   float maxEval;
   minEval =   9999999.0;
   maxEval =  -9999999.0;

   {
      int i;
      for ( i=0 ; (i<grasp->nNodesLeft) ; ++i )
      {
         const int node = grasp->nodesLeft[i];

         /* node must have conflict with all nodes inserted */
         int j=0;
         for (  ; (j<grasp->cliqueSize) ; ++j )
         {
            const int nodeClique = grasp->clique[j];

            if ( !(cgraph_conflicting_nodes( grasp->cgraph, node, nodeClique )) )
               break;
         }
         if ( j < grasp->cliqueSize )
            continue;

         const float degree = cgraph_degree( grasp->cgraph, node );

         const float evalDegree = 1.0 - ( (grasp->maxDegree - degree) / (grasp->difDegree+1.0) );

         const float evalWeight = 1.0 - (((double)frac_part_weight(grasp->w[node])) / (500.0)) ;

         const float eval = evalDegree + evalWeight;

         grasp->evalNodesLeft[ i ] = eval;

         if ( eval < minEval )
            minEval = eval;

         if ( eval > maxEval )
            maxEval = eval;
      }
   }

   const float diffEval = maxEval - minEval;

   const float filterEval = minEval + ( ( 1.0 - grasp->alpha ) * diffEval ) - 1e-5;

   grasp->nCandidates = 0;
   {
      int i;
      for ( i=0 ; (i<grasp->nNodesLeft) ; ++i )
         if ( grasp->evalNodesLeft[ i ] >=  filterEval)
            grasp->candidates[ grasp->nCandidates++ ] = grasp->nodesLeft[i];
   }
}

int grasp_select_from_candidate_list( const Grasp *grasp )
{
   assert( grasp->nCandidates >= 0 ); assert( grasp->nCandidates <= cgraph_size( grasp->cgraph) );

   const int index = INT_RANDOM( grasp->nCandidates );
      
   assert( index >= 0 ); assert( index < grasp->nCandidates );

   return grasp->candidates[ index ];
}

void grasp_fill_nodes_left( Grasp *grasp )
{
   grasp->nNodesLeft = 0;
   {
      int i;
      for ( i=0 ; (i<cgraph_size(grasp->cgraph) ) ; ++i )
      {
         grasp->nodesLeft[grasp->nNodesLeft] = i;
         grasp->nNodesLeft++;
      }
   }
}

void grasp_update_nodes_left( Grasp *grasp, const int n, const int neighbors[], const int nodeToEnter )
{
   {
      int i;
      for ( i=0 ; (i<grasp->nNodesLeft) ; ++i )
      {
         if ( grasp->nodesLeft[i] == nodeToEnter )
            goto REMOVE;

         if (vint_set_int_find(  grasp->nodesLeft[i], n, (int*)neighbors ) )
            continue;

REMOVE:
         if ( i != grasp->nNodesLeft-1 )
         {
            {
               int t = grasp->nodesLeft[i];
               grasp->nodesLeft[i] = grasp->nodesLeft[grasp->nNodesLeft-1];
               grasp->nodesLeft[grasp->nNodesLeft-1] = t;
            }
            --i;
            grasp->nNodesLeft--;
         }
         else
         {
            grasp->nNodesLeft--;
            return;
         }
      }
   }
}

void grasp_ini_reactive_info( Grasp *grasp )
{
   const int nAlphas = ( sizeof(alpha_choices) / sizeof(double) );

   if (!grasp->ra_score)
      grasp->ra_score = xmalloc( nAlphas*sizeof(int) );
   if (!grasp->ra_count)
      grasp->ra_count = xmalloc( nAlphas*sizeof(int) );
   if (!grasp->ra_q)
      grasp->ra_q = xmalloc( nAlphas*sizeof(double) );
   if (!grasp->ra_p)
      grasp->ra_p = xmalloc( nAlphas*sizeof(double) );
   if (!grasp->ra_roullete)
      grasp->ra_roullete = xmalloc( RA_ROULLETE_SIZE*sizeof(int) );

   memset( grasp->ra_score, 0, sizeof(int)*nAlphas );
   memset( grasp->ra_count, 0, sizeof(int)*nAlphas );

   grasp_fill_initial_probabilities( grasp );
   grasp_prepare_roullete( grasp );
}

CliqueSet *grasp_solution_set( Grasp *grasp )
{
   return grasp->clqSet;
}

void grasp_select_alpha( Grasp *grasp )
{
   /* prepare value of alpha */
   switch ( grasp->alpha_choice )
   {
   case GRASP_ALPHA_RANDOM:
      {
         const double nChoices = (sizeof(alpha_choices) / sizeof(double));
         grasp->ra_alphaIndex = INT_RANDOM( nChoices );
         grasp->alpha = alpha_choices[grasp->ra_alphaIndex];
         printf( "alpha %g\n", grasp->alpha );
      }
      break;
   case GRASP_ALPHA_REACTIVE:
      {
         int pRoullete = INT_RANDOM( grasp->ra_roullete_size );
         grasp->ra_alphaIndex = grasp->ra_roullete[pRoullete];
         grasp->alpha = alpha_choices[grasp->ra_alphaIndex];
      }
      break;
   }
}

void grasp_recompute_reactive_probabilities( Grasp *grasp )
{
   const int nAlphas = ( sizeof(alpha_choices) / sizeof(double) );
   int k;
   double sumQ = 0;
   double av;
   for ( k=0 ; (k<nAlphas) ; ++k )
   {
      if ( grasp->ra_count[k] > 0 )
         av = ((double)grasp->ra_score[k]) / ((double)grasp->ra_count[k]);
      else
         av = 1e-5;
      grasp->ra_q[k] = ( av / ((double)grasp->bestWeight) );
      if ( grasp->ra_q[k] < 1e-5 )
         grasp->ra_q[k] = 1e-5;
      sumQ += grasp->ra_q[k];
   }
   for ( k=0 ; (k<nAlphas) ; ++k )
   {
      grasp->ra_p[k] = grasp->ra_q[k] / sumQ;
      grasp->ra_p[k] -= 1e-5;
      if ( grasp->ra_p[k] <= 0.0 )
         grasp->ra_p[k] = 1e-5;
   }

   grasp_prepare_roullete( grasp );

   memset( grasp->ra_score, 0, sizeof(int)*nAlphas );
   memset( grasp->ra_count, 0, sizeof(int)*nAlphas );
}

void grasp_prepare_roullete( Grasp *grasp )
{
   const int nAlphas = ( sizeof(alpha_choices) / sizeof(double) );

#ifdef DEBUG
   /* checking probabilities */
   {
      int k;
      double sumP = 0.0;
      for ( k=0 ; (k<nAlphas) ; ++k )
         sumP += grasp->ra_p[k];
      /*printf("sump %g\n", sumP);*/
      assert( (fabs(sumP-1.0) <= 0.01) );
   }
#endif

   /* filling roullete */
   grasp->ra_roullete_size = 0;
   int k;
   for ( k=0 ; (k<nAlphas) ; ++k )
   {
      int nPositions = (int)(grasp->ra_p[k] * ((double)RA_ROULLETE_SIZE));
      int i;
      for ( i=0 ; (i<nPositions) ; ++i )
         grasp->ra_roullete[grasp->ra_roullete_size++] = k;
   }
}

void grasp_fill_initial_probabilities( Grasp *grasp )
{
   const int nAlphas = ( sizeof(alpha_choices) / sizeof(double) );
   int k;
   double equalP = ( 1.0 / ((double)nAlphas)) - 1e-5;
   if ( equalP <= 0.0 )
      equalP = 1e-5;

   for ( k=0 ; (k<nAlphas) ; ++k )
      grasp->ra_p[k] = equalP;
}

int grasp_get_best_weight( Grasp *grasp )
{
   return grasp->bestWeight;
}

int grasp_get_worst_weight( Grasp *grasp )
{
   return grasp->worstWeight;
}

void grasp_free( Grasp **grasp )
{
   if ( (*grasp)->ra_score )
      free( (*grasp)->ra_score );
   if ( (*grasp)->ra_count )
      free( (*grasp)->ra_count );
   if ( (*grasp)->ra_q )
      free( (*grasp)->ra_q );
   if ( (*grasp)->ra_p )
      free( (*grasp)->ra_p );
   if ( (*grasp)->ra_roullete )
      free( (*grasp)->ra_roullete );

   free( (*grasp)->evalNodesLeft );
   free( (*grasp)->conflicts );
   free( (*grasp)->nodesLeft );
   free( (*grasp)->clique );
   free( (*grasp)->candidates );
   clq_set_free( &((*grasp)->clqSet) );
   free( (*grasp) );
   (*grasp) = NULL;
}

