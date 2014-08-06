#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include "clique_enum.h"
#include "vint_set.h"
#include "memory.h"
#include "vectormgm.h"
#include "cgraph.h"
#include "macros.h"


#define SLACK_NEIGHS 100

struct _CliqueEnumerator
{
   int minW;

   IntSet *dfsNodesLeft;
   IntSet *dfsNeighs;

   int depthCap;

   int degreeCap;

   CliqueSet *clqSet;

   int *clique;
   int cliqueCap;

   const CGraph *cgraph;
};


CliqueEnumerator *clq_enum_create( const int minWeight )
{
   CliqueEnumerator *clq_enum = xmalloc( sizeof(CliqueEnumerator) );

   clq_enum->minW = minWeight;

   clq_enum->dfsNodesLeft = NULL;
   clq_enum->dfsNeighs = NULL;

   clq_enum->depthCap  = 0;
   clq_enum->degreeCap  = 0;

   clq_enum->clqSet = clq_set_create();

   clq_enum->clique = NULL;
   clq_enum->cliqueCap = 0;

   return clq_enum;
}

void clq_enum_check_space( CliqueEnumerator *clq_enum, const CGraph *cgraph, const IntSet *nodesLeft )
{
   const int maxDepth = vint_set_size( nodesLeft ) + 1;
   int lastValidDepth = 0;

   if ( maxDepth > clq_enum->depthCap )
   {
      if ( clq_enum->depthCap==0 )
      {
         clq_enum->depthCap = maxDepth;
         clq_enum->dfsNodesLeft = xmalloc( sizeof(IntSet)*clq_enum->depthCap );
         clq_enum->dfsNeighs = xmalloc( sizeof(IntSet)*clq_enum->depthCap );
      }
      else
      {
         lastValidDepth = clq_enum->depthCap;
         clq_enum->depthCap *= 2;
         if ( clq_enum->depthCap < maxDepth )
            clq_enum->depthCap = maxDepth;

         clq_enum->dfsNodesLeft = xrealloc( clq_enum->dfsNodesLeft, sizeof(IntSet)*clq_enum->depthCap );
         clq_enum->dfsNeighs = xrealloc( clq_enum->dfsNeighs, sizeof(IntSet)*clq_enum->depthCap );
      }
   }
   else
      lastValidDepth = clq_enum->depthCap;

   {
      /* initializing IntSets */
      int i;
      for ( i=lastValidDepth ; (i<clq_enum->depthCap) ; ++i )
         vint_set_init( clq_enum->dfsNodesLeft+i );
      for ( i=lastValidDepth ; (i<clq_enum->depthCap) ; ++i )
         vint_set_init( clq_enum->dfsNeighs+i );

   }

   /* checking maximum degree in nodesLeft */
   int maxDegree = -1;

   {
      int i;
      const int nlSize = vint_set_size( nodesLeft );
      const int *nl = vint_set_get_elements( nodesLeft );
      for ( i=0 ; (i<nlSize) ; ++i )
         maxDegree = MAX( maxDegree, cgraph_degree( cgraph, nl[i] ) );
   }

   clq_enum->degreeCap = maxDegree;
   int i;
   for ( i=0 ; (i<clq_enum->depthCap) ; ++i )
      vint_set_check_capacity( clq_enum->dfsNodesLeft+i, clq_enum->degreeCap*SLACK_NEIGHS );

   for ( i=0 ; (i<clq_enum->depthCap) ; ++i )
      vint_set_check_capacity( clq_enum->dfsNeighs+i, clq_enum->degreeCap*SLACK_NEIGHS );

   vmg_adjust_vector_capacity( (void**) &(clq_enum->clique), &(clq_enum->cliqueCap), cgraph_size( cgraph ), sizeof(int) );
}


/* recursive procedure */
void clq_enum_enum( CliqueEnumerator *clq_enum,
                  const int size, int clique[], const int weight,
                  const IntSet *nodesLeft, const int depth )
{
   const int nlSize = vint_set_size( nodesLeft );

   /* no more nodes to add */
   if ( nlSize == 0 )
   {
      if ( weight >= clq_enum->minW )
      {
         clq_set_add( clq_enum->clqSet, size, clique, weight );
      }

      return;
   }

   /* procesing memory for this depth in the three */
   int *dfsNeighs = vint_set_force_elements_access( clq_enum->dfsNeighs + depth );
   int *dfsNodesLeft = vint_set_force_elements_access( clq_enum->dfsNodesLeft + depth );
#ifdef DEBUG
   if ( !dfsNeighs )
   {
      printf("ERROR: clique_enum.c:no space in dfsNeighs - depth %d - depthCap %d\n", depth, clq_enum->depthCap );
      printf("\tenumerating clique [ ");
      int i;
      for ( i=0 ; (i<size) ; ++i )
         printf("%d ", clique[i] );
      printf("]\n");
      printf("\tnodes left (%d).\n", vint_set_size(nodesLeft) );
      abort();
      exit( EXIT_FAILURE );
   }
   assert( dfsNodesLeft != NULL );
#endif
   IntSet *dfsNodesLeftIS = clq_enum->dfsNodesLeft + depth;

   /* considering for insertion each one of nodes left */
   const int *nl = vint_set_get_elements( nodesLeft );
   {
      int i;
      for ( i=0 ; (i<nlSize) ; ++i )
      {
         const int node = nl[i];


         clique[size] = node;
#ifdef DEBUG
         int nn1, nn2;
         if (!clq_validate( clq_enum->cgraph, size+1, clique, &nn1, &nn2 ))
         {
            fprintf( stderr, "Error: nodes %d and %d are not neighbors.\n", nn1, nn2 );
            exit( EXIT_FAILURE);
         }
#endif

         const int nNeighs = cgraph_get_all_conflicting( clq_enum->cgraph, node, dfsNeighs, clq_enum->degreeCap*SLACK_NEIGHS );
         const int newNodesLeft = vint_set_intersection( dfsNodesLeft, nNeighs, dfsNeighs, nodesLeft );
         qsort( dfsNodesLeft, newNodesLeft, sizeof(int), vint_set_cmp_int );

         vint_set_force_size( dfsNodesLeftIS, newNodesLeft );
#ifdef DEBUG
         vint_set_force_check( dfsNodesLeftIS );
#endif

         const CGraph *cgraph = clq_enum->cgraph;
         const int nodeW = cgraph_get_node_weight( cgraph, node );

         /* going deep in recursion */
         clq_enum_enum( clq_enum, size+1, clique, weight+nodeW, dfsNodesLeftIS, depth+1 );
      }
   }
}

void clq_enum_run( CliqueEnumerator *clq_enum, const CGraph *cgraph,
                  const int size, const int clique[], const int weight,
                  const IntSet *nodesLeft )
{
   //clock_t start = clock();

   clq_enum_check_space( clq_enum, cgraph, nodesLeft );
   clq_set_clear( clq_enum->clqSet );

   clq_enum->cgraph = cgraph;

   memcpy( clq_enum->clique, clique, sizeof(int)*size );

   /* cleaning contents */
   int i;
   for ( i=0 ; (i<clq_enum->depthCap) ; ++i )
      vint_set_clear( clq_enum->dfsNeighs+i );
   for ( i=0 ; (i<clq_enum->depthCap) ; ++i )
      vint_set_clear( clq_enum->dfsNodesLeft+i );

   clq_enum_enum( clq_enum, size, clq_enum->clique, weight, nodesLeft, 0 );

   //clock_t end = clock();
   //const double time = ((double)end-start)/((double)CLOCKS_PER_SEC);
   //printf("clique enumeration ran in %.3f seconds.\n", time);
}

const CliqueSet *clq_enum_get_cliques( const CliqueEnumerator *clq_enum )
{
   return clq_enum->clqSet;
}

void clq_enum_free( CliqueEnumerator **clq_enum )
{
   const int depthCap = (*clq_enum)->depthCap;

   if ( (*clq_enum)->dfsNodesLeft )
   {
      int i;
      for ( i=0 ; (i<depthCap) ; ++i )
         vint_set_clean( (*clq_enum)->dfsNodesLeft + i );

      free( (*clq_enum)->dfsNodesLeft );
   }

   if ( (*clq_enum)->dfsNeighs )
   {
      int i;
      for ( i=0 ; (i<depthCap) ; ++i )
         vint_set_clean( (*clq_enum)->dfsNeighs + i );

      free( (*clq_enum)->dfsNeighs );
   }

   if ( (*clq_enum)->clique )
      free( (*clq_enum)->clique );

   clq_set_free( &((*clq_enum)->clqSet) );

   free(*clq_enum);

   *clq_enum = NULL;
}

void clq_enum_set_min_weight( CliqueEnumerator *clqEnum, const int minW )
{
   clqEnum->minW = minW;
}

