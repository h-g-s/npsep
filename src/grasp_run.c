#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "memory.h"
#include "cgraph.h"
#include "grasp.h"

int main( int argc, char **argv )
{
   if (argc<2)
   {
      fprintf( stderr, "Enter problem name.\n" );
      return EXIT_FAILURE;
   }

   CGraph *cgraph = cgraph_load( argv[1]);
   const int *w = cgraph_get_node_weights( cgraph );
   CGraph *ppgraph = NULL;
   const int *ppW;

   const int origNodes = cgraph_size( cgraph );
   int *iv = xmalloc( sizeof(int)*cgraph_size(cgraph) );
   int *origIdx = xmalloc( sizeof(int)*cgraph_size(cgraph) );

   /* creating induced subgraph with interesting nodes */
   {
      memset( iv, 0, sizeof(int)*cgraph_size(cgraph) );
      int i, ni = 0;
      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
         if (((cgraph_degree(cgraph,i)<=2)||(w[i]>=979)||(w[i]<21)))
            iv[i] = -1;
         else
            iv[i] = ni++;
      ppgraph = cgraph_create_induced_subgraph( cgraph, iv );
      ppW = cgraph_get_node_weights( ppgraph );
#ifdef DEBUG
      {
         cgraph_check_node_cliques( ppgraph );
         cgraph_check_neighs( ppgraph );
         cgraph_check_preproc( ppgraph, cgraph );

         cgraph_save( ppgraph, "ppgraph.txt" );
         FILE *f = fopen( "graspgraph.txt", "w" );
         int i;
         for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
         {
            fprintf( f, "%d %s\n", iv[i], cgraph_get_node_name(cgraph,i) );
         }
         fclose(f);
      }
#endif
    }
   const int ppNodes = cgraph_size( ppgraph );

   {
      int i, j=0;
      for ( i=0 ; (i<origNodes) ; i++ )
         if (iv[i] != -1)
            origIdx[j++] = i;
   }

   printf("original nodes: %d pre-processed nodes: %d\n", origNodes, ppNodes );

   Grasp *grasp = grasp_create( ppgraph, 1010 );
   grasp_set_alpha_choice( grasp, GRASP_ALPHA_REACTIVE );
   grasp_run( grasp );

   CliqueSet *clqSet = grasp_solution_set( grasp );
   /* post processing here we
    * should go back to original indexes */
   clq_set_save( ppgraph, clqSet, "cliques.txt" );

   /* validating clique set */
   {
      int i,j;
      int n1, n2;
      int *t = xmalloc( sizeof(int)*cgraph_size(cgraph) );
      for ( i=0 ; (i<clq_set_number_of_cliques(clqSet)) ; ++i )
      {
         const int n = clq_set_clique_size( clqSet,i );
         const int *el = clq_set_clique_elements( clqSet,i );
         for ( j=0 ; (j<n) ; ++j )
         {
            assert( el[j] >= 0) ;
            assert( el[j] < ppNodes );
         }
         if (!clq_validate( ppgraph, n, el, &n1, &n2 ))
         {
            fprintf( stderr, "Error: nodes %d and %d are not neighbors.\n", n1, n2 );
            exit( EXIT_FAILURE );
         }
         memcpy( t, el, sizeof(int)*n );
         for ( j=0 ; (j<n) ; ++j )
         {
            assert( t[j] >= 0 );
            assert( t[j] < ppNodes );
            assert( origIdx[t[j]] >= 0 );
            assert( origIdx[t[j]] < origNodes );
            t[j] = origIdx[t[j]];
         }
         if (!clq_validate( cgraph, n, t, &n1, &n2 ))
         {
            fprintf( stderr, "Error: nodes %d and %d are not neighbors.\n", n1, n2 );
            exit( EXIT_FAILURE );
         }
      }
      free(t);
   }

   grasp_free( &grasp );
   cgraph_free( &cgraph );
   if (ppgraph)
      cgraph_free( &ppgraph );
   free(iv);
   free(origIdx);

   return EXIT_SUCCESS;
}
