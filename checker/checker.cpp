extern "C"
{
#include "cgraph.h"
#include "memory.h"
#include "clique.h"
}

#include <cstdio>
#include <cstdlib>

using namespace std;

int main( int argc, char **argv )
{
   if ( argc < 3 )
   {
      fprintf( stderr, "graphFileName cliqueFileName.\n" );
      exit( EXIT_FAILURE );
   }


   int totalW = 0;

   CGraph *cgraph;
   cgraph = cgraph_load( argv[1] );

   const int *w;
   w = cgraph_get_node_weights( cgraph );

   CliqueSet *cset = clq_set_load( argv[2] );

   clq_set_print( cset );

   for ( int i=0 ; (i<clq_set_number_of_cliques(cset)) ; ++i )
   {
      printf("checking clique of weight %d ...", clq_set_weight(cset, i) );
      int size = clq_set_clique_size( cset, i );
      const int *el = clq_set_clique_elements( cset, i );
      int n1, n2;
      if (!clq_validate( cgraph, size, el, &n1, &n2 ))
      {
         printf("\n\tERROR: nodes %d and %d are not in conflict.\n", n1+1, n2+1 );
      }
      else
      {
         printf("OK ... maximal ... ");
         int cand[ size ];
/*         int nCand = clq_augmenting_candidates( cgraph, size, el, cand );*/
         int *candi = (int*) xmalloc( sizeof(int)*1000000 );
         int nCand = cgraph_get_candidates_clique_insertion( cgraph, size, el, cand, 1000000 );
         free(candi);
         if (nCand)
         {
            printf( "\nERROR: there ared %d candidates for clique augmentation: ", nCand );
            for ( int i=0 ; (i<nCand) ; ++i )
               printf("%d ", cand[i]+1)  ;
            printf("\n");
            printf("the clique is: \n");
            for ( int i=0 ; (i<size) ; ++i )
               printf("%d ", el[i]+1 );
            printf("\n");
         }
         else
            printf("YES. ");

         /// computing weight
         int cw = 0;
         for ( int i=0 ; (i<size) ; ++i )
            cw += w[ el[i] ];

         totalW += cw;
         printf("computed weight: %d\n", cw);
      }
   }

   printf("total weight: %d\n", totalW );


   cgraph_free( &cgraph );
   clq_set_free( &cset );

   return EXIT_SUCCESS;
}
