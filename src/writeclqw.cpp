#include <cassert>
#include <cstdio>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <cstdlib>
#include <vector>
#include <OsiClpSolverInterface.hpp>
#include <CglClique.hpp>
#include "osi_cgraph.h"

extern "C"
{
#include "memory.h"
#include "strUtils.h"
#include "cgraph.h"
#include "conflict_discover.h"
#include "clique_separation.h"
}
#include "osi_cgraph.h"

using namespace std;

#define MIN_VIOLATION 0.02

/* minimum fractional part to a variable to be considered fractional */
#define MIN_FRAC      0.001

#define MAX( a, b ) ( (a)>(b) ? (a) : (b) )
#define MIN( a, b ) ( (a)<(b) ? (a) : (b) )

void readLP( const char *fileName, OsiSolverInterface *solver );

double fracPart( const double x );


int main( int argc, char **argv )
{
   OsiSolverInterface *solver = NULL;

   OsiClpSolverInterface *realSolver = new OsiClpSolverInterface();

   solver = (OsiSolverInterface*) realSolver;

   readLP( argv[1], solver );

   printf("writeclqw\n\n");
   char problemName[ 256 ];
   getFileName( problemName, argv[1] );
   printf("loaded %s \n", problemName );
   printf("\t%d variables (%d integer) %d rows %d nz\n\n", solver->getNumCols(), solver->getNumIntegers(), solver->getNumRows(), solver->getNumElements() );

   clock_t startOG = clock();
   CGraph *cgraph = osi_build_cgraph( solver );
   clock_t endOG = clock();
   const double timeOG = ((double)(endOG-startOG)) / ((double)CLOCKS_PER_SEC);
   
   cgraph_recompute_degree( cgraph );

   unsigned long int m = 0;
   {
      int i;
      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
         m+=cgraph_degree(cgraph,i);
   }

   /*
   printf("original cgraph dimensions: n: %d m: %lu t: %.3f\n", cgraph_size(cgraph), m, timeOG );

   {
      char oname[256];
      sprintf( oname, "%s_o.clwq", problemName );

      cgraph_save( cgraph, oname );
   }*/
#ifdef DEBUG
   // check new way to iterate through nodes
   {
      int i;
      int nNodesInCliques = 0;
      int *neighs = (int*) xmalloc( cgraph_size(cgraph)*1000 );
      int *costs = (int*) xmalloc( cgraph_size(cgraph)*sizeof(int) );
      for ( i=cgraph_size(cgraph)-1 ; (i>=0) ; --i )
         costs[i] = cgraph_size(cgraph) - i + 5;
      int nneighs = 0;
      NeighIterator *nit = nit_create();
      printf("checking...\n"); fflush(stdout);
      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      {
         printf("node %d degree %d \n" , i, cgraph_degree(cgraph, i)); fflush( stdout );  fflush( stderr );
         nneighs = cgraph_get_all_conflicting( cgraph, i, neighs, cgraph_size(cgraph)*1000 );
         nit_start( nit, cgraph, i, costs );
         int *p = neighs+nneighs-1;
         int n;
         int nel = 0;
         int j=0;

         /*
         printf("\n\n");
         fflush(stdout); fflush(stderr);
         printf("NEIGHS:");
         for ( j=0 ; (j<nneighs ) ; ++j )
            printf(" %d", neighs[j]);
         fflush( stdout );  fflush( stderr );
         printf("\n\nNEIG Its:");
 */

         p = neighs+nneighs-1;
         while ( (n=nit_next(nit))!=INT_MAX )
         {
            //printf(" %d", n);
   //         fflush( stdout );  fflush( stderr );
            //printf("[%d %d %d %d ]\n", j, nneighs-j, *p, n ); fflush(stdout); fflush(stderr);
            ++j;
            assert(  (p>=neighs)  );
            assert( p<neighs+nneighs );
            assert( *p == n );
            ++nel;

            --p;
            fflush( stdout );  fflush( stderr );
         }
         assert( nel == nneighs );
         //printf(" nneighs %d nel %d \n", nneighs, nel );  fflush(stdout);
      }

      nit_free( &nit );
      free( neighs );
      free( costs );
   }
#endif

   CliqueSeparation *clqSep = NULL;

   printf("solving relaxation ... ");
   clock_t start = clock();
   solver->initialSolve();
   printf("Initial dual bound %g\n", solver->getObjValue() );
   clock_t end = clock();
   printf("solved in %.3f\n", ((double)(end-start)) / ((double)CLOCKS_PER_SEC) );
   clock_t startSep, endSep;
   double timeSep;

   cgraph_print_summary( cgraph, "original conflict graph" );

   clock_t tend = clock();
   double totalTime = ((double)(tend-start)) / ((double)CLOCKS_PER_SEC);
   printf("\nEnd of root node relaxation. Dual limit: %g time: %.3f ", solver->getObjValue(), totalTime );

   {
       int i, idx = 0;
       const double *x = solver->getColSolution();

       int *iv = new int[ solver->getNumCols() ];
       for ( i=0 ; i<solver->getNumCols() ; ++i )
           iv[i] = -1;

       for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
           if ((cgraph_degree(cgraph, i)<2)||(fracPart(x[i])<MIN_FRAC))
               continue;
           else
               iv[i] = idx++;

       CGraph *newGraph = cgraph_create_induced_subgraph( cgraph, iv );

       int nRemoved = 0;
       for ( i=0 ; (i<cgraph_size(newGraph)) ; ++i )
       {
           if (cgraph_degree(newGraph, i)<2)
           {
               iv[cgraph_get_original_node_index(newGraph,i)] = -1;
               ++nRemoved;
           }
       }
       if (nRemoved)
       {
           /* updating iv */
           idx = 0;
           for ( i=0 ; (i<cgraph_size(cgraph)); ++i )
               if (iv[i] != -1 )
                   iv[i] = idx++;

           cgraph_free( &newGraph );
           newGraph = cgraph_create_induced_subgraph( cgraph, iv );
       }

       /* removing from this graph all nodes with degree <= 1 */


       for ( i=0 ; (i<cgraph_size(newGraph)) ; ++i )
           cgraph_set_node_weight( newGraph, i, cgraph_weight(x[cgraph_get_original_node_index(newGraph,i)]) );
       cgraph_free( &cgraph );
       cgraph = newGraph;

       delete[] iv;
   }

   printf("\n");
   cgraph_print_summary( cgraph, "pre-processed conflict graph" );

   if ( cgraph_size( cgraph ) == 0 )
   {
      printf("EMPTY conflict graph. exiting...\n");
      exit(0);
   }

   char fileName[256];
   getFileName( fileName, argv[1] );
   printf("\nFile name: %s\n", fileName );
   char outName[256];
   sprintf( outName, "%s.clqw", fileName );

   cgraph_save( cgraph, outName );

   cgraph_free( &cgraph );

   delete realSolver;

   return EXIT_SUCCESS;
}

void readLP( const char *fileName, OsiSolverInterface *solver )
{
   solver->setIntParam(OsiNameDiscipline, 2);
   solver->readMps( fileName );
   solver->setIntParam(OsiNameDiscipline, 2);
   solver->messageHandler()->setLogLevel(1);
   solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
}

double fracPart( const double x )
{
   double nextInteger = ceil( x );
   double previousInteger = floor( x );

   return MIN( nextInteger-x, x-previousInteger );
}
