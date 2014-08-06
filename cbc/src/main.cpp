#include <iostream>
#include <cstdio>
#include <ctime>
#include <coin/OsiClpSolverInterface.hpp>
#include "lp_callbacks.h"
#include "cbc_hooks.hpp"

extern "C"
{
#include "strUtils.h"
#include "cgraph.h"
#include "conflict_discover.h"
}
#include "MIPSep.hpp"

using namespace std;

/// validates and returns the cost of the clique
int validateClique( CGraph *cgraph, const set<int> *clique, const vector< int > *weights );

void saveFracSol( const char *fName, OsiSolverInterface *solver );

int main( int argc, char **argv )
{
   if ( argc<2 )
   {
      fprintf( stderr, "Enter instance name.\n" );
      exit( EXIT_FAILURE );
   }

   OsiClpSolverInterface realSolver;
   OsiSolverInterface *osiSolver = (OsiSolverInterface*) &realSolver;

   osiSolver->setIntParam(OsiNameDiscipline, 2);
   osiSolver->readLp( argv[1] );
   osiSolver->setIntParam(OsiNameDiscipline, 2);
   osiSolver->messageHandler()->setLogLevel(1);
   osiSolver->setHintParam(OsiDoReducePrint,true,OsiHintTry);

   printf("cbcnps: coin branch & cut node packing separation routines\n\n");
   char problemName[ 256 ];
   getFileName( problemName, argv[1] );
   printf("loaded %s \n", problemName );
   printf("\t%d variables (%d integer) %d rows\n\n", osiSolver->getNumCols(), osiSolver->getNumIntegers(), osiSolver->getNumRows() );

   printf("solving relaxation ... ");
   clock_t start = clock();
   osiSolver->initialSolve();
   clock_t end = clock();
   printf("solved in %.3f\n", ((double)(end-start)) / ((double)CLOCKS_PER_SEC) );

   CGraph *cgraph = NULL;
   cgraph = cgraph_create( osiSolver->getNumCols() );
   LPCallbacks lpc;
   lpc.cb_get_col_info = cbc_get_col_info;
   lpc.cb_get_row_info = cbc_get_row_info;
   lpc.cb_get_num_cols = cbc_get_num_cols;
   lpc.cb_get_num_rows = cbc_get_num_rows;
   lpc.cb_get_col_name = cbc_get_col_name;
   lpc.cb_get_row_name = cbc_get_row_name;
   ConflictDiscover *cd = cd_create( cgraph, osiSolver, &lpc );

   start = clock();
   printf("\nupdating conflict graph ... "); fflush( stdout );
   cd_update_graph( cd, osiSolver );
   end = clock();
   printf("done in %.3f seconds.\n", ((double)end-start)/((double)CLOCKS_PER_SEC) ); fflush( stdout );

   printf("conflicts: %d %d - %d\n", 15, 10181, cgraph_conflicting_nodes(cgraph,15,10181));
   printf("conflicts: %d %d - %d\n", 14, 49, cgraph_conflicting_nodes(cgraph,14,49));
   printf("conflicts: %d %d - %d\n", 0, 14, cgraph_conflicting_nodes(cgraph,0,14));

#ifdef DEBUG
   /*cgraph_check_node_cliques( cgraph );
   cgraph_check_neighs( cgraph ); */
#endif

   const double *x = osiSolver->getColSolution();
   int *wOrig = new int[ osiSolver->getNumCols() ];
   for ( int i=0 ; (i<osiSolver->getNumCols()) ; ++i )
   {
      wOrig[i] = cgraph_weight( x[i] );
      cgraph_set_node_weight( cgraph, i, wOrig[i] );
   }
   /// pre processed graph
   CGraph *ppcg = NULL;
   printf("before pre proc.\n"); fflush( stdout );
   vector< int > newIndexes; newIndexes.reserve( osiSolver->getNumCols() );
   int ni = 0;
   for ( int i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      if (((cgraph_degree(cgraph,i)<=2)||(wOrig[i]>=979)||(wOrig[i]<21)))
         newIndexes[i] = -1;
      else
         newIndexes[i] = ni++;

   ppcg = cgraph_create_induced_subgraph( cgraph, &(newIndexes[0]) );
#ifdef DEBUG
   printf("Induced Subgraph:\n");
   /*cgraph_check_node_cliques( ppcg );
   cgraph_check_neighs( ppcg );*/
#endif

   vector< int > oldIndexes; oldIndexes.resize( cgraph_size(ppcg) );
   vector< int > w; w.resize( ni );
   for ( int i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      if ( newIndexes[i] != -1 )
      {
         w[ newIndexes[i] ] = wOrig[i];
         oldIndexes[ newIndexes[i] ] = i;
      }

   printf(" nindexes %d %d - %d %d\n", 15, 10181, newIndexes[15], newIndexes[10181] );
   printf("conflicts: %d %d - %d\n", 2, 32, cgraph_conflicting_nodes(ppcg,2,32));

   printf("Pre-Processed graph has %d nodes.\n", cgraph_size(ppcg) );

   char ppgFileName[ 256 ], origFileName[ 256 ];
   sprintf( ppgFileName, "%s.clqw", problemName );
   sprintf( origFileName, "%s_orig.clqw", problemName );

#ifdef DEBUG
   {
      char sName[256];
      sprintf( sName, "%s_sol.txt", problemName );
      saveFracSol( sName, osiSolver );
   }
#endif

   cgraph_save( ppcg, ppgFileName );
   /*cgraph_save( cgraph, origFileName );*/
#ifdef DEBUG
   {
      clock_t start = clock();
      // checking pre processed graph
      printf("checking consistency of pre-processed graph ... "); fflush( stdout );
      for ( int i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      {
         const int newI = newIndexes[i];
         if ( newI == -1 )
            continue;
         for ( int j=i+1 ; (j<cgraph_size(cgraph)) ; ++j )
         {
            const int newJ = newIndexes[j];
            if ( newJ == -1 )
               continue;
            // original has conflict
            int hasCO = cgraph_conflicting_nodes( cgraph, i, j );
            int hasCP = cgraph_conflicting_nodes( ppcg, newI, newJ );

            if ( hasCO != hasCP )
            {
               fprintf( stderr, "nodes %d %d ( %d %d ) have different conflict in normal (%d) and pp graph (%d)\n", i,j,newI,newJ,hasCO,hasCP);
               fprintf( stderr, "names: %s %s\n", cgraph_get_node_name(cgraph, i), cgraph_get_node_name(cgraph, j) );
               exit( EXIT_FAILURE );
            }
         } // jsp100_500_50_1.clqw
      } // i
      clock_t end = clock();
      printf("done in %.3f seconds.\n" , ((double)(end-start))/((double)CLOCKS_PER_SEC) ); fflush( stdout );
   }
#endif


   /// saving pre processor information
   FILE *fp = fopen("pp.txt", "w");
   if (!fp)
   {
      fprintf( stderr, "Could not open file pp.txt\n" );
         exit( EXIT_FAILURE );
   }
   vector< string > ppNames;
   int j = 0;
   for ( int i=0 ; (i<osiSolver->getNumCols()) ; ++i )
   {
      if ( newIndexes[i] == -1 )
         continue;

      fprintf( fp, "%d %d %s %d\n", 1+j, 1+i, osiSolver->getColName(i).c_str(), w[j] );
      ++j;
      ppNames.push_back( osiSolver->getColName(i) );
   }
   fclose( fp );

   printf("\nStarting to solve MIP separation\n\n");
   MIPSep msep( ppcg, &(w[0]), ppNames, problemName );
   strcpy( msep.problemName, problemName );
   msep.writeLP();

   msep.solve();

   const set< set<int> > *cliques = msep.cliquesFound();
   printf("\n >> %zu maximal cliques found << \n", cliques->size() );

   /// saving cliques from MIP separation
   char cliqueFileName[ 256 ];
   sprintf( cliqueFileName, "%s_cliques.txt", problemName );
   FILE *f = fopen( cliqueFileName, "w" );
   if (!f)
   {
      fprintf( stderr, "Could not open file %s.\n", &(problemName[0]) );
         exit( EXIT_FAILURE );
   }
   set< set<int> >::const_iterator sIt = cliques->begin();
   for ( ; (sIt != cliques->end()) ; ++sIt )
   {
      const set<int> &clique = *sIt;

      int cliqueW = validateClique( ppcg, &clique, &w );

      fprintf( f, "[%d] ", cliqueW );
      for ( set< int >::const_iterator nodeIt=clique.begin() ; (nodeIt!=clique.end()) ; ++nodeIt )
         fprintf( f, "[%d %s]", (*nodeIt)+1, cgraph_get_node_name( ppcg, (*nodeIt)) );
      fprintf( f, "\n");
   }
   fclose( f );

   if (cgraph)
      cgraph_free( &cgraph );
   if (ppcg)
      cgraph_free( &ppcg );
   delete[] wOrig;
   cd_free( &cd );

   return EXIT_SUCCESS;
}


int validateClique( CGraph *cgraph, const set<int> *clique, const vector< int > *weights )
{
   for ( set< int >::const_iterator iIt=clique->begin() ; (iIt!=clique->end()) ; ++iIt )
   {
      set< int >::const_iterator jIt=iIt;
      jIt++;
      for (  ; (jIt!=clique->end()) ; ++jIt )
         if (!cgraph_conflicting_nodes( cgraph, *iIt, *jIt ))
         {
            fprintf( stderr, "\nERROR: Nodes %d and %d are not in conflict.\n", *iIt, *jIt );
            exit( EXIT_FAILURE );
         }
   }

   int result = 0;
   for ( set< int >::const_iterator iIt=clique->begin() ; (iIt!=clique->end()) ; ++iIt )
      result += weights->at( *iIt );

   return result;
}

void saveFracSol( const char *fName, OsiSolverInterface *solver )
{
   FILE *f = fopen( fName, "w" );

   char name[32];

   for ( int i=0 ; (i<solver->getNumCols()) ; ++i )
   {
      if ( solver->getColSolution()[i] < 1e-5 )
         continue;

      memset( name, ' ', 32 );
      name[31] = '\0';
      strncpy( name, solver->getColName(i).c_str(), 32 );
      int len = strlen(name);
      for ( int j=len ; j<31 ; ++j )
         name[j] = ' ';
      fprintf( f, "%s %.5f\n", name, solver->getColSolution()[i] );
   }

   fclose(f);
   f = NULL;
}

