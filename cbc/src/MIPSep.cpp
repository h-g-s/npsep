#include <cstdio>
#include <ctime>
#include <coin/OsiClpSolverInterface.hpp>
#include <coin/CoinBuild.hpp>
#include <coin/CglPreProcess.hpp>
#include <coin/CbcModel.hpp>
#include <coin/CglClique.hpp>
#include <coin/CbcHeuristic.hpp>
#include <coin/CbcHeuristicFPump.hpp>
#include <coin/CbcHeuristicLocal.hpp>
#include <coin/CbcHeuristicRINS.hpp>
#include <coin/CbcHeuristicDiveGuided.hpp>
#include <coin/CbcHeuristicDINS.hpp>
#include <coin/CbcHeuristicFPump.hpp>
#include "mpIndex.hpp"
#include "MIPSep.hpp"
#include "MSepEvent.hpp"

using namespace std;

#define MAX_SECONDS 28800


// default prints during model creation
#define CALL_MODEL_CREATE_VAR( FUNCTION, DESCRIPTION, OBJTYPE, TOTALTIMEVAR )        \
  {                                                                              \
    clock_t start = clock();                                                     \
    const int maxDescrSize = 20;                                                 \
    char str[maxDescrSize];                                                      \
    memset( &(str[0]), ' ', sizeof(char)*maxDescrSize );                         \
    str[maxDescrSize-1] = '\0';                                                  \
    snprintf( &(str[0]), maxDescrSize, "\t%s", DESCRIPTION );                              \
    str[ strlen(str) ] = ' ';                                                    \
    printf( "%s ... ", &(str[0]) );                                              \
    fflush(stdout); fflush(stderr);                                              \
    int cols1 = lp->getNumCols();\
    FUNCTION;\
    int cols2 = lp->getNumCols();\
    clock_t end = clock();                                                       \
    const double seconds = ((double)end-start)/((double)CLOCKS_PER_SEC);         \
    printf("%d %s created in %.3f seconds\n", cols2-cols1, OBJTYPE, seconds );   \
    TOTALTIMEVAR += seconds;                                                     \
    fflush(stdout);                                                              \
  }

// default prints during model creation
#define CALL_MODEL_CREATE_CONS( FUNCTION, DESCRIPTION, OBJTYPE, TOTALTIMEVAR )        \
  {                                                                              \
    clock_t start = clock();                                                     \
    const int maxDescrSize = 20;                                                 \
    char str[maxDescrSize];                                                      \
    memset( &(str[0]), ' ', sizeof(char)*maxDescrSize );                         \
    str[maxDescrSize-1] = '\0';                                                  \
    snprintf( &(str[0]), maxDescrSize, "\t%s", DESCRIPTION );                              \
    str[ strlen(str) ] = ' ';                                                    \
    printf( "%s ... ", &(str[0]) );                                              \
    fflush(stdout); fflush(stderr);                                              \
    int cols1 = lp->getNumRows();\
    FUNCTION;\
    int cols2 = lp->getNumRows();\
    clock_t end = clock();                                                       \
    const double seconds = ((double)end-start)/((double)CLOCKS_PER_SEC);         \
    printf("%d %s created in %.3f seconds\n", cols2-cols1, OBJTYPE, seconds );   \
    TOTALTIMEVAR += seconds;                                                     \
    fflush(stdout);                                                              \
  }
      
MIPSep::MIPSep( CGraph *_cgraph, const int _w[], 
                const std::vector<std::string> &_colNames, const char *_probName ) :
   probName( _probName ),
   cgraph( _cgraph ),
   w( new int[ cgraph_size(_cgraph) ] ),
   clp( new OsiClpSolverInterface() ),
   colIndex( new MPIndex() )
{
   ///printf("%d\n", cgraph_size(cgraph) );

   memcpy( w, _w, sizeof(int)*cgraph_size(_cgraph) );

   colNames = _colNames;

   lp = (OsiSolverInterface*) clp;
   lp->setIntParam(OsiNameDiscipline, 2);
   lp->messageHandler()->setLogLevel(1);
   lp->setHintParam(OsiDoReducePrint,true,OsiHintTry);
   lp->setObjSense( -1.0 );
   lp->setObjName( "cliqueWeight" );

   double timeVars = 0.0;

   printf("Creating variables.\n");
   CALL_MODEL_CREATE_VAR( createXVars(), "x", "binary vars", timeVars );

   printf("Creating constraints.\n");
   CALL_MODEL_CREATE_CONS( createConsActivation(), "activation", "<= constraint", timeVars );
}

void MIPSep::createXVars()
{
   const int n = cgraph_size( cgraph );

   vector<double> obj; obj.reserve( n );
   
   for ( int i=0 ; (i<n) ; i++ )
      obj[i] = w[i];

   addBinVars( colNames, obj );
}

void MIPSep::solve()
{
   OsiSolverInterface *rootLP = lp;

   clock_t start = clock();
   printf("\nSolving LP relaxation ... "); fflush(stdout);
   rootLP->initialSolve();
   clock_t end = clock();
   timeRelax = ((double)(end-start))/((double)CLOCKS_PER_SEC);
   printf("done in %.2f seconds.\n\n", timeRelax );

   // preprocessing
   printf("Starting preprocessing ...\n");
   CglPreProcess preProc;
   int additionalPasses = MAX_SECONDS / 600;
   rootLP = preProc.preProcess( *lp, false, std::min(additionalPasses,5) );
   if (!rootLP)
   {
      fprintf( stderr, "PreProcessing says problem infeasible.\n" );
      return;
   }
   rootLP->resolve();

   /// mapping old columns to new ones
   map< int, int > ocno;
   const int *oCols = preProc.originalColumns();
   const int ncp = rootLP->getNumCols();
   for ( int i=0 ; ( i<ncp ) ; ++i )
      ocno[ oCols[i] ] = i;

   CbcModel cbc( *rootLP );

   MSepEvent mSepEvent( &cbc, cgraph, &cliques, oCols );
   cbc.passInEventHandler( &mSepEvent );
   ///MIPModelEventHandler mmEventHandler( &cbc, inst, &colIndex, V_X, fixations, activeEvents.size() );

   //CglClique cutClique;
   //cbc.addCutGenerator( &cutClique, -99, "Clique" );

   CbcRounding hRound( cbc );
   CbcHeuristicLocal hCombine( cbc );
   CbcHeuristicCrossover hCrossover( cbc );
   CbcHeuristicDINS hDins( cbc );
   CbcHeuristicRINS hRINS( cbc );

   CbcHeuristicFPump hfPump( cbc );
   hfPump.setIterationRatio( 1 );
   hfPump.setMaximumPasses( 50 );
   hfPump.setWhen( 13 );
   hfPump.setMaximumRetries( 1 );
   cbc.addHeuristic( &hfPump,  "fpump" );

   cbc.addHeuristic( &hRound,  "rounding" );
   cbc.addHeuristic( &hDins, "DINS" );
   cbc.addHeuristic( &hCrossover, "crossover" );
   cbc.addHeuristic( &hCombine,  "combine solutions" );
   cbc.addHeuristic( &hRINS,  "RINS" );

   //cbc.setNumberStrong(0);
   //cbc.setNumberBeforeTrust(0);

   cbc.setMaximumCutPassesAtRoot( 0 );
   cbc.setMaximumCutPasses( 0 );
   cbc.setMaximumSeconds( MAX_SECONDS );
   cbc.setMaximumSavedSolutions( 128 );

//   MIPHotStart::enter( lp->getSolution(), lp, &cbc, &preProc, rootLP );

   cbc.setNumberStrong( 0 );
   cbc.setNumberBeforeTrust( 0 );

   cbc.branchAndBound();

   if ( cbc.isAbandoned() )
   {
      fprintf( stderr, "ERROR: Solution process abandoned due to numerical dificulties.\n" );
      abort(); exit(1);
   }
   if ( (cbc.isProvenInfeasible()) || (cbc.isContinuousUnbounded()) )
   {
      fprintf( stderr, "ERROR: Solution process abandoned due to numerical dificulties.\n" );
      abort(); exit(1);
   }

   selectedNodes.clear();
}

void MIPSep::createConsActivation()
{
   const int n = cgraph_size( cgraph );
   const int nm1 = n-1;

   int idx[2];
   double coef[2] = { 1.0, 1.0 };

   CoinBuild cb;

   for ( int i=0 ; (i<nm1) ; ++i )
   {
      idx[0] = i;
      for ( int j=i+1 ; (j<n) ; ++j )
      {
         idx[1] = j;
         if (cgraph_conflicting_nodes( cgraph, i,j ))
            continue;

         //printf("creating for %d %d\n", i, j); fflush( stdout );
         cb.addRow( 2, idx, coef, -COIN_DBL_MAX, 1.0 );
      }
   }

   lp->addRows( cb );
}

void MIPSep::createYVars()
{
   #define MAX_COL_NAME 64
   const int n = cgraph_size( cgraph );
   const int nm1 = n-1;

   vector< string > names;
   vector< double > obj;

   const int expectedSize = n * 4;
   names.reserve( expectedSize );
   obj.reserve( expectedSize );

   for ( int i=0 ; (i<nm1) ; ++i )
   {
      for ( int j=i+1 ; (j<n) ; ++j )
      {
         if (!cgraph_conflicting_nodes( cgraph, i,j ))
            continue;
         char cname[MAX_COL_NAME];
         snprintf( cname, MAX_COL_NAME, "y_%s_%s", colNames[i].c_str(), colNames[j].c_str() );
         names.push_back( cname );
         obj.push_back( 0.0 );
      }
   }

   addBinVars( colNames, obj );
   #undef MAX_COL_NAME
}

void MIPSep::addBinVars( std::vector< std::string > &colNames, std::vector< double > &obj )
{
   const int firstIdx = lp->getNumCols();

   if ( colNames.size() == 0 )
      return;

   vector< double > lb( colNames.size(), 0.0 );
   vector< double > ub( colNames.size(), 1.0 );
   vector< int > colStarts( colNames.size()+1, 0 );

   int rows; double elements;

   lp->addCols( colNames.size(), &(colStarts[0]), &rows, &elements, &(lb[0]), &(ub[0]), &(obj[0]) );
   lp->setColNames( colNames, 0, colNames.size(), firstIdx );
   makeIntegers( firstIdx, colNames.size() );

   for ( size_t i=0 ; (i<colNames.size()) ; ++i )
      colIndex->associate( colNames[i].c_str(), firstIdx+i );
}

void MIPSep::makeIntegers( int colStart, int nCols )
{
   vector< int > integers;
   integers.reserve( nCols );
   for ( int i=0 ; (i<nCols) ; ++i )
      integers.push_back( colStart+i );

   lp->setInteger( &(integers[0]), nCols );
}

void MIPSep::saveSolutions( CbcModel *cbc, CglPreProcess *preProc )
{
   selectedNodes.reserve( cbc->getNumCols() );
   for ( int i=0 ; (i<cbc->numberSavedSolutions()) ; ++i )
   {
      const double *ppX = cbc->savedSolution( i );
      assert( ppX );

      selectedNodes.clear();

      if ( preProc )
      {
         const int *origCol = preProc->originalColumns();
         for ( int j=0 ; (j<cbc->getNumCols()) ; ++j )
            if ( ppX[j] >= 0.99 )
               selectedNodes.push_back( origCol[j] );
      }
      else
      {
         for ( int j=0 ; (j<cbc->getNumCols()) ; ++j )
            if ( ppX[j] >= 0.99 )
               selectedNodes.push_back( j );
      }

      char solName[ 256 ];
      snprintf( solName, 256, "%s_sol%03d.txt", probName.c_str(), i );

      FILE *f = fopen( solName, "w" );
      if (!f)
      {
         fprintf( stderr, "Could not open file %s.", solName );
         exit( EXIT_FAILURE );
      }

      fprintf( f, "      Cost : %g\n", cbc->savedSolutionObjective( i ) );
      if ( ( cbc->isProvenOptimal() ) && ( fabs(cbc->savedSolutionObjective( i ) - cbc->getObjValue()) < 0.001) )
         fprintf( f, "    Status : Optimal\n" );
      else
         fprintf( f, "    Status : Feasible\n" );

      fprintf( f, "Included nodes:\n" );

      if ( preProc )
      {
         const int *origCol = preProc->originalColumns();
         for ( int j=0 ; (j<cbc->getNumCols()) ; ++j )
            if ( ppX[j] >= 0.99 )
               fprintf( f, "%d\n", origCol[j] );
      }
      else
      {
         for ( int j=0 ; (j<cbc->getNumCols()) ; ++j )
            if ( ppX[j] >= 0.99 )
               fprintf( f, "%d\n", j );
      }

      fclose( f );
   }
}

void MIPSep::writeLP()
{
   lp->writeLp( problemName );
}

MIPSep::~MIPSep()
{
   delete clp;
   delete colIndex;
   delete[] w;
}
