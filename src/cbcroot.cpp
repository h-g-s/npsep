#include <cstdio>
#include <cstdlib>
#include <climits>
#include <vector>
#include <OsiClpSolverInterface.hpp>
#include <CglClique.hpp>
#include <CglEClique.hpp>
#include <CglOddHole.hpp>
#include <CglGomory.hpp>
#include <CglTwomir.hpp>
#include "osi_cgraph.h"

extern "C"
{
#include "strUtils.h"
}

using namespace std;

#define MIN_VIOLATION 0.02
#define MAX_TIME 300

/* CBC stuff */

void readLP( const char *fileName );

/* parameters */

typedef enum
{
   Default,   /* our new clique separation implementation */
   CglSepM
} SeparationMethod;

SeparationMethod sepMethod = Default;
char methodName[256];
int maxTwoMir = 0;
int maxGomory = 0;
string optFile;

void help()
{
   printf("usage:\n");
   printf("cbcroot problem [-cgl]\n");
   printf("\t-cgl changes separation routine to CGL\n");
   printf("\t-mtm=int max two mir passes for when no clique cuts are found\n");
   printf("\t-mgm=int max gomory passes for when no clique cuts are found\n");
   printf("\t-allint  considers that all variables are integers\n");
   clq_sep_params_help_cmd_line();
}

static bool transformInPI = false;

void parseParameters( int argc, char **argv );

static OsiSolverInterface *solver = NULL;

OsiClpSolverInterface *realSolver = new OsiClpSolverInterface();

/* Fills a vector with the variable names, including the binary complements */
vector<string> getVarNames(const vector<string> &colNames, int numCols);

int main( int argc, char **argv )
{
   if ( argc < 2 )
   {
      help();
      exit( EXIT_FAILURE );
   }

   realSolver = new OsiClpSolverInterface();
   /* makes CLP faster for hard instances */
   realSolver->getModelPtr()->setPerturbation(50);

   solver = (OsiSolverInterface*) realSolver;
   parseParameters( argc, argv );

   readLP( argv[1] );

   if (transformInPI)
   {
      vector< int > ints( solver->getNumCols() );
      for ( int i=0 ; (i<solver->getNumCols()) ; ++i ) ints[i] = i;
      solver->setInteger( &ints[0], solver->getNumCols() );
   }
 
   const int numCols = solver->getNumCols(), numRows = solver->getNumRows();
   char problemName[ 256 ];
   getFileName( problemName, argv[1] );

   CGraph *cgraph = osi_build_cgraph( solver );

   {
      CliqueSeparation *clqSep = clq_sep_create( cgraph );
      clq_sep_set_params_parse_cmd_line( clqSep, argc, (const char**)argv );
      clq_sep_free( &clqSep );
   }
   CliqueSeparation *clqSep = NULL;

   int pass = 0;
   int newCuts = 0, totalCuts = 0;
   vector<double> ones( numCols, 1.0 );
   const CliqueSet *clqSet = NULL;

   clock_t start = clock();
   solver->initialSolve();

   if (!solver->isProvenOptimal())
   {
      if (solver->isAbandoned())
      {
         fprintf( stderr, "LP solver abandoned due to numerical dificulties.\n" );
         exit( EXIT_FAILURE );
      }
      if (solver->isProvenPrimalInfeasible())
      {
         fprintf( stderr, "LP solver says PRIMAL INFEASIBLE.\n" );
         exit( EXIT_FAILURE );
      }
      if (solver->isProvenDualInfeasible())
      {
         fprintf( stderr, "LP solver says DUAL INFEASIBLE.\n" );
         exit( EXIT_FAILURE );
      }
      if (solver->isPrimalObjectiveLimitReached())
      {
         fprintf( stderr, "LP solver says isPrimalObjectiveLimitReached.\n" );
         exit( EXIT_FAILURE );
      }
      if (solver->isDualObjectiveLimitReached())
      {
         fprintf( stderr, "LP solver says isDualObjectiveLimitReached.\n" );
         exit( EXIT_FAILURE );
      }
      if (solver->isIterationLimitReached())
      {
         fprintf( stderr, "LP solver says isIterationLimitReached.\n" );
         exit( EXIT_FAILURE );
      }

      fprintf( stderr, "ERROR: Could not solve LP relaxation to optimality. Checking status...\n" );
      exit( EXIT_FAILURE );
   }

   double initialBound = solver->getObjValue();
   CliqueSeparation *cs = clq_sep_create( cgraph );
   clq_sep_set_params_parse_cmd_line( cs, argc, (const char**)argv );
   const int maxPasses = clq_sep_get_max_passes(cs);
   clq_sep_free( &cs );

   int passesGomory = 0;
   int passesTwoMir = 0;

   OsiCuts allCuts; //to avoid inserting cuts already inserted in previous iterations
   CoinAbsFltEq equal(1.0e-12);

   do
   {
      newCuts = 0;

      switch (sepMethod)
      {
         case Default :
            /* updating reduced costs */
            {
                CglEClique cliqueGen;
                OsiCuts cuts;
                CglTreeInfo info;
                info.level = 0;
                info.pass = 1;
                vector<string> varNames = getVarNames(solver->getColNames(), numCols);

                cliqueGen.parseParameters( argc, (const char**)argv );
                cliqueGen.setCGraph( cgraph );
                cliqueGen.setGenOddHoles( true ); //allow inserting odd hole cuts
                cliqueGen.colNames = &varNames;
                cliqueGen.generateCuts( *solver, cuts, info );
                newCuts = cuts.sizeCuts();
                solver->applyCuts( cuts );
            }


            break;
         case CglSepM :
            {
               CglClique cliqueGen;
               OsiCuts cuts;
               CglTreeInfo info;
               info.level = 0;
               info.pass = 1;
               cliqueGen.setMinViolation( MIN_VIOLATION );
               cliqueGen.setStarCliqueReport(false);
               cliqueGen.setRowCliqueReport(false);
               cliqueGen.generateCuts( *solver, cuts, info );
               newCuts = cuts.sizeCuts();

               solver->applyCuts( cuts );

            }
            break;
      }

      totalCuts += newCuts;

      clock_t pend = clock();
      double pTime = ((double)(pend-start)) / ((double)CLOCKS_PER_SEC);

      if ( pTime > MAX_TIME )
        newCuts = 0;

      if (newCuts<60)
      {
         if (passesGomory<maxGomory)
         {
            CglGomory cglGomory;
            OsiCuts cuts;
            CglTreeInfo info;
            info.level = 0;
            info.pass = 1;
            //            cglGomory.setMinViolation( MIN_VIOLATION );

            cglGomory.generateCuts( *solver, cuts, info );
            newCuts = cuts.sizeCuts();

            newCuts += cuts.sizeCuts();

            solver->applyCuts( cuts );
            ++passesGomory;
         }
         if (passesTwoMir<maxTwoMir)
         {
            CglTwomir cglTwoMir;
            OsiCuts cuts;
            CglTreeInfo info;
            info.level = 0;
            info.pass = 1;
            //            cglTwoMir.setMinViolation( MIN_VIOLATION );

            cglTwoMir.generateCuts( *solver, cuts, info );
            newCuts = cuts.sizeCuts();

            solver->applyCuts( cuts );
            newCuts += cuts.sizeCuts();

            ++passesTwoMir;
         }
      }

      ++pass;

      if (newCuts)
      {
         solver->resolve();
         if (!solver->isProvenOptimal())
         {
            if (solver->isAbandoned())
            {
               fprintf( stderr, "LP solver abandoned due to numerical dificulties.\n" );
               exit( EXIT_FAILURE );
            }
            if (solver->isProvenPrimalInfeasible())
            {
               fprintf( stderr, "LP solver says PRIMAL INFEASIBLE.\n" );
               exit( EXIT_FAILURE );
            }
            if (solver->isProvenDualInfeasible())
            {
               fprintf( stderr, "LP solver says DUAL INFEASIBLE.\n" );
               exit( EXIT_FAILURE );
            }
            if (solver->isPrimalObjectiveLimitReached())
            {
               fprintf( stderr, "LP solver says isPrimalObjectiveLimitReached.\n" );
               exit( EXIT_FAILURE );
            }
            if (solver->isDualObjectiveLimitReached())
            {
               fprintf( stderr, "LP solver says isDualObjectiveLimitReached.\n" );
               exit( EXIT_FAILURE );
            }
            if (solver->isIterationLimitReached())
            {
               fprintf( stderr, "LP solver says isIterationLimitReached.\n" );
               exit( EXIT_FAILURE );
            }

            fprintf( stderr, "ERROR: Could not solve LP relaxation. Exiting.\n" );
            exit( EXIT_FAILURE );
         }
      }
   }
   while ( (newCuts>0) && (pass<maxPasses) ) ;

   clock_t tend = clock();
   double totalTime = ((double)(tend-start)) / ((double)CLOCKS_PER_SEC);
   printf("%s %.7f %.7f %d %.3f\n", problemName, initialBound, solver->getObjValue(), totalCuts, totalTime);

   /*clq_sep_free( &clqSep );*/
   cgraph_free( &cgraph );

   //strcat(problemName, "CUT");
   //solver->writeMps(problemName);

   delete realSolver;

   return EXIT_SUCCESS;
}

void readLP( const char *fileName )
{
   solver->setIntParam(OsiNameDiscipline, 2);
   solver->messageHandler()->setLogLevel(1);
   solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);

   if ( strstr( fileName, ".lp" ) || strstr( fileName, ".LP" ) )
       solver->readLp( fileName );
    else
       solver->readMps( fileName );
}

void parseParameters( int argc, char **argv )
{
   int i;
   sprintf( methodName, "def" );
   for ( i=1 ; (i<argc) ; ++i )
   {
      if (strstr( argv[i], "-cgl"))
      {
         //printf("Using COIN CGL separation.\n");
         sprintf( methodName, "cgl" );
         sepMethod = CglSepM;
         continue;
      }

      if (strcmp( argv[i], "-allint" )==0)
          transformInPI = true;

      char pName[256];
      char pValue[256];
      getParamName( pName, argv[i] );
      getParamValue( pValue, argv[i] );


      if (strcmp( pName, "mtm" )==0)
      {
         maxTwoMir = atoi( pValue );
         printf("using at max %d two mir passes.\n", maxTwoMir );
         continue;
      }

      if (strcmp( pName, "mgm" )==0)
      {
         maxGomory = atoi( pValue );
         printf("using at max %d gomory passes.\n", maxGomory );
         continue;
      }

      if (strcmp( pName, "optFile" )==0)
      {
         optFile = pValue;
         continue;
      }

   }
}

vector<string> getVarNames(const vector<string> &colNames, int numCols)
{
   vector<string> varNames(numCols * 2);

   for(int i = 0; i < numCols; i++)
   {
      varNames[i] = colNames[i];
      varNames[i+numCols] = "¬" + colNames[i];
   }

   return varNames;
}
