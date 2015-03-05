#include <cstdio>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <cstdlib>
#include <vector>
#include <OsiClpSolverInterface.hpp>
//#include <coin/OsiCpxSolverInterface.hpp>
#include <CglClique.hpp>
#include <CglEClique.hpp>
#include <CglOddHole.hpp>
#include <CglGomory.hpp>
#include <CglTwomir.hpp>
#include <stdio.h>
extern "C"
{
#include "strUtils.h"
#include "cgraph.h"
#include "memory.h"
#include "oddhs.h"
#include "clique_separation.h"
}
#include "osi_cgraph.h"
#include "strUtils.h"

using namespace std;

#define MIN_VIOLATION 0.02
#define MAX_TIME 300
#define EPS 1e-6

/* minimum fractional part to a variable to be considered fractional */
#define MIN_FRAC      0.001

double fracPart( const double x );

/* CBC stuff */

void readLP( const char *fileName );

void printGraphSummary( CGraph *cgraph, const char *graphName );

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

int addOddHoles( OsiSolverInterface *solver, OddHoleSep *oddhs );

OsiClpSolverInterface *realSolver = new OsiClpSolverInterface();

/* Fills a vector with the variable names, including the binary complements */
vector<string> getVarNames(const vector<string> &colNames, int numCols);

void decideLpMethod();

bool differentSense( const double v1, const double v2 )
{
	if ( (v1>1e-5) && (v2<-1e-5) )
	return true;

	if ( (v1<-1e-5) && (v2>1e-5) )
	return true;

	return false;

}

const double abs_mip_gap( const double v1, const double v2 )
{
    /* are equal */
    if ( fabs(v1-v2) <= EPS )
        return 0.0;

    /* have different sense */
    if (differentSense( v1, v2 ))
        return 1.0;

    /* one of them equals zero */
    if ( (fabs(v1)<=EPS) || (fabs(v2)<=EPS) )
        return 1.0;

    double minV, maxV;
    if (v1<v2)
    {
        minV = v1;
        maxV = v2;
    }
    else
    {
        minV = v2;
        maxV = v1;
    }

    if ( minV <= -0.0001 )
    {
        minV *= -1.0;
        maxV *= -1.0;

        std::swap( minV, maxV );
    }


    double result = 1.0-(minV/maxV);

    assert( result >= -0.0001 );
    assert( result <=  1.0001 );

    result = min( 1.0, result );
    result = max( 0.0, result );

    return result;
}

map<string, double> getOptimals()
{
	map<string, double> optimals;
	FILE *file = fopen(optFile.c_str(), "r");
    if(!file)
    {
        perror("Cant open this file!\n");
        exit(EXIT_FAILURE);
    }

    char line[128];
    if(fgets(line, 128, file))
    {
    	while(fgets(line, 128, file) != NULL)
    	{
    		char *instance, *cOpt;
    		instance = strtok(line, ",;\n");
    		cOpt = strtok(NULL, ",;\n");
    		optimals.insert(pair<string, double>(instance, atof(cOpt)));
    	}
    }

    fclose(file);
    return optimals;
}

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

   //decideLpMethod();

   //printf("cbcroot: root node relaxation and clique cuts\n\n");
   char problemName[ 256 ];
   getFileName( problemName, argv[1] );
   //printf("loaded %s \n", problemName );
   //printf("\t%d variables (%d integer) %d rows\n\n", numCols, solver->getNumIntegers(), numRows );


   map<string, double> optimals = getOptimals();
   assert(optimals.find(problemName) != optimals.end());
   double instOpt = optimals[problemName];

   CGraph *cgraph = osi_build_cgraph( solver );

   {
      CliqueSeparation *clqSep = clq_sep_create( cgraph );
      //clq_sep_set_verbose( clqSep, 1 );
      clq_sep_set_params_parse_cmd_line( clqSep, argc, (const char**)argv );
      //clq_sep_params_print( clqSep );
      clq_sep_free( &clqSep );
   }
   CliqueSeparation *clqSep = NULL;

   int pass = 0;
   int newCuts = 0, totalCuts = 0;
   vector<double> ones( numCols, 1.0 );
   const CliqueSet *clqSet = NULL;

   //printf(">>> solving relaxation ... ");
   clock_t start = clock();
   //solver->setHintParam( OsiDoDualInInitial, true, OsiHintDo );
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

   //printf("Initial dual bound %g\n", solver->getObjValue() );
   clock_t end = clock();
   //printf("solved in %.3f\n", ((double)(end-start)) / ((double)CLOCKS_PER_SEC) );
   printf("%.2lf %d %d %.7lf %.7lf %.7lf\n", ((double)(end-start)) / ((double)CLOCKS_PER_SEC), pass, 0, solver->getObjValue(),
   												instOpt, abs_mip_gap(solver->getObjValue(), instOpt));
   //clock_t startSep = 0, endSep;
   double timeSep;

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

      //const double *x = solver->getColSolution();

      clock_t startSep = clock();

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

                //removing duplicate cuts (inserted in previous iterations)
                bool cutsToRemove[cuts.sizeRowCuts()];
                fill(cutsToRemove, cutsToRemove + cuts.sizeRowCuts(), false);
                for(int i = 0; i < cuts.sizeRowCuts(); i++)
                {
                    const OsiRowCut &newOrc = cuts.rowCut(i);
                    for(int j = 0; j < allCuts.sizeRowCuts(); j++)
                    {
                        const OsiRowCut &orc = allCuts.rowCut(j);
                        if(newOrc == orc)
                        {
                            //printf("A duplicate cut has been detected!\n");
                            cutsToRemove[i] = true;
                            break;
                        }
                    }
                }
                for(int i = 0; i < cuts.sizeRowCuts(); i++)
                {
                    if(!cutsToRemove[i]) //if cut i was not inserted in allCuts we have to insert it
                        allCuts.insertIfNotDuplicate(cuts.rowCut(i), equal);
                    else cuts.eraseRowCut(i);
                }

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

               cliqueGen.generateCuts( *solver, cuts, info );
               newCuts = cuts.sizeCuts();

               solver->applyCuts( cuts );

            }

            break;
      }

      totalCuts += newCuts;

      clock_t pend = clock();
      double pTime = ((double)(pend-start)) / ((double)CLOCKS_PER_SEC);
      double sepTime = ((double)(pend-startSep)) / ((double)CLOCKS_PER_SEC);


      if ( pTime > MAX_TIME )
      {
         //printf("time limit reached. discarding cuts.\n");
         newCuts = 0;
      }
      //printf( "round %d : %s separated %d cuts in %.4f seconds. dual limit now: %g time: %.3f\n", pass+1, methodName, newCuts, sepTime, solver->getObjValue(), pTime );

      if (newCuts<60)
      {
         if (passesGomory<maxGomory)
         {
            printf("\n- generating gomory cuts.\n");
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
            printf("\n- generating two mir cuts.\n");
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
         //printf("resolving relaxation ... ");
         fflush( stdout );
         clock_t start = clock();
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
         clock_t end = clock();
         //printf("solved in %.3f after adding cuts of pass %d. dual limit now %g\n", ((double)(end-start)) / ((double)CLOCKS_PER_SEC), pass, solver->getObjValue() );
         //printf("\n");
         fflush( stdout );
         printf("%.2lf %d %d %.7lf %.7lf %.7lf\n", sepTime, pass, newCuts, solver->getObjValue(),
         								instOpt, abs_mip_gap(solver->getObjValue(), instOpt));
      }

      fflush( stdout );

   }
   while ( (newCuts>0) && (pass<maxPasses) ) ;

   clock_t tend = clock();
   double totalTime = ((double)(tend-start)) / ((double)CLOCKS_PER_SEC);
   //printf("\nend of root node relaxation. initial dual limit: %.7f final: %.7f time: %.3f total cuts: %d\n", initialBound, solver->getObjValue(), totalTime, totalCuts );

   //printf("min violation: %g\n", cliqueGen.getMinViolation() );

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

void printGraphSummary( CGraph *cgraph, const char *graphName )
{
   vector< int > neighs;
   neighs.reserve( cgraph_size(cgraph)*100 );
   printf("\n");
   printf("Summary for graph: %s\n", graphName );

   int minD = INT_MAX;
   int maxD = -1;
   double minF = DBL_MAX;
   double maxF = -1;
   double mostF = 0.0;
   const int numCol = solver->getNumCols();
   int nEdges = 0;
   for ( int i=0 ; (i<cgraph_size(cgraph)) ; ++i )
   {
      int nConf = cgraph_get_all_conflicting( cgraph, i, &(neighs[0]), numCol*100 );
      nEdges += nConf;

      if ( nConf > maxD )
         maxD = nConf;
      if ( nConf < minD )
         minD = nConf;

      double f = ((double)cgraph_get_node_weight(cgraph,i)) / 1000.0;

      if ( f < minF )
         minF = f;
      if ( f > maxF )
         maxF = f;

      //const double *lb = solver->getColLower();
      //const double *ub = solver->getColUpper();
      if ( solver->isBinary(cgraph_get_original_node_index(cgraph, i)) )
      {
         double fracPart = std::min( 1.0-f, f );
         if ( fracPart > mostF )
            mostF = fracPart;
      }

   }
   nEdges /= 2;

   printf("\tnodes : %d\n", cgraph_size(cgraph) );
   printf("\tedges : %d\n", nEdges );
   printf("\tminimum degree: %d\n", minD);
   printf("\tmaximum degree: %d\n", maxD);
   printf("\tminF: %.4f\n", minF);
   printf("\tmaxF: %.4f\n", maxF);
   printf("\tmostF: %.4f\n", mostF);
   printf("\n");
}

double fracPart( const double x )
{
   double nextInteger = ceil( x );
   double previousInteger = floor( x );

   return std::min( nextInteger-x, x-previousInteger );
}

void parseParameters( int argc, char **argv )
{
   int i;
   sprintf( methodName, "def" );
   for ( i=1 ; (i<argc) ; ++i )
   {
      if (strstr( argv[i], "-cgl"))
      {
         printf("Using COIN CGL separation.\n");
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

int addOddHoles( OsiSolverInterface *solver, OddHoleSep *oddhs )
{
   int r = 0;
   const int numCol = solver->getNumCols();
   vector< double > row( numCol, 1.0);
   OsiCuts oc;
   for ( int i=0; (i<oddhs_get_odd_hole_count(oddhs)) ; i++ )
   {
      const int *s = oddhs_get_odd_hole( oddhs, i );
      const int *e = oddhs_get_odd_hole( oddhs, i+1 );
      const int size = e-s;

      double lhs = 0.0;
      for ( int j=0; (j<size) ; j++ )
         lhs +=  solver->getColSolution()[s[j]];

      const double viol = lhs - ((double)(size/2));
      if (viol < 1e-5)
         continue;

      printf("lhs %g viol %g rhs %g\n", lhs, viol, (double)(size/2) );

      OsiRowCut orc;
      orc.setRow( size, s, &row[0] );
      orc.setUb( ((double)(size/2)) );
      orc.setLb( -COIN_DBL_MIN );
      orc.setGloballyValid( true );

      oc.insertIfNotDuplicate( orc );
      r++;
   }

   solver->applyCuts( oc );

   return r;
}

void decideLpMethod()
{
   if ( realSolver->getNumRows() > 100000 )
   {
      ClpSolve::SolveType method = ClpSolve::useBarrier;
      ClpSolve::PresolveType presolveType = ClpSolve::presolveOn;
      int numberPasses = 5;
#ifndef UFL_BARRIER
      int options[] = {0,0,0,0,0,0};
#else
      // we can use UFL code
      int options[] = {0,0,0,0,4,0};
#endif
      int extraInfo[] = {-1,-1,-1,-1,-1,-1};
      int independentOptions[] = {0,0,3};
      ClpSolve clpSolve(method,presolveType,numberPasses,
            options,extraInfo,independentOptions);
      // =======================
      // now pass options in
      realSolver->setSolveOptions(clpSolve);
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
