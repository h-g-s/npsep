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
#define MAX_TIME 600
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
map<string, double> optimals;

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

void getOptimals()
{
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
            double opt;
            if(strcmp(cOpt, "Infeasible")==0) opt = DBL_MAX;
            else opt = atof(cOpt);
            optimals.insert(pair<string, double>(instance, opt));
        }
    }

    fclose(file);
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

    double opt;
    if(!optFile.empty())
    {
        getOptimals();
        if(optimals.find(problemName) == optimals.end())
        {
        	fprintf(stderr, "ERROR: optimal value not found!\n");
        	exit(EXIT_FAILURE);
        }
        opt = optimals[problemName];
    }

    double initialBound = solver->getObjValue();
    clock_t end = clock();
    printf("%.2lf %d %d %.7lf", ((double)(end-start)) / ((double)CLOCKS_PER_SEC), pass, 0, solver->getObjValue());
    if(!optFile.empty())
    {
    	printf(" %.7lf %.7lf", opt, abs_mip_gap(solver->getObjValue(), opt));
    }
    printf("\n");
    clock_t startSep = 0, endSep;
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
        clock_t startSep = clock();

        switch (sepMethod)
        {
        case Default :
        {
            CglEClique cliqueGen;
            OsiCuts cuts;
            CglTreeInfo info;
            info.level = 0;
            info.pass = 1;
            vector<string> varNames = getVarNames(solver->getColNames(), numCols);

            cliqueGen.parseParameters( argc, (const char**)argv );
            cliqueGen.setCGraph( cgraph );
            cliqueGen.setGenOddHoles( true ); //allow (or not) inserting odd hole cuts
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
        double sepTime = ((double)(pend-startSep)) / ((double)CLOCKS_PER_SEC);


        if ( pTime > MAX_TIME )
            newCuts = 0;

        if (newCuts<60)
        {
            if (passesGomory<maxGomory)
            {
                //printf("\n- generating gomory cuts.\n");
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
                //printf("\n- generating two mir cuts.\n");
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
            fflush( stdout );
            printf("%.2lf %d %d %.7lf", sepTime, pass, newCuts, solver->getObjValue());
            if(!optFile.empty())
            	printf(" %.7lf %.7lf", opt, abs_mip_gap(solver->getObjValue(), opt));
            printf("\n");
        }

        fflush( stdout );

    }
    while ( (newCuts>0) && (pass<maxPasses) ) ;

    //clock_t tend = clock();
    //double totalTime = ((double)(tend-start)) / ((double)CLOCKS_PER_SEC);
    //printf("\nend of root node relaxation. initial dual limit: %.7f final: %.7f time: %.3f total cuts: %d\n", initialBound, solver->getObjValue(), totalTime, totalCuts );

    /*clq_sep_free( &clqSep );*/
    cgraph_free( &cgraph );
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
