#include <OsiClpSolverInterface.hpp>
#include <CglClique.hpp>
#include <CglEClique.hpp>
#include "osi_cgraph.h"
extern "C"
{
    #include "strUtils.h"
}

using namespace std;

#define MIN_VIOLATION 0.02
#define MAX_TIME 600
#define MAX_PASSES 999
#define EPS 1e-6

typedef enum
{
    Default,   /* our new clique separation implementation */
    CglSepM
} SeparationMethod;

SeparationMethod sepMethod = Default;
string optFile;
map<string, double> optimals;

void readLP(OsiSolverInterface *solver, const char *fileName)
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
    for ( i=1 ; (i<argc) ; ++i )
    {
        if (strstr( argv[i], "-cgl"))
        {
            sepMethod = CglSepM;
            continue;
        }

        char pName[256];
        char pValue[256];
        getParamName( pName, argv[i] );
        getParamValue( pValue, argv[i] );

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

/* Fills a vector with the variable names, including the binary complements */
vector<string> getVarNames(const vector<string> &colNames, int numCols);
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
        printf("Invalid number of parameters!\n");
        exit( EXIT_FAILURE );
    }

    OsiClpSolverInterface solver;
    solver.getModelPtr()->setPerturbation(50); /* makes CLP faster for hard instances */
    parseParameters( argc, argv );
    readLP( &solver, argv[1] );

    const int numCols = solver.getNumCols(), numRows = solver.getNumRows();
    char problemName[ 256 ];
    getFileName( problemName, argv[1] );

    CGraph *cgraph = osi_build_cgraph( &solver );
    int pass = 0, newCuts = 0, totalCuts = 0;
    double pTime, opt;

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

    clock_t start = clock();
    solver.initialSolve();

    if (!solver.isProvenOptimal())
    {
        if (solver.isAbandoned())
        {
            fprintf( stderr, "LP solver abandoned due to numerical dificulties.\n" );
            exit( EXIT_FAILURE );
        }
        if (solver.isProvenPrimalInfeasible())
        {
            fprintf( stderr, "LP solver says PRIMAL INFEASIBLE.\n" );
            exit( EXIT_FAILURE );
        }
        if (solver.isProvenDualInfeasible())
        {
            fprintf( stderr, "LP solver says DUAL INFEASIBLE.\n" );
            exit( EXIT_FAILURE );
        }
        if (solver.isPrimalObjectiveLimitReached())
        {
            fprintf( stderr, "LP solver says isPrimalObjectiveLimitReached.\n" );
            exit( EXIT_FAILURE );
        }
        if (solver.isDualObjectiveLimitReached())
        {
            fprintf( stderr, "LP solver says isDualObjectiveLimitReached.\n" );
            exit( EXIT_FAILURE );
        }
        if (solver.isIterationLimitReached())
        {
            fprintf( stderr, "LP solver says isIterationLimitReached.\n" );
            exit( EXIT_FAILURE );
        }

        fprintf( stderr, "ERROR: Could not solve LP relaxation to optimality. Checking status...\n" );
        exit( EXIT_FAILURE );
    }

    double initialBound = solver.getObjValue();
    printf("%.2lf %d %d %.7lf", ((double)(clock()-start)) / ((double)CLOCKS_PER_SEC), pass, 0, solver.getObjValue());
    if(!optFile.empty())
    {
        printf(" %.7lf %.7lf", opt, abs_mip_gap(solver.getObjValue(), opt));
    }
    printf("\n");

    do
    {
        clock_t startSep = clock();
        newCuts = 0;

        switch (sepMethod)
        {
            case Default :
            {
                CglEClique cliqueGen;
                OsiCuts cuts;
                CglTreeInfo info;
                info.level = 0;
                info.pass = 1;
                vector<string> varNames = getVarNames(solver.getColNames(), numCols);
                cliqueGen.parseParameters( argc, (const char**)argv );
                cliqueGen.setCGraph( cgraph );
                cliqueGen.setGenOddHoles( true ); //allow (or not) inserting odd hole cuts
                cliqueGen.colNames = &varNames;
                cliqueGen.generateCuts( solver, cuts, info );
                newCuts = cuts.sizeCuts();
                solver.applyCuts( cuts );
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
                cliqueGen.generateCuts( solver, cuts, info );
                newCuts = cuts.sizeCuts();
                solver.applyCuts( cuts );
            }
            break;
        }

        pTime = ((double)(clock()-start)) / ((double)CLOCKS_PER_SEC);
        if(pTime > MAX_TIME) break;

        totalCuts += newCuts;
        ++pass;

        if (newCuts)
        {
            solver.resolve();
            if (!solver.isProvenOptimal())
            {
                if (solver.isAbandoned())
                {
                    fprintf( stderr, "LP solver abandoned due to numerical dificulties.\n" );
                    exit( EXIT_FAILURE );
                }
                if (solver.isProvenPrimalInfeasible())
                {
                    fprintf( stderr, "LP solver says PRIMAL INFEASIBLE.\n" );
                    exit( EXIT_FAILURE );
                }
                if (solver.isProvenDualInfeasible())
                {
                    fprintf( stderr, "LP solver says DUAL INFEASIBLE.\n" );
                    exit( EXIT_FAILURE );
                }
                if (solver.isPrimalObjectiveLimitReached())
                {
                    fprintf( stderr, "LP solver says isPrimalObjectiveLimitReached.\n" );
                    exit( EXIT_FAILURE );
                }
                if (solver.isDualObjectiveLimitReached())
                {
                    fprintf( stderr, "LP solver says isDualObjectiveLimitReached.\n" );
                    exit( EXIT_FAILURE );
                }
                if (solver.isIterationLimitReached())
                {
                    fprintf( stderr, "LP solver says isIterationLimitReached.\n" );
                    exit( EXIT_FAILURE );
                }

                fprintf( stderr, "ERROR: Could not solve LP relaxation. Exiting.\n" );
                exit( EXIT_FAILURE );
            }

            pTime = ((double)(clock()-start)) / ((double)CLOCKS_PER_SEC);
            if(pTime > MAX_TIME) break;

            double sepTime = ((double)(clock()-startSep)) / ((double)CLOCKS_PER_SEC);
            printf("%.2lf %d %d %.7lf", sepTime, pass, newCuts, solver.getObjValue());
            if(!optFile.empty())
                printf(" %.7lf %.7lf", opt, abs_mip_gap(solver.getObjValue(), opt));
            printf("\n");
        }
    }
    while ( (newCuts>0) && (pass<MAX_PASSES) ) ;

    cgraph_free( &cgraph );

    return EXIT_SUCCESS;
}
