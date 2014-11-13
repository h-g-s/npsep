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

using namespace std;

#define MIN_VIOLATION 0.02

/* minimum fractional part to a variable to be considered fractional */
#define MIN_FRAC      0.001

#define MAX( a, b ) ( (a)>(b) ? (a) : (b) )
#define MIN( a, b ) ( (a)<(b) ? (a) : (b) )

void readLP( const char *fileName, OsiSolverInterface *solver );

int main( int argc, char **argv )
{
    clock_t start;
    double pairAnalysis, cliqueAnalysis;

    OsiSolverInterface *solver = NULL;
    OsiClpSolverInterface *realSolver = new OsiClpSolverInterface();
    solver = (OsiSolverInterface*) realSolver;
    readLP( argv[1], solver );
    char problemName[ 256 ];
    getFileName( problemName, argv[1] );

    start = clock();
    CGraph *cgraphPair = osi_build_cgraph_pairwise( solver );
    pairAnalysis = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    if ( cgraph_size( cgraphPair ) == 0 )
    {
        printf("EMPTY conflict graph. exiting...\n");
        exit(0);
    }
    double conflictsPair = 0.0;
    for(int i = 0; i < cgraph_size( cgraphPair ); i++)
        conflictsPair += (double)cgraph_degree(cgraphPair, i);
    conflictsPair /= 2.0;
    //cgraph_save(cgraphPair, "pair");
    cgraph_free( &cgraphPair );

    start = clock();
    CGraph *cgraphClique = osi_build_cgraph( solver );
    cliqueAnalysis = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    if ( cgraph_size( cgraphClique ) == 0 )
    {
        printf("EMPTY conflict graph. exiting...\n");
        exit(0);
    }
    double conflictsClique = 0.0;
    for(int i = 0; i < cgraph_size( cgraphClique ); i++)
        conflictsClique += (double)cgraph_degree(cgraphClique, i);
    conflictsClique /= 2.0;
    //cgraph_save(cgraphClique, "clique");
    cgraph_free( &cgraphClique );

    printf("%s %.1lf %.3lf %.1lf %.3lf\n", problemName, conflictsPair, pairAnalysis, conflictsClique, cliqueAnalysis);

    delete realSolver;

    return EXIT_SUCCESS;
}

void readLP( const char *fileName, OsiSolverInterface *solver )
{
    solver->messageHandler()->setLogLevel(0);
    solver->setIntParam(OsiNameDiscipline, 2);
    solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    solver->readMps( fileName );
    //solver->readLp( fileName );
    solver->setIntParam(OsiNameDiscipline, 2);
}
