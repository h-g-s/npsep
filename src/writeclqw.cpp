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
    //clock_t start;
    //double cliqueAnalysis;

    OsiSolverInterface *solver = NULL;
    OsiClpSolverInterface *realSolver = new OsiClpSolverInterface();
    solver = (OsiSolverInterface*) realSolver;
    readLP( argv[1], solver );
    char problemName[ 256 ];
    getFileName( problemName, argv[1] );

    //start = clock();
    CGraph *cgraphClique = osi_build_cgraph( solver );
    //cliqueAnalysis = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    if ( cgraph_size( cgraphClique ) == 0 )
    {
        printf("EMPTY conflict graph. exiting...\n");
        exit(0);
    }
    unsigned long int conflictsClique = 0;
    for(int i = 0; i < cgraph_size( cgraphClique ); i++)
        conflictsClique += cgraph_degree(cgraphClique, i);
    conflictsClique /= 2;
    cgraph_free( &cgraphClique );

    //printf("%s %.1lu %.3lf\n", problemName, conflictsClique, cliqueAnalysis);

    delete realSolver;

    return EXIT_SUCCESS;
}

void readLP( const char *fileName, OsiSolverInterface *solver )
{
    solver->messageHandler()->setLogLevel(0);
    solver->setIntParam(OsiNameDiscipline, 2);
    solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    solver->setIntParam(OsiNameDiscipline, 2);
	
	if(strstr(fileName, ".mps") != NULL)
    	solver->readMps( fileName );
    else if(strstr(fileName, ".lp") != NULL)
    	solver->readLp( fileName );
    else
    {
    	perror("File not recognized!\n");
    	exit(EXIT_FAILURE);
    }
}