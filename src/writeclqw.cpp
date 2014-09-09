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

int main( int argc, char **argv )
{
	clock_t start = clock(), end, test;
	double osiTime, recompTime, writeTime, totalTime, cpuTime;

    OsiSolverInterface *solver = NULL;

    OsiClpSolverInterface *realSolver = new OsiClpSolverInterface();

    solver = (OsiSolverInterface*) realSolver;

    readLP( argv[1], solver );

    //printf("writeclqw\n\n");
    char problemName[ 256 ];
    getFileName( problemName, argv[1] );
    /*printf("loaded %s \n", problemName );
    printf("\t%d variables (%d integer) %d rows %d nz\n\n", solver->getNumCols(), solver->getNumIntegers(), solver->getNumRows(), solver->getNumElements() );*/

    CGraph *cgraph = osi_build_cgraph( solver );
    osiTime = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    printf("osi_cgraph took %.3f seconds.\n", osiTime);

    printf("Recomputing degree...");
    test = clock();
    cgraph_recompute_degree( cgraph );
    recompTime = ((double(clock() - test))/((double)CLOCKS_PER_SEC));
    printf("Done in %.3f seconds\n", recompTime);
    
	//cgraph_print_summary( cgraph, "Conflict graph" );

    if ( cgraph_size( cgraph ) == 0 )
    {
        printf("EMPTY conflict graph. exiting...\n");
        exit(0);
    }

    printf("Writing conflict graph...");
    test = clock();
    char fileName[256];
    getFileName( fileName, argv[1] );
    //printf("\nFile name: %s\n", fileName );
    char outName[256];
    sprintf( outName, "%s.clqw", fileName );
    cgraph_save( cgraph, outName );
    writeTime = ((double(clock() - test))/((double)CLOCKS_PER_SEC));
    printf("Done in %.3f seconds\n", writeTime);

    end = clock();
    cpuTime = ((double(end - start))/((double)CLOCKS_PER_SEC));

    printf("%s %d %d %d %d %d\n",fileName, activePairwise, activeReadyToUse
    							, inactivePairwise, inactiveReadyToUse, mixedPairwise);

	printf("Total time: %.3f seconds\n", cpuTime);
	printf("Row: %s \t time: %.3f seconds\n", solver->getRowName(maxRowTimeIdx).c_str(), maxRowTime);

    if(!success) printf("Time limit exceeded!\n");

    cgraph_free( &cgraph );
    delete realSolver;

    return EXIT_SUCCESS;
}

void readLP( const char *fileName, OsiSolverInterface *solver )
{
    solver->setIntParam(OsiNameDiscipline, 2);
    solver->readMps( fileName );
    //solver->readLp( fileName );
    solver->setIntParam(OsiNameDiscipline, 2);
    solver->messageHandler()->setLogLevel(1);
    solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
}
