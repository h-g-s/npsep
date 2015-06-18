#include "constraint_propagation.h"
#include <OsiClpSolverInterface.hpp>
#include "osi_cgraph.h"

extern "C"
{
#include "strUtils.h"
}

using namespace std;

void readLP( const char *fileName, OsiSolverInterface *solver )
{
    solver->setIntParam(OsiNameDiscipline, 2);
    solver->messageHandler()->setLogLevel(0);
    solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    
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

int main( int argc, char **argv )
{
	char problemName[ 256 ];
	OsiSolverInterface *solver = new OsiClpSolverInterface();
	
	getFileName( problemName, argv[1] );
	readLP( argv[1], solver );

	CPropagation *cp = cpropagation_create(solver);
    int nindexes[solver->getNumCols()];
    OsiSolverInterface* preProcSolver = cpropagation_preProcess(cp, nindexes);

    /* saving preprocessed lp */
    char output[256];
    sprintf(output, "%s_PP", problemName);
    preProcSolver->writeLp(output);

    delete solver;
    cpropagation_free(cp);

	return EXIT_SUCCESS;
}