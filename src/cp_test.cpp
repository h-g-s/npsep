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

int main( int argc, char **argv )
{
	char problemName[ 256 ];
	OsiSolverInterface *solver = new OsiClpSolverInterface();
	CPropagation *cpropagation;
	
	getFileName( problemName, argv[1] );
	readLP( argv[1], solver );
	cpropagation = cpropagation_create(solver);

	test(cpropagation);

    delete solver;
    cpropagation_free(cpropagation);

	return EXIT_SUCCESS;
}