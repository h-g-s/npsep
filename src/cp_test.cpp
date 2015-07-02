#include <omp.h>
#include <OsiClpSolverInterface.hpp>
#include "osi_cgraph.h"
#include "constraint_propagation.h"

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

	double start = omp_get_wtime();
	CPropagation *cp = cpropagation_create(solver);
    int nindexes[solver->getNumCols()];
    cpropagation_get_vars_to_fix(cp);
    
    /* preprocessing and saving preprocessed lp */
    /* char output[256];
    sprintf(output, "%s_PP", problemName);
    OsiSolverInterface* preProcSolver = cpropagation_preprocess(cp, nindexes);
    preProcSolver->writeLp(output); */

    printf("%s %d %d %.2lf\n", problemName, cpropagation_get_num_vars_to_fix(cp), cpropagation_get_num_skipped(cp), omp_get_wtime() - start);

    delete solver;
    cpropagation_free(cp);

	return EXIT_SUCCESS;
}