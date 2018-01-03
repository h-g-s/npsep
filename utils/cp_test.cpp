#include <OsiClpSolverInterface.hpp>
#include "constraint_propagation.h"
#include "preprocess.h"

extern "C"
{
    #include "strUtils.h"
    #include "cgraph.h"
    #include "build_cgraph.h"
    #include "problem.h"
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
	OsiClpSolverInterface *realSolver = new OsiClpSolverInterface();
    realSolver->getModelPtr()->setPerturbation(50); /* makes CLP faster for hard instances */
    OsiSolverInterface *solver = (OsiSolverInterface*) realSolver;
	
	getFileName( problemName, argv[1] );
	readLP( argv[1], solver );


	printf("%s ", problemName);
	clock_t start = clock();
    Problem *p = problem_create_using_osi(solver);
    Preprocess *preproc = preprocess_create(p);
    Problem *preProcProblem = preprocess_basic_preprocessing(preproc);
 //    OsiSolverInterface* preProcSolver = problem_convert_to_osi(preProcProblem);
    
    clock_t end = clock();
    printf(" %.2lf\n", ((double)(end-start)) / ((double)CLOCKS_PER_SEC));
 //    char output[256];
 //    strncpy(output, problemName, 256);
 //    strcat(output, "_PP");
 //    preProcsolver->writeMps(output);
    
    problem_free(&preProcProblem);
    preprocess_free(&preproc);
    problem_free(&p);
    // delete preProcSolver;
 //    delete realSolver;

	// Problem *problem = problem_create_using_osi(solver);
	// CGraph *cgraph = build_cgraph(problem);

 //    Preprocess *preproc = preprocess_create(problem);
    //cpropagation_get_vars_to_fix(preproc, cgraph);

    // printf("%s %d\n", problemName, preprocess_get_num_gub_constraints(preproc));
    // for(int i = 0; i < preprocess_get_num_gub_constraints(preproc); i++)
    // {
    //     int rowSize = problem_row_size(problem, i);
    //     const int *idxs = problem_row_idxs(problem, i);
    //     const double *coefs = problem_row_coefs(problem, i);
    //     for(int j = 0; j < rowSize; i++)
    //         printf("%lf%s ", coefs[j], problem_var_name(problem, idxs[j]));
    //     printf("%c %lf\n", problem_row_sense(problem, i), problem_row_rhs(problem, i));
    // }


    // problem_free(&problem);
    // preprocess_free(&preproc);
    // cgraph_free(&cgraph);



	// clock_t start = clock();
	// CPropagation *cp = cpropagation_create(problem);
 //    int nindexes[solver->getNumCols()];
 //    cpropagation_get_vars_to_fix(cp, cgraph);
 //    double cpTime = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
 //    printf("toFix: %d\n", cpropagation_get_num_vars_to_fix(cp));
 //    int toFix = 0, activate = 0, deactivate = 0;
 //    for(int i = 0; i < solver->getNumCols(); i++)
 //    {
 //        if(cpropagation_var_is_to_fix(cp, i) == DEACTIVATE)
 //        {
 //            printf("1- %s\n", solver->getColName(i).c_str());
 //            toFix++;
 //            deactivate++;
 //            solver->setColBounds(i, 1.0, 1.0);
 //            solver->initialSolve();
 //            solver->branchAndBound();
 //            assert(solver->isBinary(i));
 //            if(solver->isProvenOptimal())
 //            	printf("%lf\n", solver->getObjValue());
 //            assert(!solver->isProvenOptimal());
 //            solver->setColBounds(i, 0.0, 1.0);
 //        }
        
 //        else if(cpropagation_var_is_to_fix(cp, i) == ACTIVATE)
 //        {
 //            toFix++;
 //            activate++;
 //            solver->setColBounds(i, 0.0, 0.0);
 //            solver->initialSolve();
 //            solver->branchAndBound();
 //            assert(solver->isBinary(i));
 //            assert(!solver->isProvenOptimal());
 //            solver->setColBounds(i, 0.0, 1.0);
 //        }
 //    }
 //    assert(toFix == cpropagation_get_num_vars_to_fix(cp));
 //    printf("%s %d %d %d %d %d %d %.2lf\n", problemName, solver->getNumCols(), solver->getNumRows(), solver->getNumElements(),
 //                                     toFix, activate, deactivate, cpTime);
    
    /* preprocessing and saving preprocessed lp */
    /* char output[256];
    sprintf(output, "%s_PP", problemName);
    OsiSolverInterface* preProcSolver = cpropagation_preprocess(cp, nindexes);
    preProcsolver->writeLp(output); */
    // cpropagation_free(cp);
    // problem_free(&problem);

    delete realSolver;

	return EXIT_SUCCESS;
}