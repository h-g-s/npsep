#include <cassert>
#include <cstdio>
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <vector>
#include <OsiClpSolverInterface.hpp>
#include "build_cgraph.h"
#include "problem.h"

extern "C"
{
#include "strUtils.h"
}

using namespace std;

void readLP( const char *fileName, OsiSolverInterface *solver );

int main( int argc, char **argv )
{
    clock_t start;
    double osiTime, problemTime, pairwiseTime;
    OsiClpSolverInterface solver;
    readLP( argv[1], &solver );
    char problemName[ 256 ];
    getFileName( problemName, argv[1] );

    /*-------------------------------------------------------------------------------------------------------------------*/
    start = clock();
    CGraph *cgraphOsi = build_cgraph( &solver );
    osiTime = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    if ( cgraph_size( cgraphOsi ) == 0 )
    {
        printf("EMPTY conflict graph. exiting...\n");
        exit(0);
    }
    unsigned long int conflictsOsi = 0;
    for(int i = 0; i < cgraph_size( cgraphOsi ); i++)
        conflictsOsi += cgraph_degree(cgraphOsi, i);
    conflictsOsi /= 2;
    /*-------------------------------------------------------------------------------------------------------------------*/
    start = clock();
    Problem *problem = problem_create_using_osi(&solver);
    CGraph *cgraphProblem = build_cgraph(problem);
    problemTime = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    if ( cgraph_size( cgraphProblem ) == 0 )
    {
        printf("EMPTY conflict graph. exiting...\n");
        exit(0);
    }
    unsigned long int conflictsProblem = 0;
    for(int i = 0; i < cgraph_size( cgraphProblem ); i++)
        conflictsProblem += cgraph_degree(cgraphProblem, i);
    conflictsProblem /= 2;
    /*-------------------------------------------------------------------------------------------------------------------*/
    start = clock();
    Problem *problem2 = problem_create_using_osi(&solver);
    CGraph *cgraphPairwise = build_cgraph_pairwise(problem2);
    pairwiseTime = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    if ( cgraph_size( cgraphPairwise ) == 0 )
    {
        printf("EMPTY conflict graph. exiting...\n");
        exit(0);
    }
    unsigned long int conflictsPairwise = 0;
    for(int i = 0; i < cgraph_size( cgraphPairwise ); i++)
        conflictsPairwise += cgraph_degree(cgraphPairwise, i);
    conflictsPairwise /= 2;
    /*-------------------------------------------------------------------------------------------------------------------*/

    printf("%s %.1lu %.3lf %.1lu %.3lf %.1lu %.3lf\n", problemName, conflictsOsi, osiTime, conflictsProblem,
                                                       problemTime, conflictsPairwise, pairwiseTime);

    assert(cgraph_size(cgraphOsi) == cgraph_size(cgraphProblem) && cgraph_size(cgraphOsi) == cgraph_size(cgraphPairwise));
    for(int i = 0; i < cgraph_size( cgraphProblem ); i++)
        for(int j = i+1; j < cgraph_size( cgraphProblem ); j++)
            assert(cgraph_conflicting_nodes(cgraphProblem, i, j) == cgraph_conflicting_nodes(cgraphOsi, i, j) &&
                   cgraph_conflicting_nodes(cgraphPairwise, i, j) == cgraph_conflicting_nodes(cgraphOsi, i, j));

    cgraph_free( &cgraphProblem );
    cgraph_free( &cgraphOsi );
    cgraph_free( &cgraphPairwise );
    problem_free( &problem );
    problem_free( &problem2 );

















    // clock_t start;
    // double cliqueAnalysis, pairwiseAnalysis;

    // OsiSolverInterface *solver = NULL;
    // OsiClpSolverInterface *realSolver = new OsiClpSolverInterface();
    // solver = (OsiSolverInterface*) realSolver;
    // readLP( argv[1], solver );
    // char problemName[ 256 ];
    // getFileName( problemName, argv[1] );

    // start = clock();
    // CGraph *cgraphClique = osi_build_cgraph( solver );
    // cliqueAnalysis = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    // if ( cgraph_size( cgraphClique ) == 0 )
    // {
    //     printf("EMPTY conflict graph. exiting...\n");
    //     exit(0);
    // }
    // unsigned long int conflictsClique = 0;
    // for(int i = 0; i < cgraph_size( cgraphClique ); i++)
    //     conflictsClique += cgraph_degree(cgraphClique, i);
    // conflictsClique /= 2;
    
    // start = clock();
    // CGraph *cgraphPairwise = osi_build_cgraph_pairwise( solver );
    // pairwiseAnalysis = ((double(clock() - start))/((double)CLOCKS_PER_SEC));
    // if ( cgraph_size( cgraphPairwise ) == 0 )
    // {
    //     printf("EMPTY conflict graph. exiting...\n");
    //     exit(0);
    // }
    // unsigned long int conflictsPairwise = 0;
    // for(int i = 0; i < cgraph_size( cgraphPairwise ); i++)
    //     conflictsPairwise += cgraph_degree(cgraphPairwise, i);
    // conflictsPairwise /= 2;

    // printf("%s %.1lu %.3lf %.1lu %.3lf\n", problemName, conflictsPairwise, pairwiseAnalysis, conflictsClique, cliqueAnalysis);

    // assert(cgraph_size(cgraphClique) == cgraph_size(cgraphPairwise));
    // for(int i = 0; i < cgraph_size( cgraphClique ); i++)
    //     for(int j = i+1; j < cgraph_size( cgraphClique ); j++)
    //         if(cgraph_conflicting_nodes(cgraphPairwise, i, j) != cgraph_conflicting_nodes(cgraphClique, i, j))
    //         {
    //             printf("%d %d\n", cgraph_conflicting_nodes(cgraphPairwise, i, j), cgraph_conflicting_nodes(cgraphClique, i, j));
    //             printf("%d %d\n", i, j);
    //             printf("%s %s\n", i < solver->getNumCols() ? solver->getColName(i).c_str() : ("¬" + solver->getColName(i - solver->getNumCols())).c_str(),
    //                               j < solver->getNumCols() ? solver->getColName(j).c_str() : ("¬" + solver->getColName(j - solver->getNumCols())).c_str() );
    //             i = cgraph_size( cgraphClique );
    //             break;
    //         }


    // cgraph_free( &cgraphClique );
    // cgraph_free( &cgraphPairwise );

    // delete realSolver;

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