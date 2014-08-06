#include <iomanip>
#include <cstring>
#include "BKGraph.hpp"

extern "C"
{
#include "clique_enum.h"
#include "grasp.h"
};

using namespace std;

#define MIN_VIOL      1010
#define TIME_LIMIT      60
#define MAX_GRASP_IT 32768

int main(int argc, char *argv[])
{
    BKGraph graph(argv[1]);
    vector<BKVertex> vertices;

    CGraph *cgraph = cgraph_load(argv[1]);

    BKGraph graph2(cgraph);

    int stat1 = graph.execute(MIN_VIOL, TIME_LIMIT);
    int tWeight1 = graph.writeSolutions("output.txt");
    if(stat1) cout<<"Time Over!"<<endl;
    cout<<"Total weight1: "<<tWeight1<<endl;

    int stat2 = graph2.execute(MIN_VIOL, TIME_LIMIT);
    int tWeight2 = graph2.writeSolutions("output2.txt");
    if(stat2) cout<<"Time Over!"<<endl;
    cout<<"Total weight2: "<<tWeight2<<endl;


    /* basic result check. grasp cannot find more cliques than BK... */
    Grasp *grasp = grasp_create( cgraph, MIN_VIOL );
    grasp_set_max_no_improvement( grasp, MAX_GRASP_IT );
    grasp_set_alpha_choice( grasp, GRASP_ALPHA_REACTIVE );
    grasp_run( grasp );
    grasp_solution_set( grasp );

    char errorMsg[256] = "";

    printf("number of cliques:\n");
    CliqueSet *clqS1 = graph.convertToClqSet();
    CliqueSet *clqS2 = graph2.convertToClqSet();

    printf("\tbk1: %d bk2: %d \n", clq_set_number_of_cliques(clqS1), clq_set_number_of_cliques(clqS2) );
    if ( clq_set_number_of_cliques(clqS1) != clq_set_number_of_cliques(clqS2) )
        sprintf( errorMsg, "BK1 and BK2 computed different number of cliques: %d %d.\n", clq_set_number_of_cliques(clqS1), clq_set_number_of_cliques(clqS2) );

    printf("\tgrasp: %d  \n", clq_set_number_of_cliques(grasp_solution_set(grasp)) );
    if ( clq_set_number_of_cliques(clqS1) < clq_set_number_of_cliques(grasp_solution_set(grasp)) )
        sprintf( errorMsg, "Heuristic found more cliques than BK: %d %d.\n",  clq_set_number_of_cliques(grasp_solution_set(grasp)), clq_set_number_of_cliques(clqS1) );
    clq_set_free( &clqS1 );
    clq_set_free( &clqS2 );
    grasp_free( &grasp );
    cgraph_free( &cgraph );

    if (strlen(errorMsg)>0)
    {
        printf("ERROR\n");
        fprintf( stderr, "ERROR\n");
        printf("%s\n", errorMsg );
        fprintf( stderr, "%s\n", errorMsg );
        exit( EXIT_FAILURE );
    }

    return EXIT_SUCCESS;
}
