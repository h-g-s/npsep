#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <math.h>
#include <float.h>
#include <map>
#include <vector>

extern "C"
{
    #include "lp.h"
    #include "strutils.h"
    #include "build_cgraph.h"
}

#define MIN_VIOLATION 0.02
#define MAX_TIME 600
#define MAX_PASSES 999
#define EPS 1e-8

using namespace std;

string optFile, output;
map<string, double> optimals;
char aggrClique = 0, oddHole = 0, knapsack = 0;

void parseParameters( int argc, char **argv)
{
    int i;
    for ( i=1 ; (i<argc) ; ++i )
    {
        char pName[256];
        char pValue[256];
        getParamName( pName, argv[i] );
        getParamValue( pValue, argv[i] );

        if (strcmp( pName, "optFile" )==0)
        {
            optFile = pValue;
            continue;
        }
        if (strcmp( pName, "log" )==0)
        {
            output = pValue;
            continue;
        }
        if (strcmp( pName, "aggrClique" )==0)
        {
            if(strcmp(pValue, "off") == 0)
                aggrClique = 0;
            else if(strcmp(pValue, "on") == 0)
                aggrClique = 1;
            else
                fprintf(stderr, "Invalid value for parameter <aggrClique>: %s\n", pValue);
            continue;
        }
        if (strcmp( pName, "oddHole" )==0)
        {
            if(strcmp(pValue, "off") == 0)
                oddHole = 0;
            else if(strcmp(pValue, "on") == 0)
                oddHole = 1;
            else
                fprintf(stderr, "Invalid value for parameter <oddHole>: %s\n", pValue);
            continue;
        }
        if (strcmp( pName, "knapsack" )==0)
        {
            if(strcmp(pValue, "off") == 0)
                knapsack = 0;
            else if(strcmp(pValue, "on") == 0)
                knapsack = 1;
            else
                fprintf(stderr, "Invalid value for parameter <knapsack>: %s\n", pValue);
            continue;
        }
    }
}

/* Fills a vector with the variable names, including the binary complements */
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
    clock_t start = clock();

    if(argc < 2)
    {
        printf("Invalid number of parameters!\n");
        exit( EXIT_FAILURE );
    }

    FILE *log = NULL;
    if(!output.empty())
    {
        log = fopen(output.c_str(), "a");
        if(!log)
        {
            printf("Could not open the file!\n");
            exit(EXIT_FAILURE);
        }
    }

    char problemName[256], name[256];
    getFileName( problemName, argv[1] );
    parseParameters( argc, argv );

    LinearProgram *mip = lp_create();
    lp_read(mip, argv[1]);
    lp_set_print_messages(mip, 0);
    const int numCols = lp_cols(mip), numRows = lp_rows(mip);
    int pass = 0, newCuts = 0, totalCuts = 0;
    double opt, pTime;

    clock_t graphTime = clock();
    CGraph *cgraph = build_cgraph_lp(mip);
    printf("Conflict graph built in %.2lf seconds.\n", ((double)(clock()-graphTime)) / ((double)CLOCKS_PER_SEC));

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

    printf("numcols %d numrows %d\n", numCols, numRows);

    int status = lp_optimize_as_continuous(mip);

    switch (status)
    {
        case LP_OPTIMAL:
            break;
        case LP_UNBOUNDED:
            fprintf(stderr, "LP_UNBOUNDED");
            exit(EXIT_FAILURE);
        case LP_INFEASIBLE:
            fprintf(stderr, "LP_INFEASIBLE");
            exit(EXIT_FAILURE);
        case LP_FEASIBLE:
            break;
        case LP_INTINFEASIBLE:
            fprintf(stderr, "LP_INTINFEASIBLE");
            exit(EXIT_FAILURE);
        case LP_NO_SOL_FOUND:
            fprintf(stderr, "LP_NO_SOL_FOUND");
            exit(EXIT_FAILURE);
        case LP_ERROR:
            fprintf(stderr, "LP_ERROR");
            exit(EXIT_FAILURE);
        default:
            fprintf(stderr, "lp status not recognized: %d\n", status );
            exit(EXIT_FAILURE);
    }

    double initialBound = lp_obj_value(mip);
    printf("%.2lf %d %d %.7lf", ((double)(clock()-start)) / ((double)CLOCKS_PER_SEC), pass, 0, initialBound);
    if(!optFile.empty())
    {
        printf(" %.7lf %.7lf", opt, abs_mip_gap(initialBound, opt));
    }
    printf("\n");

    do
    {
        pTime = ((double)(clock()-start)) / ((double)CLOCKS_PER_SEC);
        if(pTime > MAX_TIME) break;

        CutPool *cutPool = cut_pool_create(numCols);

        if(aggrClique)
            lp_generate_clique_cuts(mip, cgraph, cutPool, argc, (const char**)argv);

        if(oddHole)
            lp_generate_odd_hole_cuts(mip, cgraph, cutPool);

        if(knapsack)
            lp_generate_knapsack_cuts(mip, cgraph, cutPool);

        cut_pool_update(cutPool);
        newCuts = cut_pool_size(cutPool);

        for(int j = 0; j < newCuts; j++)
        {
            const Cut *cut = cut_pool_get_cut(cutPool, j);
            sprintf(name, "cut_%d", totalCuts++);
            lp_add_row(mip, cut_size(cut), cut_get_idxs(cut), cut_get_coefs(cut), name, 'L', cut_get_rhs(cut));
        }

        cut_pool_free(&cutPool);

        pTime = ((double)(clock()-start)) / ((double)CLOCKS_PER_SEC);
        if(pTime > MAX_TIME) break;

        ++pass;

        if(newCuts)
        {
            status = lp_optimize_as_continuous(mip);

            switch (status)
            {
                case LP_OPTIMAL:
                    break;
                case LP_UNBOUNDED:
                    fprintf(stderr, "LP_UNBOUNDED");
                    exit(EXIT_FAILURE);
                case LP_INFEASIBLE:
                    fprintf(stderr, "LP_INFEASIBLE");
                    exit(EXIT_FAILURE);
                case LP_FEASIBLE:
                    break;
                case LP_INTINFEASIBLE:
                    fprintf(stderr, "LP_INTINFEASIBLE");
                    exit(EXIT_FAILURE);
                case LP_NO_SOL_FOUND:
                    fprintf(stderr, "LP_NO_SOL_FOUND");
                    exit(EXIT_FAILURE);
                case LP_ERROR:
                    fprintf(stderr, "LP_ERROR");
                    exit(EXIT_FAILURE);
                default:
                    fprintf(stderr, "lp status not recognized: %d\n", status );
                    exit(EXIT_FAILURE);
            }

            double sepTime = ((double)(clock()-start)) / ((double)CLOCKS_PER_SEC);
            printf("%.2lf %d %d %.7lf", sepTime, pass, newCuts, lp_obj_value(mip));
            if(!optFile.empty())
                printf(" %.7lf %.7lf", opt, abs_mip_gap(lp_obj_value(mip), opt));
            printf("\n");
        }
    }
    while ( (newCuts>0) && (pass<MAX_PASSES) ) ;

    double totalTime = ((double)(clock()-start)) / ((double)CLOCKS_PER_SEC);
    printf("%s %.2lf %d %d %.7lf", problemName, totalTime, pass - 1, totalCuts, lp_obj_value(mip));

    if(log)
    {
        fprintf(log, "%s %.2lf %d %d %.7lf", problemName, totalTime, pass - 1, totalCuts, lp_obj_value(mip));
        if(!optFile.empty())
            fprintf(log, " %.7lf", abs_mip_gap(lp_obj_value(mip), opt));
        fprintf(log, "\n");
    }

    sprintf(name, "%s_with_cuts.lp", problemName);
    lp_write_lp(mip, name);

   	lp_free(&mip);
    cgraph_free(&cgraph);

    return EXIT_SUCCESS;
}
