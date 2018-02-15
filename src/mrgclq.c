#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <float.h>
#include "clique_merge.h"
#include "lp.h"
#include "build_cgraph.h"

extern int clqMergeVerbose;
extern double clqMergeSecsCheckClique;
extern double clqMergeSecsExtendAndDominate;
extern double clqMergeSecsExtend;
extern double clqMergeSecsAddAndRemove;
extern int clqMergeNExtended;
extern int clqMergeNDominatedFull;
extern int clqMergeNDominatedEqPart;
double nlb = DBL_MAX;



int main( int argc, const char **argv )
{
    if (argc<6)
    {
        fprintf( stderr, "usage: lpFileName maxExtensions threads verbose maxItBk\n");
        exit( EXIT_FAILURE );
    }

    LinearProgram *mip = lp_create();
    lp_read( mip, argv[1] );

    printf("lp read with dimensions (cols, rows, nzs): (%d, %d, %d)\n", lp_cols(mip), lp_rows(mip), lp_nz(mip) );
    printf("creating cgraph ... ");

    clock_t startcg = clock();
    CGraph *cgraph = build_cgraph_lp( mip );
    printf("done in %.4f seconds\n", ((double)clock()-startcg) / ((double)CLOCKS_PER_SEC) );

    int maxExt = atoi(argv[2]);

    printf("each clique will be extended to at most %d new cliques\n", maxExt );
    
    int threads = atoi( argv[3] );

    omp_set_num_threads( threads );
    
    printf("processing will use %d threads.\n", threads );

    clqMergeVerbose = atoi(argv[4]);

    int maxItBk = atoi(argv[5]);

    double lb = DBL_MAX;
    int status = lp_optimize_as_continuous( mip );
    if (status == LP_OPTIMAL)
        lb = lp_obj_value( mip );

    if (cgraph)
    {
        merge_cliques( mip, cgraph, maxExt, maxItBk );
        
        lp_write_lp( mip, "pp.lp" );
        
        status = lp_optimize_as_continuous( mip );
        if (status == LP_OPTIMAL)
            nlb = lp_obj_value( mip );
         
    }
    char first = 1;
    {
        FILE *ff = fopen( "summary.csv", "r" );
        if (ff)
        {
            fclose(ff);            
            first = 0;        
        }
    }        
        
    FILE *f = fopen( "summary.csv", "a" );
    if (first)
        fprintf( f, "time discover cliques,time ext merge and dominate,time ext,time add and remove,nr. ext. constraints,nr. dom. constraints leq,nr. dom. constraints eq,lb,nlb\n" );
            
    fprintf( f, "%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%g,%g\n", clqMergeSecsCheckClique, clqMergeSecsExtendAndDominate, clqMergeSecsExtend, clqMergeSecsAddAndRemove, clqMergeNExtended, clqMergeNDominatedFull, clqMergeNDominatedEqPart, lb, nlb );
    fclose( f );

    cgraph_free( &cgraph );
    lp_free( &mip );

    return EXIT_SUCCESS;
}
