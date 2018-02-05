#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "clique_merge.h"
#include "lp.h"
#include "build_cgraph.h"

extern int clqMergeVerbose;

int main( int argc, const char **argv )
{
    if (argc<5)
    {
        fprintf( stderr, "usage: lpFileName maxExtensions threads verbose\n");
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

    merge_cliques( mip, cgraph, maxExt );
    
    lp_write_lp( mip, "pp.lp" );

    cgraph_free( &cgraph );
    lp_free( &mip );

    return EXIT_SUCCESS;
}
