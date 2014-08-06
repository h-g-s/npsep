/* 
 * File:   eclq.c
 * Author: haroldo
 *
 * Created on 16 de Março de 2011, 19:56
 */

#include <stdio.h>
#include <stdlib.h>
#include "cgraph.h"
#include "clique_enum.h"
#include "memory.h"


/*
 * 
 */
int main(int argc, char** argv) {
    if (argc<2)
    {
        fprintf( stderr, "Enter problem file name.\n" );
        exit( EXIT_FAILURE );
    }

    char probName[256];
    getFileName( probName, argv[1] );

    CGraph *cgraph = cgraph_load( argv[1] );
    CliqueEnumerator *clqe = clq_enum_create( 1010 );
    
    int *nodesLeft = xmalloc( sizeof(int)*cgraph_size(cgraph));
    int i;
    for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
        nodesLeft[i] = i;
    IntSet is;
    vint_set_init( &is );
    vint_set_add( &is, nodesLeft, cgraph_size(cgraph) );
    
    clq_enum_run( clqe, cgraph, 0, NULL, 0, &is );

    const CliqueSet *clqSet = clq_enum_get_cliques( clqe );

    char clqFName[256];
    sprintf( clqFName, "%s_clq.txt", probName );
    clq_set_save( cgraph, clqSet, clqFName );

    vint_set_clean( &is );
    free( nodesLeft );
    clq_enum_free( &clqe );
    cgraph_free( &cgraph );

    return (EXIT_SUCCESS);
}

