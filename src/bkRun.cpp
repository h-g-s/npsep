#include <cstdio>
#include <cstdlib>
extern "C"
{
#include "cgraph.h"
};
#include "bron_kerbosch.h"
#include "clique.h"

int main( int argc, char ** argv )
{
   CGraph *cgraph = cgraph_load( argv[1] );

   BronKerbosch *bk = bk_create( cgraph );

   bk_run( bk, 1010, 300 );

   const CliqueSet *clqS = bk_get_clq_set( bk );

   clq_set_save( cgraph, clqS, "c1.clq" );

   bk_free( bk );
}

