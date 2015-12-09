#include <cstdio>
#include <cstdlib>
extern "C"
{
#include "cgraph.h"
#include "strUtils.h"
};
#include "bron_kerbosch.h"
#include "clique.h"

int main( int argc, char ** argv )
{
   CGraph *cgraph = cgraph_load( argv[1] );
   BronKerbosch *bk = bk_create( cgraph );

   char problemName[ 256 ];
   getFileName( problemName, argv[1] );

   bk_run( bk, 1010, 300 );

   const CliqueSet *clqS = bk_get_clq_set( bk );

   clq_set_save( cgraph, clqS, "c1.clq" );

   bk_free( bk );
}
