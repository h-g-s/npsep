#include "bron_kerbosch.h"
#include "BKGraph.hpp"
#include "clique.h"

struct _BronKerbosch
{
    BKGraph *bkg;
    CliqueSet *clqSet;
};

extern "C" {
   BronKerbosch *bk_create( const CGraph *cgraph );

BronKerbosch *bk_create( const CGraph *cgraph )
{
    BronKerbosch *result = (BronKerbosch *) xmalloc( sizeof(BronKerbosch) );
    result->bkg = new BKGraph( cgraph );
    result->clqSet = NULL;

    return result;
}

int bk_run( BronKerbosch *bk, const int minViol, const int timeLimit )
{
    int status = bk->bkg->execute( minViol, timeLimit );

    if (bk->clqSet)
    {
        clq_set_free( &(bk->clqSet) );
        bk->clqSet = NULL;
    }

    bk->clqSet = bk->bkg->convertToClqSet();

    return status;
}

CliqueSet *bk_get_clq_set( BronKerbosch *bk )
{
    return bk->clqSet;
}

void bk_free( BronKerbosch *bk )
{
    if ((bk)->clqSet)
    {
        clq_set_free( &((bk)->clqSet) );
        (bk)->clqSet = NULL;
    }

    if (bk->bkg)
    {
       delete bk->bkg;
       bk->bkg = NULL;
    }

    free( bk );
}

}


