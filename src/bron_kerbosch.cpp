#include "BKGraph.hpp"

extern "C" {
    #include "bron_kerbosch.h"
}

struct _BronKerbosch
{
    BKGraph *bkg;
    CliqueSet *clqSet;
    int minWeight;
    size_t maxIt;
};

BronKerbosch *bk_create(const CGraph *cgraph)
{
    BronKerbosch *result = (BronKerbosch *) xmalloc( sizeof(BronKerbosch) );
    result->bkg = new BKGraph(cgraph);
    result->clqSet = NULL;
    return result;
}

int bk_run(BronKerbosch *bk)
{
    if (bk->clqSet)
    {
        clq_set_free( &(bk->clqSet) );
        bk->clqSet = NULL;
    }

    int status = bk->bkg->execute();
    bk->clqSet =bk->bkg->getCliqueSet();

    return status;
}

void bk_set_min_weight(BronKerbosch *bk, int minWeight)
{
    bk->minWeight = minWeight;
    bk->bkg->setMinWeight(minWeight);
}

void bk_set_max_it(BronKerbosch *bk, size_t maxIt)
{
    bk->maxIt = maxIt;
    bk->bkg->setMaxIt(maxIt);
}

const CliqueSet *bk_get_clq_set( BronKerbosch *bk )
{
    return bk->clqSet;
}

int bk_get_max_weight( BronKerbosch *bk )
{
    return bk->bkg->getMaxWeight();
}

void bk_free( BronKerbosch *bk )
{
    delete bk->bkg;
	bk->bkg = NULL;
    free( bk );
}
