/* just a C wrapper for BK */

#ifndef BRON_KERBOSCH_H_INCLUDED
#define BRON_KERBOSCH_H_INCLUDED

#include "cgraph.h"
#include "clique.h"

typedef struct _BronKerbosch BronKerbosch;

BronKerbosch *bk_create(const CGraph *cgraph);

int bk_run(BronKerbosch *bk );

const CliqueSet *bk_get_clq_set( BronKerbosch *bk );

int bk_get_max_weight( BronKerbosch *bk );

void bk_free( BronKerbosch *bk );

void bk_set_min_weight(BronKerbosch *bk, int minWeight);
void bk_set_max_it(BronKerbosch *bk, size_t maxIt);


#endif
