/* just a C wrapper for BK */

#ifndef BRON_KERBOSCH_H_INCLUDED
#define BRON_KERBOSCH_H_INCLUDED

typedef struct _BronKerbosch BronKerbosch;

#ifdef __cplusplus
extern "C"
{
#include "cgraph.h"
#include "clique.h"

BronKerbosch *bk_create(const CGraph *cgraph );

int bk_run(BronKerbosch *bk, const int minViol, const int timeLimit );

CliqueSet *bk_get_clq_set( BronKerbosch *bk );

void bk_free( BronKerbosch *bk );

}
#else

BronKerbosch *bk_create(const CGraph *craph);

int bk_run(BronKerbosch *bk, const int minViol, const int timeLimit );

CliqueSet *bk_get_clq_set( BronKerbosch *bk );

void bk_free( BronKerbosch *bk );

#endif




/* returns 0 if complete search was done, 1 otherwise */


#endif
