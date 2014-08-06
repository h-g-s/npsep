#ifndef CLIQUE_ENUM_H_INCLUDED
#define CLIQUE_ENUM_H_INCLUDED

#include "clique.h"

typedef struct _CliqueEnumerator CliqueEnumerator;

CliqueEnumerator *clq_enum_create( const int minWeight );

/* receives a clique (possible empty) and tries to extend using nodesLeft */
void clq_enum_run( CliqueEnumerator *clq_enum, const CGraph *cgraph,
                  const int size, const int clique[], const int weight,
                  const IntSet *nodesLeft );

void clq_enum_set_min_weight( CliqueEnumerator *clqEnum, const int minW );

const CliqueSet *clq_enum_get_cliques( const CliqueEnumerator *clq_enum );

void clq_enum_free( CliqueEnumerator **clq_enum );

#endif 

