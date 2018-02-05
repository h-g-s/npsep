#ifndef CLIQUE_EXTENDER_H
#define CLIQUE_EXTENDER_H

/**
 * starting from a clique found in a preprocessed graph,
 * extends it to include more nodes from the original graph
 **/

#include "cgraph.h"
#include "vint_set.h"
#include "clique.h"

typedef struct _CliqueExtender CliqueExtender;

typedef enum    /* zero refers to no extension at all */
{
   CLQEM_RANDOM = 1,
   CLQEM_MAX_DEGREE = 2,
   CLQEM_PRIORITY_GREEDY = 3,
   CLQEM_EXACT = 4
} CliqueExtendingMethod;

CliqueExtender *clqe_create();

/**
 * extends clique "clique". extended
 * cliques are stores in internal cliqueSet
 **/

int clqe_extend( CliqueExtender *clqe, const CGraph *cgraph, const IntSet *clique,
                 const int weight, const CliqueExtendingMethod clqem );

const CliqueSet *clqe_get_cliques( CliqueExtender *clqe );

/* sets up costs for n variables */
void clqe_set_costs( CliqueExtender *clqe, const int costs[], const int n );

const int *clqe_get_costs( CliqueExtender *clqe );

void clqe_set_clear( CliqueExtender *clqe );

void clqe_free( CliqueExtender **clqe );


/* parameters */
void clqe_set_max_candidates( CliqueExtender *clqe, const int max_size );
int clqe_get_max_candidates( CliqueExtender *clqe );
void clqe_set_max_cost( CliqueExtender *clqe, const int maxCost );
int clqe_get_max_cost( CliqueExtender *clqe );
void clqe_set_max_clq_gen( CliqueExtender *clqe, const int maxClqGen );
int clqe_get_max_clq_gen( CliqueExtender *clqe );
void clqe_set_rc_percentage( CliqueExtender *clqe, const double rcPercentage );
double clqe_get_rc_percentage( CliqueExtender *clqe );

#endif
