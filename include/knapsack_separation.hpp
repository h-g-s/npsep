#ifndef KNAPSACK_SEPARATION_H_INCLUDED
#define KNAPSACK_SEPARATION_H_INCLUDED

extern "C"
{
    #include "cut.h"
    #include "cgraph.h"
    #include "lp.h"
}

#include <OsiSolverInterface.hpp>

typedef struct _KnapsackSeparation KnapsackSeparation;

KnapsackSeparation* knapsack_sep_create(const OsiSolverInterface &si);
KnapsackSeparation* knapsack_sep_create(const LinearProgram *lp);
void knapsack_sep_free(KnapsackSeparation **_ksep);
int knapsack_sep_separate(KnapsackSeparation *ksep, const double *x, const double *rc, CutPool *cutPool);
void knapsack_sep_set_graph(KnapsackSeparation *ksep, const CGraph *cg);
void knapsack_sep_set_params_parse_cmd_line( KnapsackSeparation *ksep, const int argc, const char **argv );

int knapsack_sep_get_num_cols(const KnapsackSeparation *ksep);
int knapsack_sep_get_num_rows(const KnapsackSeparation *ksep);

#endif // KNAPSACK_SEPARATION_H_INCLUDED
