#ifndef BUILDCGRAPH_H_INCLUDED
#define BUILDCGRAPH_H_INCLUDED


#include "cgraph.h"
#include "problem.h"

CGraph* build_cgraph_osi(const void *_solver);
CGraph* build_cgraph_pairwise_osi(const void *_solver);

CGraph *build_cgraph(const Problem *problem);
CGraph *build_cgraph_pairwise(const Problem *problem);
CGraph *build_cgraph_pairwise_no_gub(const Problem *problem);

#endif