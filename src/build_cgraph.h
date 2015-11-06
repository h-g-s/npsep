#ifndef BUILDCGRAPH_H_INCLUDED
#define BUILDCGRAPH_H_INCLUDED

extern "C"
{
	#include "cgraph.h"
}

#include "problem.h"
#include <OsiSolverInterface.hpp>

CGraph *build_cgraph(const OsiSolverInterface *solver);
CGraph *build_cgraph_pairwise(const OsiSolverInterface *solver);

CGraph *build_cgraph(const Problem *problem);
CGraph *build_cgraph_pairwise(const Problem *problem);

#endif