#ifndef OSICGRAPH_H_INCLUDED
#define OSIGRAPH_H_INCLUDED

extern "C"
{
#include "cgraph.h"
}

#include <vector>

extern double maxRowTime;
extern int maxRowTimeIdx;
extern int success;

/* lp should be a pointer to a OsiSolverInterface object */
//CGraph *osi_build_cgraph( void *_lp, bool greedyClqExt = true );
CGraph *osi_build_cgraph( void *_lp );
CGraph *osi_build_cgraph_pairwise( void *_lp );

#endif