#ifndef OSICGRAPH_H_INCLUDED
#define OSIGRAPH_H_INCLUDED

extern "C"
{
#include "cgraph.h"
}

#include <OsiSolverInterface.hpp>

/* lp should be a pointer to a OsiSolverInterface object */
CGraph *osi_build_cgraph( const OsiSolverInterface *solver );
CGraph *osi_build_cgraph_pairwise( const OsiSolverInterface *solver );

#endif