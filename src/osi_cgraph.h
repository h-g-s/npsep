
#ifndef OSI_CGRAPH

extern "C"
{
#include "cgraph.h"
}

extern double maxRowTime;
extern int maxRowTimeIdx;
extern int success;

/* lp should be a pointer to a OsiSolverInterface object */
CGraph *osi_build_cgraph( void *_lp );

#endif
