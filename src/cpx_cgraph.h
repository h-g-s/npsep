
#ifndef CPX_CGRAPH

extern "C"
{
#include "cgraph.h"
}

/* lp should be a pointer to a CPLEX problem object */
CGraph *cpx_build_cgraph( void *_lp );

#endif

