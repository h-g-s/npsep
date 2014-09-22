
#ifndef OSI_CGRAPH

extern "C"
{
#include "cgraph.h"
}

extern double maxRowTime;
extern int maxRowTimeIdx;
extern int success;
extern unsigned long int completeClique, incompleteClique, completeCliqueComplement, incompleteCliqueComplement,
					activePairwise, inactivePairwise, mixedPairwise;

/* lp should be a pointer to a OsiSolverInterface object */
CGraph *osi_build_cgraph( void *_lp );

#endif
