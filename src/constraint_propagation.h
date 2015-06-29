#ifndef CONSTRAINT_PROPAGATION_H
#define CONSTRAINT_PROPAGATION_H

extern "C"
{
	#include "cgraph.h"
}

#include <OsiSolverInterface.hpp>

typedef struct _CPropagation CPropagation;

CPropagation *cpropagation_create(const OsiSolverInterface *solver);
void cpropagation_free(CPropagation *cp);
void cpropagation_get_vars_to_fix(CPropagation *cp);
int cpropagation_get_num_vars_to_fix(CPropagation *cp);
OsiSolverInterface* cpropagation_preprocess(CPropagation *cp, int nindexes[]);

#endif