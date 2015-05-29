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

void test(CPropagation *cp);

#endif