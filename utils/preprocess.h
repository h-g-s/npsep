#ifndef PREPROCESS_H
#define PREPROCESS_H

extern "C"
{
    #include "cgraph.h"
    #include "problem.h"
}

typedef struct _Preprocess Preprocess;

Preprocess* preprocess_create(const Problem *problem);
void preprocess_free(Preprocess **pp);
Problem* preprocess_basic_preprocessing(Preprocess *pp);

void cpropagation_get_vars_to_fix(Preprocess *pp, const CGraph* cgraph);

#endif