#ifndef PREPROCESS_H
#define PREPROCESS_H

#include "problem.h"

typedef struct _Preprocess Preprocess;

Preprocess* preprocess_create(const Problem *problem);
void preprocess_free(Preprocess **pp);
Problem* preprocess_basic_preprocessing(Preprocess *pp);

#endif