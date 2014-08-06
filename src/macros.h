#ifndef CGRAPH_MACROS
#define CGRAPH_MACROS

#include <stdlib.h>
#include <stdio.h>

#define INT_RANDOM( n ) \
   ((int)( ((double)n) * (((double)rand())/(((double)RAND_MAX)+((double)1.0))) ))

#define MAX( a, b ) ( (a)>(b) ? (a) : (b) )
#define MIN( a, b ) ( (a)<(b) ? (a) : (b) )

#endif

