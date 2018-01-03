#ifndef VECTORMGM_H_INCLUDED
#define VECTORMGM_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>

/**
 * allocates a vector of strings
 **/
#define CREATE_STRING_VECTOR( svector, lines, columns ) \
   svector = malloc( sizeof(char*)*lines ); \
   if (!svector) { \
      fprintf( stderr, "no memory for string vector.\n" ); abort(); } \
   svector[0] = malloc( sizeof(char)*lines*columns ); \
   if (!svector[0]) { \
      fprintf( stderr, "no memory for string vector.\n" ); abort(); } \
   { \
      int i; \
      for ( i=1 ; (i<lines) ; ++i ) \
         svector[i] = svector[i-1] + columns; \
   };

#define FREE_STRING_VECTOR( svector ) free( svector[0] ); free( svector ); svector = NULL;

/**
 * set of macros for the
 * management of vectors
 **/

#define ALLOCATE_INT_VECTOR( vector, size ) \
   vector = malloc( sizeof(int)*size ); \
   if (!vector)  { \
      fprintf( stderr, "Error: at: %s:%d . No memory to allocate a vector with %d integers.\n", __FILE__, __LINE__, size ); \
      exit( EXIT_FAILURE ); \
   }

#define RESIZE_INT_VECTOR( vector, newSize ) \
{ \
   int *temp = realloc( vector, (sizeof(int)*newSize) ); \
   if (!temp)  { \
      fprintf( stderr, "Error: at: %s:%d . No memory to allocate a vector with %d integers.\n", __FILE__, __LINE__, newSize ); \
      exit( EXIT_FAILURE ); \
   } \
   vector = temp; \
}

/**
 * adjusts capacity if necessary,
 * resizing vector and preserving its contents
 * expands to at least twice the size at time
 **/
#define ADJUST_INT_VECTOR_CAPACITY( vector, capacity, required ) \
{\
   if ( (capacity) < (required) ) \
   { \
      (capacity) *= 2; \
      if ( (capacity) < (required) ) \
         (capacity) = (required); \
      int *temp = realloc( vector, (capacity)*sizeof(int) ); \
      if (!temp) \
      { \
         fprintf( stderr, "Error trying to allocate a vector with %d integers. No more memory available. Error at %s:%d.\n", (required), __FILE__, __LINE__ ); \
         exit( EXIT_FAILURE ); \
      } \
      vector = temp; \
   } \
}

#define ADJUST_VECTOR_CAPACITY( vector, capacity, required ) \
{\
   if ( (capacity) < (required) ) \
   { \
      (capacity) *= 2; \
      if ( (capacity) < (required) ) \
         (capacity) = (required); \
      typeof(vector) temp = realloc( vector, (capacity)*sizeof(typeof( (*vector) )) ); \
      if (!temp) \
      { \
         fprintf( stderr, "Error trying to allocate a vector. No more memory available. Error at %s:%d.\n", __FILE__, __LINE__ ); \
         exit( EXIT_FAILURE ); \
      } \
      vector = temp; \
   } \
}

#define VECTOR_ADD( vector, element, size, capacity ) \
   if ( ( capacity==0 ) || ( vector==NULL ) ) \
   { \
      capacity = 512; \
      vector = malloc( sizeof(typeof(*vector))*capacity ); \
      if (!vector) \
      { \
         fprintf( stderr, "ERROR: no memory. Source: %s Line: %d\n", __FILE__, __LINE__ ); \
         exit( EXIT_FAILURE ); \
      } \
   } else \
   { \
      ADJUST_VECTOR_CAPACITY( vector, capacity, size+1 ); \
      vector[size] = element; \
      ++size; \
   }

void vmg_adjust_vector_capacity( void **v, int *cap, const int required, const size_t elementSize );

#endif

