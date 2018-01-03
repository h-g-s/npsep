/** set of integers implemented as a vector
    fast for queries
    slow for multiple modifications  **/

#ifndef VINT_SET_H_INCLUDED
#define VINT_SET_H_INCLUDED

#include <stddef.h>

typedef struct
{
   int *elements;
   int size;
   int capacity;
} IntSet;


void vint_set_init( IntSet *iset );

const int *vint_set_get_elements( const IntSet *iset );

const int *vint_set_find( const IntSet *iset, int key );

const int vint_set_intersection( int intersection[], const int size, const int elements[], const IntSet *is2 );

int vint_set_size( const IntSet *iset );

void vint_set_add( IntSet *iset, const int elements[], int size );

/* faster function to add elements:
 * WARNING: considers that all elements in elements are sorted and
 * that there are no repeated entries */
void vint_set_add_opt( IntSet *iset, const int elements[], int size );

void vint_set_add_using_original_indexes( IntSet *iset, const int elements[], int size, const int orig[] );

const int *vint_set_int_find( const int key, const int n, int numbers[] );

/**
 * returns 1 if these int_set are equal, 0 otherwise
 **/
int vint_set_equals( const IntSet *is1, const IntSet *is2 );

int vint_set_cmp_int( const void *e1, const void *e2 );

/* be sure that is is ready to receive "required" number of elements.
   if necessary pre-allocate space */
void vint_set_check_capacity( IntSet *iset, const size_t required );

/**
 * frees used memory
 */
void vint_set_clean( IntSet *iset );

/**
 * removes all elements
 */
void vint_set_clear( IntSet *iset );


/** "private" functions
     be CAREFULL to use **/
void vint_set_force_size( IntSet *iset, const int size );

int *vint_set_force_elements_access( IntSet *iset );

/* checks consistency: sorting */
void vint_set_force_check( IntSet *iset );

const int *bsearch_int( const int v[], const size_t n, const int key );

void qsort_int( int v[], const size_t n );

void vint_insert_sort( const int key, int v[], const size_t n );

#endif

