#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stddef.h>
#include "memory.h"
#include "macros.h"
#include "vint_set.h"

#define MIN_CAP 128
#define LARGE_VECTOR 4096

void vint_set_init( IntSet *iset )
{
   iset->elements = NULL;
   iset->size = 0;
   iset->capacity = 0;
}

int vint_set_cmp_int( const void *e1, const void *e2 )
{
   int i1 = (*((const int*) e1));
   int i2 = (*((const int*) e2));

   return (i1-i2);
}

void vint_set_add_using_original_indexes( IntSet *iset, const int elements[], int size, const int orig[] )
{
   const int newSize = iset->size + size;

   vint_set_check_capacity( iset, newSize );

   memcpy( iset->elements + iset->size, elements, sizeof(int)*size );
   /* going to original indexes */
   {
      int *ptr = iset->elements + iset->size;
      int *end = ptr + size;
      for ( ; (ptr<end) ; ++ptr )
         *ptr = orig[*ptr];
   }

   qsort_int( iset->elements, newSize );

   iset->size = newSize;

   /** removing duplicate entries **/
   int lastOne = INT_MAX;
   int removals = 0;
   {
      int i;
      for ( i=0 ; ( i<newSize ) ; ++i )
      {
         if ( iset->elements[i] == lastOne )
         {
            iset->elements[i] = INT_MAX;
            removals++;
         }
         else
         {
            lastOne = iset->elements[i];
         }
      }
   }

   if ( removals )
   {
      qsort_int( iset->elements, newSize );
      iset->size -= removals;
   }
}

void vint_set_add( IntSet *iset, const int elements[], int size )
{
   const int newSize = iset->size + size;

   vint_set_check_capacity( iset, newSize );

   /* if there are few new elements, adding as insetion sort */
   if ( size<=3 )
   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
      {
         const int key = elements[i];
         if ( !bsearch_int( iset->elements, iset->size, key  ) )
         {
            vint_insert_sort( key, iset->elements, iset->size );
            iset->size++;
         }
      }
      return;
   }

   memcpy( iset->elements + iset->size, elements, sizeof(int)*size );
   qsort_int( iset->elements, newSize );

   iset->size = newSize;

   /** removing duplicate entries **/
   int lastOne = INT_MAX;
   int removals = 0;
   {
      int i;
      for ( i=0 ; ( i<newSize ) ; ++i )
      {
         if ( iset->elements[i] == lastOne )
         {
            iset->elements[i] = INT_MAX;
            removals++;
         }
         else
         {
            lastOne = iset->elements[i];
         }
      }
   }

   if ( removals )
   {
      qsort_int( iset->elements, newSize );
      iset->size -= removals;
   }
}

void vint_set_add_opt( IntSet *iset, const int elements[], int size )
{
#ifdef DEBUG
   {
      int i;
      for ( i=1 ; (i<size) ; ++i )
      {
         if ( elements[i-1] > elements[i] )
         {
            int j;
            printf("\n");
            for ( j=0 ; (j<size) ; ++j )
               printf("%d ", elements[j]);
            printf("\n");
            fprintf( stderr, "ERROR: vint_set_add_opt - elements are not sorted in the expected order.\n" );
            fflush(stdout); fflush(stderr);
            abort();
            exit(EXIT_FAILURE);
         }
         if ( elements[i-1] == elements[i] )
         {
            fprintf( stderr, "ERROR: vint_set_add_opt - repeated elements are included.\n" );
            exit(EXIT_FAILURE);
         }
      }
   }
#endif
   int wasEmpty = iset->size == 0;
   const int newSize = iset->size + size;

   vint_set_check_capacity( iset, newSize );

   if (wasEmpty)
   {
      memcpy( iset->elements, elements, sizeof(int)*size );
      iset->size = size;
      return;
   }

   /* if there are few new elements, adding as insetion sort */
   if ( size<=3 )
   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
      {
         const int key = elements[i];
         if ( !bsearch_int( iset->elements, iset->size, key  ) )
         {
            vint_insert_sort( key, iset->elements, iset->size );
            iset->size++;
         }
      }
      return;
   }

   memcpy( iset->elements + iset->size, elements, sizeof(int)*size );

   qsort_int( iset->elements, newSize );

   iset->size = newSize;

   /** removing duplicate entries **/
   int lastOne = INT_MAX;
   int removals = 0;
   {
      int i;
      for ( i=0 ; ( i<newSize ) ; ++i )
      {
         if ( iset->elements[i] == lastOne )
         {
            iset->elements[i] = INT_MAX;
            removals++;
         }
         else
         {
            lastOne = iset->elements[i];
         }
      }
   }

   if ( removals )
   {
      qsort_int( iset->elements, newSize );
      iset->size -= removals;
   }
}

const int *vint_set_get_elements( const IntSet *iset )
{
   return iset->elements;
}

int vint_set_size( const IntSet *iset )
{
   return iset->size;
}

const int *vint_set_find( const IntSet *iset, int key )
{
   if (!iset->size)
      return NULL;

   return bsearch_int( iset->elements, iset->size, key );
}

const int vint_set_intersection( int intersection[], const int size, const int elements[], const IntSet *is2 )
{
   int result = 0;

   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
      {
         if ( vint_set_find( is2, elements[i] ) )
            intersection[result++] = elements[i];
      }
   }

   return result;
}

int *vint_set_int_find( const int key, const int n, int numbers[] )
{
   return bsearch_int( numbers, n, key );
}

int vint_set_equals( const IntSet *is1, const IntSet *is2 )
{
   if ( vint_set_size( is1 ) != vint_set_size( is2 ) )
      return 0;

   const size_t size = vint_set_size( is1 );
   const int *el1 = vint_set_get_elements( is1 );
   const int *el2 = vint_set_get_elements( is2 );

   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
         if ( el1[i] != el2[i] )
            return 0;
   }

   return 1;
}

void vint_set_clear( IntSet *iset )
{
   iset->size = 0;
}

inline void vint_set_check_capacity( IntSet *iset, const size_t required )
{
#ifdef DEBUG
   assert( required > 0 );
   assert( iset->capacity >=0 );
   if ( iset->size > iset->capacity )
   {
      fprintf( stderr, "Error: vector with size %d and capacity %d.\n", iset->size, iset->capacity );
      fflush( stderr );
      abort();
   }
#endif
   if ( required <= iset->capacity )
      return;

   /* need to increase capacity */
   if ( iset->capacity == 0 )
   {
      iset->capacity = MAX( required, MIN_CAP );

      iset->elements = xmalloc( sizeof(int)*(iset->capacity) );
   }
   else
   {
      if ( required > LARGE_VECTOR )
         iset->capacity = (int)(1.5*((double)iset->capacity));
      else
         iset->capacity *= 2;

      if ( required > iset->capacity )
         iset->capacity = required;

      iset->elements = xrealloc( iset->elements, sizeof(int)*iset->capacity );
   } /* xrealloc */
}

void vint_set_force_size( IntSet *iset, const int size )
{
   iset->size = size;
}

void vint_set_clean( IntSet *iset )
{
   if ( iset->elements )
      free( iset->elements );
   iset->elements = NULL;
   iset->size = 0;
   iset->capacity = 0;
}

void vint_set_force_check( IntSet *iset )
{
   int i;
   const int lastN = vint_set_size(iset)-1;
   for ( i=0 ; (i<lastN) ; i++ )
   {
      if (iset->elements[i]>iset->elements[i+1])
      {
         fprintf( stderr, "ERROR: IntSet should have sorted elements.\n" );
         exit( EXIT_FAILURE );
      }
   }
}

int *vint_set_force_elements_access( IntSet *iset )
{
   return iset->elements;
}

inline void qsort_int( int v[], const size_t n )
{
#define LOW  0
#define HIGH 1
#define SWAP(a,b) { int t=a; a=b; b=t; }
#define MAX_STACK_HEIGHT 64
#define PIVOT_MEDIAN_MIN_SIZE 16 /* minimum size so that pivot selection is made using median */
  int stack[MAX_STACK_HEIGHT][2];
  int l,r,i,j,p;

  /* initialize our stack */
  stack[0][LOW] = 0;
  stack[0][HIGH] = n-1;
  int depth = 1;

  int selements[3];

PARTITION:
  --depth;

  l = stack[depth][LOW];
  r = stack[depth][HIGH];

  int psize = r - l + 1;

  if (psize < 2)
    goto NEXT_PARTITION;

  /* pivot selection: picking the median of three elements
     if the partition is larger than five */
  if ( psize >= PIVOT_MEDIAN_MIN_SIZE )
  {
    selements[0] = v[l+rand()%psize];
    selements[1] = v[l+rand()%psize];
    selements[2] = v[l+rand()%psize];
    if ( selements[0] > selements[1] )
      SWAP(selements[0],selements[1]);
    if ( selements[1] > selements[2] )
      SWAP(selements[1],selements[2]);
    if ( selements[0] > selements[2] )
      SWAP(selements[0],selements[2]);
    p = selements[1];
  }
  else
    p = v[l+rand()%psize];

  i = l;
  j = r;

  do
  {
    while (v[i]<p) ++i;
    while (v[j]>p) --j;

    if (i<=j)
    {
      SWAP( v[i], v[j] );
      ++i;
      --j;
    }
  } while( i<= j);

  if (l<j)
  {
    stack[depth][LOW] = l;
    stack[depth][HIGH] = j;
    depth++;
  }
  if (i<r)
  {
    stack[depth][LOW] = i;
    stack[depth][HIGH] = r;
    depth++;
  }

NEXT_PARTITION:
  if (depth > 0)
    goto PARTITION;

#undef MAX_STACK_HEIGHT
#undef PIVOT_MEDIAN_MIN_SIZE
#undef SWAP
#undef LOW
#undef HIGH
}

inline int *bsearch_int( int v[], const size_t n, const int key )
{
  register int l = 0;
  register int r = n-1;
  register int m;
  while ( l<= r )
  {
    m = (l+r) / 2;
    if ( v[m] == key )
      return v+m;
    else
    {
      if ( key<v[m] )
        r = m-1;
      else
        l = m+1;
    }
  }

  return NULL;
}

/* inserts given element in an appropriated position in a
 * vector keepking everything sorted
 * one should reserve the appropriated memory space before calling this */
inline void vint_insert_sort( const int key, int v[], const size_t n )
{
#ifdef DEBUG
   {
      int i;
      for ( i=1 ; (i<n) ; ++i )
      {
         if ( v[i-1] > v[i] )
         {
            fprintf( stderr, "ERROR: passing an unsorted vector to function vint_insert_sort\n" );
            exit( EXIT_FAILURE );
         }
      }
   }
#endif
   /* doing a binary search */
   register int l = 0;
   register int r = n-1;
   register int m;
   register int ip = -1;  /* insertion pos */
   while ( l<= r )
   {
      m = (l+r) / 2;
      if ( v[m] == key )
      {
         ip = m;
         while ( (ip>=1) && (v[ip-1]==v[ip]) )
            --ip;
         goto FOUND_POS;
      }
      else
      {
         if ( key<v[m] )
            r = m-1;
         else
            l = m+1;
      }
   }
FOUND_POS:
   if ( ip == -1 )
      ip = l;

   size_t msize = sizeof(int) * (n - ip);
   memmove( (v+ip+1), (v+ip), msize );
   v[ip] = key;
}

