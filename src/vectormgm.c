#ifndef VECTORMGM_C_INCLUDED
#define VECTORMGM_C_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include "memory.h"

void vmg_adjust_vector_capacity( void **v, int *cap, const int required, const size_t elementSize )
{
   if ( required <= (*cap) )
      return;

   if ( ((*v)==NULL) || ((*cap)==0) )
   {
      (*v) = malloc( elementSize*required );
      if ( (*v)==NULL )
         goto ERROR;

      (*cap) = required;

      return;
   }
   else
   {
      (*cap) *= 2;
      if ( (*cap) < required )
         (*cap) = required;

      void *other = realloc( (*v), elementSize*(*cap) );
      if ( other==NULL )
         goto ERROR;

      (*v) = other;
   }

   return;

ERROR:
   fprintf( stderr, "ERROR: No memory left. Could not allocate/resize vector\n" );
   exit( EXIT_FAILURE );
}

#endif

