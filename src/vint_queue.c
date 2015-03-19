#include <stdlib.h>
#include <stdio.h>
#include "vint_queue.h"
#include "memory.h"

void vint_queue_init( IntQueue *iqueue, int cap )
{
	iqueue->capacity = cap + 1;
	iqueue->front = 0;
	iqueue->back = 0;
	iqueue->elements = xmalloc(sizeof(int) * iqueue->capacity);
}

int vint_queue_is_empty( const IntQueue *iqueue )
{
	return (iqueue->front == iqueue->back);
}

void vint_queue_push( IntQueue *iqueue, int value )
{
	int backPosition = (iqueue->back + 1) % iqueue->capacity;
	if(backPosition == iqueue->front)
	{
		fprintf(stderr, "Queue is full!\n");
      	exit(1);
	}
	iqueue->elements[iqueue->back] = value;
	iqueue->back = (iqueue->back + 1) % iqueue->capacity;
}

void vint_queue_pop( IntQueue *iqueue, int *value )
{
	if(vint_queue_is_empty(iqueue))
	{
		fprintf(stderr, "Queue is empty!\n");
      	exit(1);
	}

	*value = iqueue->elements[iqueue->front];
	iqueue->front = (iqueue->front + 1) % iqueue->capacity;
}

void vint_queue_clean( IntQueue *iqueue )
{
   if ( iqueue->elements )
      free( iqueue->elements );
   iqueue->elements = NULL;
   iqueue->capacity = 0;
   iqueue->front = 0;
   iqueue->back = 0;
}
