/** queue of integers implemented as a circular vector **/

#ifndef VINT_QUEUE_H_INCLUDED
#define VINT_QUEUE_H_INCLUDED

typedef struct
{
   int *elements;
   int capacity;
   int front;
   int back;
} IntQueue;

void vint_queue_init( IntQueue *iqueue, int cap );

int vint_queue_is_empty( const IntQueue *iqueue ); //returns 1 if queue is empty

void vint_queue_push( IntQueue *iqueue, int value );

void vint_queue_pop( IntQueue *iqueue, int *value );

void vint_queue_clean( IntQueue *iqueue );

#endif

