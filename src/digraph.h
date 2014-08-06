/*
 * digraph.h
 * Developed by Haroldo Gambini Santos
 * hsantos@ic.uff.br
 */

#ifndef DIGRAPH_H
#define DIGRAPH_H

typedef struct
{
   int tail;
   int head;
   int distance;
} Arc;

typedef struct
{
   int node;
   int distance;
} Neighbor;

#endif /* ifndef DIGRAPH_H */
