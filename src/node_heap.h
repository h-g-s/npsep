/*
 * NODE_HEAP.h
 *
 * Monotone heap
 * Updates MUST always decrease costs
 *
 * Developed by Haroldo Gambini Santos
 * haroldo.santos@gmail.com
 */

#ifndef NODE_HEAP_H
#define NODE_HEAP_H

typedef struct _NodeHeap NodeHeap;

/* creates the heap with space
 * for nodes {0,...,nodes-1}
 * infinity specifies the constant which will be used to indicate
 * that a node has infinity cost or was removed form the heap
 */
NodeHeap *nh_create( const size_t nodes, const int infinity );

/* updates, always in decreasing order,
 * the cost of a node
 */
void nh_update( NodeHeap * npq, const int node, const int cost );

// removes the next element in priority queue npq,
// returns the cost and fills node
int nh_remove_first( NodeHeap * npq, int *node );

/* returns the current
 * cost of a node
 */
int nh_get_cost( NodeHeap * npq, const int node );

/* sets all costs to
 * infinity again
 */
void nh_reset( NodeHeap * npq );

/* frees
 * the entire memory
 */
void nh_free( NodeHeap **nh );

#endif /* ifndef NODE_HEAP_H */

