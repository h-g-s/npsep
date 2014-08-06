#ifndef VNODE_DEGREE_H
#define VNODE_DEGREE_H

typedef struct _VNodeDegree VNodeDegree;

VNodeDegree *vndg_create( int iniCap );

void vndg_add( VNodeDegree *vndg, int node, int degree );

int vndg_size( VNodeDegree *vndg );

int vndg_get_node( VNodeDegree *vndg, int pos );

int vndg_get_degree( VNodeDegree *vndg, int pos );

/* sorts from the smallest to the largest */
void vndg_sort( VNodeDegree *vndg );

void vndg_clear( VNodeDegree *vndg );

void vndg_free( VNodeDegree **vndg );

#endif

