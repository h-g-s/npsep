#include <stdlib.h>
#include <assert.h>
#include "vnode_degree.h"
#include "memory.h"

typedef struct
{
   int node;
   int degree;
} INodeDegree;

struct _VNodeDegree
{
   int capacity;
   int size;

   INodeDegree *items;
};

int cmp_inode_degree( const void *v1, const void *v2 )
{
   const INodeDegree *ind1 = v1;
   const INodeDegree *ind2 = v2;

   return ( ind1->degree - ind2->degree );
}

VNodeDegree *vndg_create( int iniCap )
{
   VNodeDegree *result = xmalloc( sizeof(VNodeDegree) );
   
   result->capacity = iniCap;
   result->size = 0;
   result->items = xmalloc( sizeof(INodeDegree)*iniCap );
   
   return result;
}

int vndg_size( VNodeDegree *vndg )
{
   return vndg->size;
}

void vndg_add( VNodeDegree *vndg, int node, int degree )
{
   if ( (vndg->size+1) >= vndg->capacity )
   {
      vndg->capacity *= 2;
      vndg->items = xrealloc( vndg->items,  sizeof(INodeDegree)*vndg->capacity );
   }

   vndg->items[vndg->size].node = node;
   vndg->items[vndg->size].degree = degree;

   vndg->size++;
}

int vndg_get_node( VNodeDegree *vndg, int pos )
{
#ifdef DEBUG
   assert( pos >= 0);
   assert( pos < vndg->size );
#endif
   return vndg->items[pos].node;
}

int vndg_get_degree( VNodeDegree *vndg, int pos )
{
#ifdef DEBUG
   assert( pos >= 0);
   assert( pos < vndg->size );
#endif

   return vndg->items[pos].degree;
}

void vndg_sort( VNodeDegree *vndg )
{
   qsort( vndg->items, vndg->size, sizeof(INodeDegree), cmp_inode_degree );
}

void vndg_clear( VNodeDegree *vndg )
{
   vndg->size = 0;
}

void vndg_free( VNodeDegree **vndg )
{
   free( (*vndg)->items );
   free( (*vndg) );
   (*vndg) = NULL;
}

