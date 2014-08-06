/*
 * spaths.c
 * Developed by Haroldo Gambini Santos
 * haroldo@iceb.ufop.br
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <ctype.h>
#include "digraph.h"
#include "memory.h"
#include "spaths.h"
#include "node_heap.h"
#include "strUtils.h"
#include "vectormgm.h"

struct _ShortestPathsFinder
{
   int capnodes;
   int caparcs;

   int nodes;
   int arcs;

   // all neighbors
   Neighbor *neighs;
   // start of neighbors for node i
   // the neighbor ends at startn[i+1]
   Neighbor **startn;

   // solution
   int *dist;
   int *previous;
   int *path;     // temporary storage for path

   NodeHeap *nh;

   // Floyd Warshall data
   // for computing data
   int fwCapNodes;
   int fwCapArcs;

   int **fwDist;
   int **fwPrev;
};

// lexicographical comparison
int compArcsLX( void *a1, void *a2 )
{
   Arc *pa1 = a1;
   Arc *pa2 = a2;

   if (pa1->tail !=  pa2->tail)
      return (pa1->tail - pa2->tail);

   return ( pa1->head - pa2->head );
}

int compNeighs( const void *n1, const void *n2 )
{
   const Neighbor *pn1 = n1;
   const Neighbor *pn2 = n2;

   return pn1->node - pn2->node;
}

/*
 * returns a pointer to the start of neighbors list of a node
 */
Neighbor *spf_start_n( ShortestPathsFinder* spf, const int node );

/*
 * returns a pointer to the end of neighbors list of a node
 */
Neighbor *spf_end_n( ShortestPathsFinder* spf, const int node );


/*
 * updates the working graph, the new graph
 * can have a different number of nodes/arcs
 */
void spf_update_digraph( ShortestPathsFinder* spf, const int nodes, const int narcs, Arc *arcs );

// Floyd Warshall computing
// space
/*
 * allocates (if needed) and initializes all fw data
 */
void allocateFWSpace( ShortestPathsFinder* spf );
/*
 * frees (if needed)all fw data
 */
void freeFWSpace( ShortestPathsFinder* spf );

void spf_proccessFWLoop( ShortestPathsFinder* spf ) __attribute__((hot));

ShortestPathsFinderPtr spf_create( )
{
   ShortestPathsFinderPtr result;

   result = (ShortestPathsFinderPtr) xmalloc( sizeof( ShortestPathsFinder ) );

   result->capnodes   = 0;
   result->caparcs    = 0;
   result->nodes      = 0;
   result->arcs       = 0;
   result->neighs     = NULL;
   result->startn     = NULL;
   result->previous   = NULL;
   result->nh         = NULL;
   result->dist       = NULL;
   result->path       = NULL;

   // fw data
   result->fwCapNodes = 0;
   result->fwCapArcs  = 0;
   result->fwDist     = NULL;
   result->fwPrev     = NULL;

   return result;
}

void spf_find( ShortestPathsFinder* spf, const int origin )
{
   NodeHeap *nh = spf->nh;

   nh_reset( nh );
   {
      int i;
      for ( i=0 ; (i<spf->nodes) ; i++ )
         spf->dist[i] = SP_INFTY_DIST;
      
      for ( i=0 ; (i<spf->nodes) ; i++ )
         spf->previous[i] = NULL_NODE;
   }

   spf->dist[origin] = 0;

   nh_update( nh, origin, 0 );

   int topCost, topNode;
   while ( (topCost=nh_remove_first( nh, &topNode )) < SP_INFTY_DIST )
   {
      //printf("top node: %d\n", topNode );
      // updating neighbors distances
      // by iterating in all neighbors
      Neighbor *n    = spf->startn[topNode];
      Neighbor *endN = spf->startn[topNode+1];
      for ( ; (n<endN) ; n++ )
      {
         const int toNode  = n->node;
         const int dist    = n->distance;
         const int newDist = topCost + dist;
         if ( spf->dist[ toNode ] > newDist )
         {
            spf->previous[ toNode ] = topNode;
            spf->dist[ toNode ]     = newDist;
            nh_update( spf->nh, toNode, newDist );
         } // updating heap if necessary
      } // going through node neighbors
   } // going through all nodes in priority queue
}

void spf_update_digraph( ShortestPathsFinder* spf, const int nodes, const int narcs, Arc *arcs )
{
   assert( narcs );
   spf->nodes = nodes;
   spf->arcs  = narcs;

   if (nodes > spf->capnodes )
   {
      if ( spf->startn )
         free( spf->startn );
      if ( spf->previous )
         free( spf->previous );
      if ( spf->dist )
         free( spf->dist );
      if ( spf->nh )
         nh_free( &(spf->nh) );

      spf->capnodes = nodes;
      if (nodes<1000)
         spf->capnodes = (int) ((double)nodes)*1.5;
      else
         spf->capnodes = (int) ((double)nodes)*1.1;

      spf->startn   = (Neighbor**) xmalloc( sizeof(Neighbor*)*(spf->capnodes+1) );
      spf->nh       = nh_create( spf->capnodes, SP_INFTY_DIST );
      spf->previous = (int*) xmalloc( sizeof(int)*spf->capnodes );
      spf->dist     = (int*) xmalloc( sizeof(int)*spf->capnodes );
      spf->path     = (int*) xmalloc( sizeof(int)*spf->capnodes );
   }

   if ( narcs > spf->caparcs )
   {
      if ( spf->caparcs )
         free( spf->neighs );

      spf->caparcs = narcs;
      if (narcs<1000)
         spf->caparcs = (int) ((double)narcs)*1.5;
      else
         spf->caparcs = (int) ((double)narcs)*1.1;

      spf->neighs = (Neighbor*) xmalloc( sizeof(Neighbor)*spf->caparcs );
   }

   qsort( (void *)arcs, narcs, sizeof(Arc), (void*)compArcsLX );

   if ( arcs[0].tail )  // not starting in zero ... fixing
   {
      int shift = arcs[0].tail;
      Arc *ptr = arcs;
      Arc *ptrEnd = arcs + narcs;
      for ( ; (ptr<ptrEnd) ; ptr++ )
      {
         ptr->tail -= shift;
         ptr->head -= shift;
      }
   }

   Arc *currArc        = arcs;
   Arc *endArcs        = currArc + narcs;

   //validating contents (only in debug)
#ifdef DEBUG
   for (  ; (currArc<endArcs) ; currArc++ )
   {
      assert( (currArc->head >= 0) );
      assert( (currArc->head < narcs) );
      assert( (currArc->tail >= 0) );
      assert( (currArc->tail < narcs) );
   }
#endif

   // filling default values for startn
   {
      int i;
      for ( i=0 ; (i<=nodes) ; i++ )
         spf->startn[i] = NULL;
   }

   int currNode = -1;
   Neighbor *currNeigh = spf->neighs;
   currArc = arcs;
   endArcs = arcs + narcs;
   Neighbor *lastValid = NULL;
PROCCESS_ARC:
   if ( currArc>=endArcs )
      goto END_ARCS;
   if ( currArc->tail != currNode )
   {
      currNode = currArc->tail;
      spf->startn[currNode] = currNeigh;
   }

   currNeigh->node = currArc->head;
   currNeigh->distance = currArc->distance;

   ++currArc;
   ++currNeigh;
   goto PROCCESS_ARC;
END_ARCS:

   lastValid = currNeigh;
   {
      int i;
      for ( i=nodes ; (i>=0) ; i-- )
      {
         if ( !spf->startn[i] )
            spf->startn[i] = lastValid;
         else
            lastValid = spf->startn[i];
      }
   }
}

int spf_nodes( ShortestPathsFinder* spf )
{
   return spf->nodes;
}

int spf_arcs( ShortestPathsFinder* spf )
{
   return spf->arcs;
}

Neighbor *spf_start_n( ShortestPathsFinder* spf, const int node )
{
   return spf->startn[node];
}

Neighbor *spf_end_n( ShortestPathsFinder* spf, const int node )
{
   return spf->startn[node+1];
}

int spf_get_dist( const ShortestPathsFinderPtr spf, const int node )
{
   assert( node < spf->nodes );
   return spf->dist[ node ];
}

int spf_get_previous( const ShortestPathsFinderPtr spf, const int node )
{
   assert( node < spf->nodes );
   return spf->previous[ node ];
}

int *spf_previous( const ShortestPathsFinder *spf )
{
   return spf->previous;
}

int spf_get_path( const ShortestPathsFinder *spf, const int toNode, int indexes[] )
{
   // filling first in path
   int currNode = spf->previous[toNode];
   if ( currNode == NULL_NODE )
      return 0;

   int *ptrIdx = spf->path;
   *ptrIdx = toNode;
   ++ptrIdx;
   do
   {
      *ptrIdx = currNode;
      ++ptrIdx;
   } while ( (currNode=spf->previous[currNode]) != NULL_NODE );

   const int n = ptrIdx-spf->path;
   --ptrIdx;
   {
      int i;
      for ( i=0 ; (i<n) ; ++i,--ptrIdx )
         indexes[i] = *ptrIdx;
   }

   return ( n );
}

int spf_get_path_fw( const ShortestPathsFinder *spf, const int fromNode, const int toNode, int indexes[] )
{
   // filling first in path
   int currNode = spf->fwPrev[fromNode][toNode];
   if ( currNode == NULL_NODE )
      return 0;

   int *ptrIdx = spf->path;
   *ptrIdx = toNode;
   ++ptrIdx;
   do
   {
      *ptrIdx = currNode;
      ++ptrIdx;
   } while ( (currNode=spf->fwPrev[fromNode][currNode]) != NULL_NODE );

   const int n = ptrIdx-spf->path;
   --ptrIdx;
   {
      int i;
      for ( i=0 ; (i<n) ; ++i,--ptrIdx )
         indexes[i] = *ptrIdx;
   }

   return ( n );
}

void spf_free( ShortestPathsFinderPtr *spf )
{
   freeFWSpace( *spf );

   if ( (*spf)->neighs )
      free ( (*spf)->neighs );
   if ( (*spf)->startn )
      free ( (*spf)->startn );
   if ( (*spf)->nh )
      nh_free( &((*spf)->nh) );
   if ( (*spf)->previous )
      free( (*spf)->previous );
   if ( (*spf)->dist )
      free ( (*spf)->dist );
   if ( (*spf)->path )
      free ( (*spf)->path );

   free( (*spf) );
   (*spf) = NULL;
}

void allocateFWSpace( ShortestPathsFinder* spf )
{
   if ( ( spf->nodes > (spf->fwCapNodes) ) || ( spf->arcs > (spf->fwCapArcs) ) )
   {
      freeFWSpace( spf );

      spf->fwCapNodes = spf->nodes;
      spf->fwCapArcs  = spf->arcs;

      const int cells = spf->fwCapNodes*spf->fwCapNodes;

      // allocating
      spf->fwDist = xmalloc( sizeof(int*)*spf->fwCapNodes );
      spf->fwDist[0] = xmalloc( sizeof(int)*cells );
      {
         int i;
         for ( i=1 ; (i<spf->fwCapNodes) ; ++i )
            spf->fwDist[i] = spf->fwDist[i-1]+spf->fwCapNodes;
      }

      spf->fwPrev = xmalloc( sizeof(int*)*spf->fwCapNodes );
      spf->fwPrev[0] = xmalloc( sizeof(int)*cells );
      {
         int i;
         for ( i=1 ; (i<spf->fwCapNodes) ; ++i )
            spf->fwPrev[i] = spf->fwPrev[i-1]+spf->fwCapNodes;
      }
   }

   {
      int i,j;
      for ( i=0 ; (i<spf->nodes) ; ++i )
         for ( j=0 ; (j<spf->nodes) ; ++j )
            spf->fwPrev[i][j] = -1;
   }
   
   {
      int i,j;
      for ( i=0 ; (i<spf->nodes) ; ++i )
         for ( j=0 ; (j<spf->nodes) ; ++j )
            spf->fwDist[i][j] = SP_INFTY_DIST;
   }
   
   {
      int i;
      for ( i=0 ; (i<spf->nodes) ; ++i )
         spf->fwDist[i][i] = 0;
   }

   // filling info of neighbors
   {
      int i;
      for ( i=0 ; (i<spf->nodes) ; ++i )
      {
         Neighbor *n = spf->startn[i];
         Neighbor *e = spf->startn[i+1];

         for ( ; (n<e) ; ++n )
         {
            spf->fwDist[i][n->node] = n->distance;
            spf->fwPrev[i][n->node] = i;
         }
      }
   }
}

void spf_proccessFWLoop( ShortestPathsFinder* spf )
{
   int **dist = spf->fwDist;
   int **prev = spf->fwPrev;

   const int nodes = spf->nodes;
   
   {
      int i,k,j;
      for ( k=0 ; (k<nodes) ; ++k )
         for ( i=0 ; (i<nodes) ; ++i )
         {
            const int distIK = dist[i][k];
            for ( j=0 ; (j<nodes) ; ++j )
            {
               const int distWithK = distIK + dist[k][j];
               if ( dist[i][j]>distWithK )
               {
                  dist[i][j] = distWithK;
                  prev[i][j] = k;
               }
            }

         }
   }
}

void spf_fw_find( ShortestPathsFinder* spf )
{
   allocateFWSpace( spf );

   spf_proccessFWLoop( spf );
}

int spf_fw_get_dist( ShortestPathsFinder* spf, const int i, const int j )
{
#ifdef DEBUG
   assert( spf!=NULL );
   assert( i>=0 );
   assert( j>=0 );
   assert( i<spf->nodes );
   assert( j<spf->nodes );
#endif

   return spf->fwDist[i][j];
}


int lenAlphaNumChars( char *str )
{
   int result = 0;
   char *pStr = str;
   for ( ; ( *pStr!='\0' ) ; pStr++ )
      result += isalnum(*pStr);
   return result;
}

ShortestPathsFinder *spf_load_gr( const char *fileName )
{
#define STR_SIZE 256
   int nnodes = -1, narcs = -1;

   FILE *f = fopen( fileName, "r" );
   if (!f)
   {
      fprintf( stderr, "Error reading instance.\n" );
      exit(EXIT_FAILURE);
   }

   Arc *ptrArc  = NULL;
   Arc *endArcs = NULL;
   Arc *arcs    = NULL;
   char c;
   int nreads;


   char line[STR_SIZE];

   char **hColumns;
   CREATE_STRING_VECTOR( hColumns, 10, STR_SIZE );

   while ( (fgets( &(line[0]), STR_SIZE, f )) )
   {
      int len = strlen( &(line[0]) );
      if (!len)
         continue;
      if (!lenAlphaNumChars(&(line[0])))
         continue;
      switch ( line[0] )
      {
      case 'p':
         splitString( hColumns, line, ' ', 10, STR_SIZE, 1 );
         nnodes = atoi( hColumns[2] );
         narcs = atoi( hColumns[3] );

         ptrArc = arcs = (Arc*) xmalloc( sizeof(Arc)*narcs );
         endArcs = ptrArc + narcs;

         break;
      case 'a':
         if (!ptrArc)
         {
            fprintf( stderr, "Number of arcs not defined.\n" );
            exit( EXIT_FAILURE );
         }
         if ( ptrArc >= endArcs )
         {
            fprintf( stderr, "Number of arcs wrongly defined.\n" );
            exit( EXIT_FAILURE );
         }

         nreads = sscanf( &(line[0]), "%c %d %d %d", &c, &(ptrArc->tail), &(ptrArc->head), &(ptrArc->distance)  );
         if (nreads<4)
         {
            fprintf( stderr, "Incomplete arc line.\n" );
            exit(EXIT_FAILURE);
         }
         ++ptrArc;

         break;
      default:
         continue;
      }
   }

   if (ptrArc != endArcs)
   {
      fprintf( stderr, "Not all arcs informed.\n" );
      exit(EXIT_FAILURE);
   }

   fclose(f);

   //printf("Successfully read graph with %d nodes and %d arcs.\n", nnodes, narcs );

   ShortestPathsFinder *result  = spf_create();
   spf_update_digraph( result, nnodes, narcs, arcs );
   free( arcs );

   FREE_STRING_VECTOR( hColumns );

   return result;

#undef STR_SIZE
}

void spf_update_arc( ShortestPathsFinder* spf, const int tail, const int head, const int cost )
{
   const Neighbor *start  = spf->startn[ tail ];
   const Neighbor *end    = spf->startn[ tail+1 ];
   const Neighbor key     = { head, 0 };
   Neighbor *result = bsearch( &key, start, end-start, sizeof(Neighbor), &compNeighs );
   assert( ( (result) && (result->node==head) ) );
   result->distance = cost;
}

int spf_get_arc( ShortestPathsFinder* spf, const int head, const int tail )
{
   const Neighbor *start  = spf->startn[ head ];
   const Neighbor *end    = spf->startn[ head+1 ];
   const Neighbor key     = { tail, 0 };
   const Neighbor *result = bsearch( &key, start, end-start, sizeof(Neighbor), &compNeighs );
#ifdef DEBUG
   if (!result)
   {
      fprintf( stderr, "Could not find arc (%d->%d) in graph.\n", head, tail );
      abort();
      exit(EXIT_FAILURE);
   }
   if ( result->node )
   assert( ( (result) && (result->node==tail) ) );
#endif
   return result->distance;
}

void spf_temp_remove_arc( ShortestPathsFinder* spf, const int tail, const int head )
{
   spf_update_arc( spf, tail, head, SP_INFTY_DIST );
}

void spf_restore_arc( ShortestPathsFinder* spf, const int tail, const int head )
{
   spf_update_arc( spf, tail, head, SP_INFTY_DIST );
}

void spf_update_graph( ShortestPathsFinder* spf, const int nodes, const int arcs, const int *arcStart, const int *toNode, const int *dist )
{
   //printf("i %d %d %d\n", nodes, arcs, spf->capnodes ); fflush(stdout);
   spf->nodes = nodes;
   spf->arcs  = arcs;

   if (nodes > spf->capnodes )
   {
      if ( spf->startn )
         free( spf->startn );
      if ( spf->previous )
         free( spf->previous );
      if ( spf->dist )
         free( spf->dist );
      if ( spf->nh )
         nh_free( &(spf->nh) );

      spf->capnodes = nodes;
      if (nodes<1024)
         spf->capnodes = (int) ((double)nodes)*1.5;
      else
         spf->capnodes = (int) ((double)nodes)*1.25;

      spf->startn   = (Neighbor**) xmalloc( sizeof(Neighbor*)*(spf->capnodes+1) );
      spf->nh       = nh_create( spf->capnodes, SP_INFTY_DIST );
      spf->previous = (int*) xmalloc( sizeof(int)*spf->capnodes );
      spf->dist     = (int*) xmalloc( sizeof(int)*spf->capnodes );
      spf->path     = (int*) xmalloc( sizeof(int)*spf->capnodes );
   }

   if ( arcs > spf->caparcs )
   {
      if ( spf->caparcs )
         free( spf->neighs );

      spf->caparcs = arcs;
      if (arcs<1000)
         spf->caparcs = (int) ((double)arcs)*1.5;
      else
         spf->caparcs = (int) ((double)arcs)*1.1;

      spf->neighs = (Neighbor*) xmalloc( sizeof(Neighbor)*spf->caparcs );
   }
   
   {
      int n;
      for ( n=0 ; (n<=spf->nodes) ; ++n )
         spf->startn[n] = spf->neighs + arcStart[n];
   }

   Neighbor *ptrNeigh = spf->neighs;
   int idx = 0;
   const Neighbor *ptrEndNeigh = spf->neighs + arcs;
   for ( ; (ptrNeigh<ptrEndNeigh) ; ++ptrNeigh,++idx )
   {
      ptrNeigh->node     = toNode[ idx ];
      assert( ptrNeigh->node<spf->nodes );
      ptrNeigh->distance = dist[ idx ];
   }
}

int spf_fw_ran( ShortestPathsFinder* spf )
{
   return (spf->fwDist!=NULL);

}

void freeFWSpace( ShortestPathsFinder* spf )
{
   if ( spf->fwCapNodes )
   {
      spf->fwCapNodes = 0;
      spf->fwCapArcs  = 0;
      free ( spf->fwDist[0] );
      free ( spf->fwDist );
      free ( spf->fwPrev[0] );
      free ( spf->fwPrev );
   }
}
