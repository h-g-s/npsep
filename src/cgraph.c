#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "cgraph.h"
#include "memory.h"
#include "vint_set.h"
#include "vectormgm.h"
#include "clique.h"
#include "macros.h"

#define INI_CLQ_CAP 512

#define MAX_NAME_SIZE 64

typedef struct
{
   int node;
   unsigned long int degree;
} NodeDegree;

int cmp_by_node_degree( const void *v1, const void *v2 );

struct _CGraph
{
   /* per node stored conflicts
      at position i there is a set with
      all node conflicts */
   IntSet *nodeConflicts;
   IntSet *nodeCliques;  /* all cliques in which a node appears */
   int nodeSize;         /* number of nodes considered */

   CliqueSet *clqSet;

   /* degree of each node */
   unsigned long int *degree;

   unsigned long int minDegree;
   unsigned long int maxDegree;

   unsigned long int lowDegree;

   int *iv;  /* incidence vector */

   int *origIdx; /* if it is a preprocessed graph,
                    indicates for each node i its original node */

   /* node names and weights are optional information */
   char **nodeNames;

   int *w;

   /* private data, used for processing */
   NodeDegree *nd;
};

typedef struct
{
   int node;
   int cost;
} NodeCost;

struct _NeighIterator
{
   /* binary heap vector */
   NodeCost *bhv;
   int bhvSize;
   int bhvCap;
   int node;

   /* previously returned node */
   int previousNode;
};

typedef struct CGArc
{
   int tail;
   int head;
} CGArc;

void cgraph_check_node_size( CGraph *cgraph, int requiredNodeSize );
void cgraph_check_clique_size( CGraph *cgraph, int requiredCliqueSize );
void cgraph_add_clique_as_normal_conflicts( CGraph *cgraph, const int nEl, const int elClique[] );

/** for every tail add all heads as neighbors. vector arcs should have capacit of nArcs+1 **/
void cgraph_add_vector_arcs( CGraph *cgraph, const CGArc *arcs, const int nArcs, int neighs[] );

CGraph *cgraph_preprocess( const CGraph *cgraph, int nindexes[] )
{
   const int gSize = cgraph_size( cgraph );
   /* pre processed graph size */
   int ppgSize = 0;

   {
      int i;
      for ( i=0 ; (i<gSize) ; ++i )
         nindexes[i] = -1;
   }

   int nConflicts;
   const int confSpace = gSize*100;
   int *conflicts = malloc( sizeof(int)*confSpace );

   /* checking how many columns will remain */
   {
      int j;
      for ( j=0 ; ( j<gSize ) ; ++j )
      {
         nConflicts = cgraph_get_all_conflicting( cgraph, j, conflicts, confSpace );
         if ( nConflicts<3 )
            continue;

         nindexes[ j ] = ppgSize++;
      }
   }

   CGraph *result = NULL;
   result = cgraph_create( ppgSize );

   /* adding new conflicts */
   {
      int j;
      for ( j=0 ; (j<gSize) ; ++j )
      {
         const int newNode = nindexes[j];
         if ( newNode == -1 )
            continue;

         int nConflicts = 0;

         const int *gConf = vint_set_get_elements( cgraph->nodeConflicts+j );
         const int gConfSize = vint_set_size( cgraph->nodeConflicts+j );

         {
            int k;
            for ( k=0 ; (k<gConfSize) ; ++k )
            {
               const int newNodeConf = nindexes[ gConf[k] ];
               if ( newNodeConf == -1 )
                  continue;
               conflicts[nConflicts++] = newNodeConf;
            }
         }

         vint_set_add( result->nodeConflicts+newNode, conflicts, nConflicts );
      }
   }

   /* adding new cliques */
   {
      int i;

      for ( i=0 ; (i<clq_set_number_of_cliques(cgraph->clqSet)) ; ++i )
      {
         nConflicts = 0;

         const int nEl = clq_set_clique_size( cgraph->clqSet,i );
         const int *elClique = clq_set_clique_elements( cgraph->clqSet,i );
         {
            int j;
            for ( j=0 ; (j<nEl) ; ++j )
            {
               const int newNode = nindexes[ elClique[j] ];
               if ( newNode == -1 )
                  continue;

               conflicts[nConflicts++] = newNode;
            }
         }

         if ( nConflicts < 2 )
            continue;

         if ( nConflicts < 50 )
         {
            /* not adding as clique anymore */
            cgraph_add_clique_as_normal_conflicts( result, nEl, elClique );

            continue;
         }

         /* large clique, adding as it is */

         cgraph_add_clique( result, conflicts, nConflicts );
      }
   }

#ifdef DEBUG
   assert( result != NULL );
   assert( cgraph_size(result) == ppgSize );
   {
      int i;
      for ( i=0 ; (i<cgraph_size(result)) ; ++i )
      {
         nConflicts = cgraph_get_all_conflicting( result, i, conflicts, confSpace );
         assert( nConflicts >= 3 );
         int j;
         for ( j=0 ; (j<nConflicts) ; ++j )
         {
            assert( conflicts[j] < ppgSize );
            assert( conflicts[j] >= 0 );
         }
         assert( nConflicts >= 3);
      }
   }
#endif

   free(conflicts);

   cgraph_update_min_max_degree( result );

   return result;
}

int cgraph_cmp_int( const void *e1, const void *e2 )
{
   int i1 = (*((const int*) e1));
   int i2 = (*((const int*) e2));

   return (i1-i2);
}

CGraph *cgraph_create( int initialColumns )
{
   CGraph *result = xmalloc( sizeof(CGraph) );

   result->nodeSize = initialColumns;
   result->nodeConflicts = xmalloc( sizeof(IntSet)*result->nodeSize );
   result->nodeCliques   = xmalloc( sizeof(IntSet)*result->nodeSize );
   {
      int i;
      for ( i=0 ; (i<result->nodeSize) ; ++i )
         vint_set_init( result->nodeConflicts+i );
      for ( i=0 ; (i<result->nodeSize) ; ++i )
         vint_set_init( result->nodeCliques+i );
   }

   result->degree = xmalloc( sizeof(unsigned long int)*initialColumns );
   memset( result->degree, 0, sizeof(unsigned long int)*initialColumns );
   result->iv = xmalloc( sizeof(int)*initialColumns );
   memset( result->iv, 0, sizeof(int)*initialColumns );

   result->nd = xmalloc( sizeof(NodeDegree)*initialColumns );

   result->clqSet = clq_set_create();

   result->nodeNames = NULL;

   result->w = NULL;

   result->origIdx = NULL;

   result->lowDegree = 15;

   return result;
}

void cgraph_add_node_conflicts( CGraph *cgraph, const int node, const int conflicts[], const int size )
{
   /* checking node capacity */
   int lastNode = node;
   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
         lastNode = MAX( lastNode, conflicts[i] );
   }
   ++lastNode;
   cgraph_check_node_size( cgraph, lastNode );

   /* increasing degree if necessary */
   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
      {
         if (!cgraph_conflicting_nodes( cgraph, node, conflicts[i] ))
            cgraph->degree[node]++;
         if (!cgraph_conflicting_nodes( cgraph, conflicts[i], node ))
            cgraph->degree[conflicts[i]]++;
      }
   }
   vint_set_add( cgraph->nodeConflicts+node, conflicts, size );

   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
         vint_set_add( cgraph->nodeConflicts+conflicts[i], &node , 1 );
   }
}

void cgraph_add_node_conflicts_no_sim( CGraph *cgraph, const int node, const int conflicts[], const int size )
{
   int lastNode = node;
   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
         lastNode = MAX( lastNode, conflicts[i] );
   }
   ++lastNode;

   cgraph_check_node_size( cgraph, lastNode );

   cgraph->degree[node] += size;
   /*
      {
      int i;
      for ( i=0 ; (i<size) ; ++i )
      if (!cgraph_conflicting_nodes( cgraph, node, conflicts[i] ))
      cgraph->degree[node]++;
      }*/

   if (size>0)
      vint_set_add_opt( cgraph->nodeConflicts+node, conflicts, size );
}

void cgraph_add_clique( CGraph *cgraph, int clique[], const int size )
{
   /* checking node capacity */
   int lastNode = -1;
   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
         lastNode = MAX( lastNode, clique[i] );
   }
   ++lastNode;

   cgraph_check_node_size( cgraph, lastNode );

   if (!(clq_set_add( cgraph->clqSet, size, clique, 0 )))
      return;

   /* updating degrees if necessary
      const int sizeM1 = size - 1;
      for ( int i=0 ; (i<sizeM1) ; ++i )
      {
      for ( int j=i+1 ; (j<size) ; ++j )
      {
      if (!cgraph_conflicting_nodes( cgraph, clique[i], clique[j] ) )
      cgraph->degree[clique[i]]++;
      if (!cgraph_conflicting_nodes( cgraph, clique[j], clique[i] ) )
      cgraph->degree[clique[j]]++;
      }
      } */

   /* aproximated degree */
   {
      int i;
      for ( i=0 ; (i<size) ; ++i )
         cgraph->degree[ clique[i] ] += size-1;
   }

   /* this clique will be related to all nodes inside it */
   {
      int i;
      const int idxClique = clq_set_number_of_cliques(cgraph->clqSet) - 1;
      for ( i=0 ; (i<size) ; ++i )
      {
         cgraph_check_node_size( cgraph, clique[i]+1 );
         vint_set_add( cgraph->nodeCliques + clique[i], &idxClique, 1 );
      }
   }
}

int cgraph_size( const CGraph *cgraph )
{
   return cgraph->nodeSize;
}

void cgraph_check_node_size( CGraph *cgraph, int requiredNodeSize )
{
   if ( cgraph->nodeSize >= requiredNodeSize )
      return;

   int newNodeSize = (int)(((double)cgraph->nodeSize) * 1.5);
   if ( newNodeSize < requiredNodeSize )
      newNodeSize = requiredNodeSize;

   /* printf(">>>> node size changed to %d\n", cgraph->nodeSize ); */

   cgraph->nodeConflicts = xrealloc( cgraph->nodeConflicts, sizeof(IntSet)*newNodeSize );
   cgraph->nodeCliques = xrealloc( cgraph->nodeCliques, sizeof(IntSet)*newNodeSize );
   cgraph->degree = xrealloc( cgraph->degree, sizeof(unsigned long int)*newNodeSize );
   cgraph->nd = xrealloc( cgraph->nd, sizeof(NodeDegree)*newNodeSize );
   memset( cgraph->degree+cgraph->nodeSize, 0, newNodeSize-cgraph->nodeSize );
   cgraph->iv = xrealloc( cgraph->iv, sizeof(int)*newNodeSize );
   memset( cgraph->iv+cgraph->nodeSize, 0, newNodeSize-cgraph->nodeSize );
   if (cgraph->w)
      cgraph->w = xrealloc( cgraph->w, sizeof(int)*newNodeSize );
   if (cgraph->nodeNames)
   {
      cgraph->nodeNames = xrealloc( cgraph->nodeNames, sizeof(char*)*newNodeSize );
      cgraph->nodeNames[0] = xrealloc( cgraph->nodeNames[0], sizeof(char)*newNodeSize*MAX_NAME_SIZE );
      int i;
      for ( i=1 ; (i<newNodeSize); ++i )
         cgraph->nodeNames[i] = cgraph->nodeNames[i-1] + MAX_NAME_SIZE;
   }

   {
      int i;
      for ( i=cgraph->nodeSize ; (i<newNodeSize) ; ++i )
         vint_set_init( cgraph->nodeConflicts+i );
      for ( i=cgraph->nodeSize ; (i<newNodeSize) ; ++i )
         vint_set_init( cgraph->nodeCliques+i );
   }

   cgraph->nodeSize = newNodeSize;
}

int cgraph_conflicting_nodes( const CGraph *cgraph, const int i, const int j )
{
   const IntSet *cliques;
#ifdef DEBUG
   assert( i>=0 );
   assert( j>=0 );
   assert( i<cgraph_size(cgraph) );
#endif
   if ( vint_set_size( cgraph->nodeConflicts+i ) == 0 )
      goto FIND_IN_CLIQUES;

   if ( vint_set_find( cgraph->nodeConflicts+i, j ) )
      return 1;

FIND_IN_CLIQUES:
   cliques = cgraph->nodeCliques+i;

   const int nCliques = vint_set_size( cliques );
   const int *idxCliques = vint_set_get_elements( cliques );

   {
      int idx;
      for ( idx=0 ; (idx<nCliques ) ; ++idx )
         if ( clq_set_clique_has_element( cgraph->clqSet, idxCliques[idx], j ) )
         {
#ifdef DEBUG
            if (!clq_set_clique_has_element( cgraph->clqSet, idxCliques[idx], i ))
            {
               fprintf( stderr, "Error:  node %d does not appears in clique %d\n", i, idxCliques[idx] );
            }
            assert( clq_set_clique_has_element( cgraph->clqSet, idxCliques[idx], i ) );
            assert( clq_set_clique_has_element( cgraph->clqSet, idxCliques[idx], j ) );
#endif
            return 1;
         }
   }

   return 0;
}

int cgraph_get_all_conflicting( const CGraph *cgraph, int node, int neighs[], int maxSize )
{
#ifdef DEBUG
   assert( cgraph!=NULL );
   assert( (node>=0) && (node<cgraph_size(cgraph)) );
   assert( maxSize>0 );
#endif

   const IntSet *isnc = cgraph->nodeConflicts+node;
   int result = isnc->size;
   int currClique=0, clqSize=0;

   if ( result > maxSize )
      goto NO_SPACE;

   memcpy( neighs, vint_set_get_elements(isnc), sizeof(int)*vint_set_size(isnc) );

   /* now filling information from cliques, i.e., implicitly stored conflicts */
   const int nCliques = vint_set_size( cgraph->nodeCliques+node );
   const int *idxCliques = vint_set_get_elements( cgraph->nodeCliques+node );

   {
      for ( currClique=0 ; (currClique<nCliques) ; ++currClique )
      {
         clqSize = clq_set_clique_size( cgraph->clqSet, idxCliques[currClique] );
         const int *clqElements = clq_set_clique_elements( cgraph->clqSet, idxCliques[currClique] );

         int *ptr = neighs + result;
         result += clqSize-1;
         if ( result > maxSize )
            goto NO_SPACE;
         
         /* copying all elements except itself */
         {
            int j=0,i=0;
            for ( ; (i<clqSize) ; i++ )
            {
               if (clqElements[i] == node) 
                  continue;
               ptr[j++] = clqElements[i];
            }
         }
         
         qsort( neighs, result, sizeof(int), cgraph_cmp_int );
         
         /* conflicts in this clique may be already added, removing */
         int rep = 0;
         {
            int j;
            for ( j=1 ; (j<result) ; ++j )
               if ( neighs[j] == neighs[j-1] )
               {
                  neighs[j] = INT_MAX;
                  ++rep;
               }
         }
         
         qsort( neighs, result, sizeof(int), cgraph_cmp_int );

         result -= rep;
      }
   }

#ifdef DEBUG
   {
      assert( result >= 0);
      assert( result < cgraph_size(cgraph) );
      assert( result <= maxSize );
      int i;
      /* all inside conflict with node */
      for ( i=0 ; (i<result) ; i++ ) 
      {
         assert( neighs[i]>=0 );
         assert( neighs[i]<cgraph_size(cgraph) );
         assert( neighs[i]!=node );
         assert( cgraph_conflicting_nodes(cgraph, node, neighs[i] ) );
      }
      /* all which conflict with node are inside */
      for ( i=0 ; (i<cgraph_size(cgraph)) ; i++ ) 
      {
         if ( (i!=node) && (cgraph_conflicting_nodes( cgraph, node, i )) )
            assert( bsearch( &i, neighs, result, sizeof(int), cgraph_cmp_int ) != NULL );
      }
    }
#endif
   
   return result;

NO_SPACE:
   fprintf( stderr, "ERROR: cgraph_get_all_conflicting:: Not enough space specified in maxSize.\n" );
   fprintf( stderr, "Working with node %d, which appears in %d cliques. When adding clique %d size %d. Result %d. MaxSize %d.\n", node, vint_set_size(cgraph->nodeCliques+node), currClique, clqSize, result, maxSize );
   fprintf( stderr, "at: %s:%d\n", __FILE__, __LINE__ );
   abort();
   exit( EXIT_FAILURE );
}

int cgraph_get_candidates_clique_insertion( CGraph *cgraph, const int n, const int clique[], int candidates[], const int capacityCandidates )
{
   NodeDegree *nd = cgraph->nd;

   /* clique sorted by degree (smallest to greatest) */
   {
      int i;
      for ( i=0 ; (i<n) ; ++i )
      {
         nd[i].node = clique[i];
         nd[i].degree = cgraph_degree( cgraph, clique[i] );
      }
   }

   qsort( nd, n, sizeof(NodeDegree), cmp_by_node_degree );

#ifdef DEBUG
   {
      /* checking if sorting is ok */
      int i;
      for ( i=0 ; (i<(n-1)) ; ++i )
      {
         assert( clique[i]<=clique[i+1] );
         assert( nd[i].degree<=nd[i+1].degree );
      }
   }
#endif

   int nCandidates = cgraph_get_all_conflicting( cgraph, nd[0].node, candidates, capacityCandidates );

   /* only conflicts which are not inserted will be considered */
   {
      int removals = 0;
      int *cand = candidates;
      const int *endCand = candidates + nCandidates;
      for ( ; (cand<endCand) ; ++cand )
      {
         if (bsearch( cand, clique, n, sizeof(int), vint_set_cmp_int ))
         {
            ++removals;
            *cand = INT_MAX;
         }
      }

      if ( removals )
      {
         qsort( candidates, nCandidates, sizeof(int), vint_set_cmp_int );
         nCandidates -= removals;
      }
   }

#ifdef DEBUG
   assert( nCandidates >= 0 ); /* non negativity */
#endif

   /* now we have all conflicts of the node with the smallest degree which are not
      inserted in the clique. checking how many of these can expand the clique.
      checking if they have conflict with other clique. */

   if ( !nCandidates )
      return 0;

   {
      int i, removals = 0;
      for ( i=0 ; (i<nCandidates) ; ++i )
      {
         const int node = candidates[i];
         int j;
         for ( j=1 ; (j<n) ; ++j )
         {
            if (!cgraph_conflicting_nodes(cgraph, node, nd[j].node))
            {
               candidates[i] = INT_MAX;
               ++removals;
               break;
            }
         } /* all other clique already in the clique */
      } /* all candidates for insertion */

      if ( removals )
      {
         qsort( candidates, nCandidates, sizeof(int), vint_set_cmp_int );
         nCandidates -= removals;
      }
   }

#ifdef DEBUG
   assert( nCandidates >= 0 ); /* non negativity */
#endif

   return nCandidates;
}

void cgraph_save( CGraph *cgraph, const char *fileName )
{
   int *w = cgraph->w;

   FILE *f = fopen( fileName, "w" );
   if (!f)
   {
      fprintf( stderr, "Could not open file %s.", &(fileName[0]) );
      exit( EXIT_FAILURE );
   }

   /* counting edges */
   int nEdges = 0;
   int nodes = cgraph_size( cgraph );

   int *degree = xmalloc( sizeof(unsigned long int)*nodes );
   memset( degree, 0, sizeof(unsigned long int)*nodes );

   {
      int i,j;
      for ( i=0 ; (i<nodes-1) ; ++i )
         for ( j=i+1 ; (j<nodes) ; ++j )
            if ( cgraph_conflicting_nodes( cgraph, i, j ) )
            {
               degree[i]++;
               degree[j]++;
            }
      for ( i=0 ; (i<nodes-1) ; ++i )
         for ( j=i+1 ; (j<nodes) ; ++j )
            if ( ( cgraph_conflicting_nodes( cgraph, i, j ) )  )
               nEdges++;

   }


   fprintf( f, "p edges %d %d\n", cgraph_size( cgraph ), nEdges );

   {
      int i, j;
      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      {
         const IntSet *is = cgraph->nodeConflicts + i;
         const int nsize = vint_set_size( is );
         const int *neighs = vint_set_get_elements( is );
         for ( j=0 ; (j<nsize) ; ++j )
            if ( neighs[j] > i )
               fprintf( f, "e %d %d\n", i+1, neighs[j]+1 );
      }
   }

   if (w)
   {
      {
         int i;
         for ( i=0 ; (i<nodes) ; ++i )
            fprintf( f, "w %d %d\n", i+1, w[i] );
      }
   }

   {
      int i,j;
      for ( i=0 ; (i<clq_set_number_of_cliques( cgraph->clqSet)) ; ++i )
      {
         const int size = clq_set_clique_size( cgraph->clqSet, i );
         const int *elements = clq_set_clique_elements( cgraph->clqSet, i );
         fprintf( f, "c %d\n", size );
         for ( j=0 ; (j<size ) ; ++j )
            fprintf( f, "%d\n", elements[j]+1 );
      }
   }

   if ( cgraph->nodeNames )
   {
      int i;
      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      {
         if ( cgraph_get_node_name( cgraph, i ) )
         {
            char name[256];
            strcpy( name, cgraph_get_node_name( cgraph, i ) );
            if ( name[ strlen(name)-1 ] == '\n' )
               name[ strlen(name)-1 ] = '\0';
            fprintf( f, "n %d %s\n", i+1, name );
         }
      }
   }

   if ( cgraph->origIdx )
   {
      int i;
      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
         fprintf( f, "o %d %d\n", i+1, cgraph->origIdx[i]+1 );
   }

   free( degree );

   fclose( f );
}

int compArc( const void *p1, const void *p2 )
{
   if ( ((const CGArc*)p1)->tail != ((const CGArc*)p2)->tail )
      return ((const CGArc*)p1)->tail - ((const CGArc*)p2)->tail;

   return ((const CGArc*)p1)->head - ((const CGArc*)p2)->head;
}

void cgraph_add_vector_arcs( CGraph *cgraph, const CGArc *arcs, const int nArcs,int neighs[] )
{
   if ( nArcs <= 0 )
      return;

   neighs[ 0 ] = arcs[0].head;
   int nNeighs = 1;
   int lastNode = arcs[0].tail;
   const int nArcsP1 = nArcs+1;

   {
      int i;
      for ( i=1 ; (i<(nArcsP1)) ; ++i )
      {
         if ( lastNode != arcs[i].tail )
         {
            assert( lastNode >= 0 ); assert( lastNode < cgraph_size( cgraph ) );
            cgraph_add_node_conflicts_no_sim( cgraph, lastNode, neighs, nNeighs );
            lastNode = arcs[i].tail;
            nNeighs = 0;
            if ( i == nArcs )
               break;
         }

         assert( arcs[ i ].head >= 0 ); assert( arcs[ i ].head < cgraph_size( cgraph ) );
         neighs[nNeighs++] = arcs[ i ].head;
      }
   }
}

CGraph *cgraph_load( const char *fileName )
{
#define LINE_SIZE 256
#define STR_SIZE 128
   int nodes = -1;
   int edges = -1;

   FILE *f = fopen( fileName, "r" );
   if (!f)
   {
      fprintf( stderr, "Could not open file %s.\n", fileName );
         exit( EXIT_FAILURE );
   }

   char line[LINE_SIZE];
   while ( fgets( line, LINE_SIZE, f ) )
   {
      if ( (line[0] == 'p') ||  (line[0] == 'P') )
      {
         char *ptr = line;
         char *end = ptr + LINE_SIZE;
         while ( (!isdigit(*ptr)) && (ptr<end) )
            ++ptr;
         sscanf( ptr, "%d %d", &nodes, &edges );
         break;
      }
   }

   if (nodes==-1)
   {
      fprintf( stderr, "Invalid file format.\n" );
      exit( EXIT_FAILURE );
   }

   CGArc *arcs = NULL;
   int arcCap = 8192;
   int nArcs = 0;
   arcs = xmalloc( sizeof(CGArc)*(arcCap+1) );
   int *clqElement = xmalloc( sizeof(int)*nodes );
   CliqueSet *clqSet = clq_set_create();
   int weightsLoaded = 0; /* lineNumber = 1; */

   char **names = NULL;
   int *tw = NULL;

   while ( fgets( line, LINE_SIZE, f ) )
   {
      if ( (line[0] == 'e') || (line[0] == 'E') )
      {
         int tail, head;
         sscanf( line+2, "%d %d", &tail, &head );
         --tail; --head;
         assert( tail >= 0 ); assert( tail < nodes );
         assert( head >= 0 ); assert( head < nodes );
         if ( tail > head )
         {
            int t = tail;
            tail = head;
            head = t;
         }

         /*         ADJUST_VECTOR_CAPACITY( arcs, arcCap, nArcs+2 );   */
         vmg_adjust_vector_capacity( (void**) (&arcs), &arcCap, nArcs+2, sizeof(CGArc) );

         arcs[nArcs].tail = tail;
         arcs[nArcs].head = head;
         nArcs++;
      }
      else
      {
         if ( ( (line[0] == 'w') || (line[0] == 'W') ) )
         {
            if (!tw)
            {
               tw = xmalloc( sizeof(int)*nodes );
            }

            int nodeIdx, nodeWeight;
            sscanf( line+2, "%d %d", &nodeIdx, &nodeWeight );
            nodeIdx--;

            assert( ( nodeIdx >=0 ) && ( nodeIdx < nodes ) );
            tw[nodeIdx] = nodeWeight;

            weightsLoaded = 1;
         }
         else
         {
            if ( ( (line[0] == 'c') || (line[0] == 'C') ) )
            {
               assert( weightsLoaded );
               int clqSize, totalWeight = 0;
               sscanf( line+2, "%d", &clqSize );
               {
                  int i;
                  for ( i=0 ; (i<clqSize) ; ++i )
                  {
                     if ( fscanf( f, "%d", clqElement+i ) != 1 )
                     {
                        fprintf( stderr, "\nERROR: Expected clique element not found.\n" );
                        exit( EXIT_FAILURE );
                     }
                     clqElement[i]--;
                     assert( clqElement[i]>=0 ); assert( clqElement[i]<nodes );
                     totalWeight += tw[ clqElement[i] ];
                  }
               }
               clq_set_add( clqSet, clqSize, clqElement, totalWeight );
            } /* reading clique */
            else
            {
               /* reading column names */
               if ( ( (line[0] == 'n') || (line[0] == 'N') ) )
               {
                  if ( names == NULL )
                  {
                     names = xmalloc( sizeof(char*)*nodes );
                     names[0] = xmalloc( sizeof(char)*MAX_NAME_SIZE*nodes );
                     int i;
                     for ( i=1 ; (i<nodes) ; ++i )
                        names[i] = names[i-1] + MAX_NAME_SIZE;
                  }

                  int idxNode;
                  sscanf( line+2, "%d", &idxNode );
                  int digits = 1;
                  idxNode--;
                  if ( idxNode != 0 )
                     digits = (int) ( log10( ((float)idxNode) ) + 1.0 );

                  strncpy( names[idxNode], line + 2 + digits + 1, MAX_NAME_SIZE );
                  int len = strlen( names[idxNode] );
                  if ( len > 0 )
                     if ( names[idxNode][len-1]=='\n' )
                        names[idxNode][len-1] = '\0';
               }
            } /* not clique */
         } /* not weight */
      } /* not edge */
   }

   qsort( arcs, nArcs, sizeof(CGArc), compArc );

   /* removing repetitions */
   int nRep = 0;
   CGArc lastArc = { -1, -1 };
   {
      int i;
      for ( i=0 ; (i<nArcs) ; ++i )
         if ( compArc( arcs+i, &lastArc ) == 0 )
         {
            arcs[i].tail = INT_MAX;
            arcs[i].head = INT_MAX;
            nRep++;
         }
         else
            lastArc = arcs[i];
   }

   if (nRep)
      qsort( arcs, nArcs, sizeof(CGArc), compArc );

   nArcs -= nRep;

   CGraph *cgraph = cgraph_create( nodes );

   if (names)
      cgraph->nodeNames = names;
   if (tw)
   {
      cgraph->w = xmalloc( sizeof(int)*nodes );
      memcpy( cgraph->w, tw, sizeof(int)*nodes );
      free( tw );
   }

   int *neighs = malloc( sizeof(int)*nArcs );

   arcs[nArcs].tail = INT_MAX;
   arcs[nArcs].head = INT_MAX;

   cgraph_add_vector_arcs( cgraph, arcs, nArcs, neighs );

   /* inverting head and tail */
   {
      int i;
      for ( i=0 ; (i<nArcs) ; ++i )
      {
         const int t = arcs[i].tail;
         arcs[i].tail = arcs[i].head;
         arcs[i].head = t;
      }
   }

   qsort( arcs, nArcs, sizeof(CGArc), compArc );

   cgraph_add_vector_arcs( cgraph, arcs, nArcs, neighs );

   {
      int *t = xmalloc( sizeof(int)*nodes );
      int i;
      for ( i=0 ; (i<clq_set_number_of_cliques(clqSet)) ; ++i )
      {
         const int size = clq_set_clique_size(clqSet, i);
         const int *el = clq_set_clique_elements(clqSet, i);
         memcpy( t, el, sizeof(int)*size );
         cgraph_add_clique( cgraph, t, size );
      }
      free(t);
   }

   fclose( f );
   free( arcs );
   free( neighs );
   free(clqElement);
   clq_set_free( &clqSet );

   cgraph_update_min_max_degree( cgraph );

   return cgraph;
#undef STR_SIZE
#undef LINE_SIZE
}

void cgraph_print( CGraph *cgraph, const int w[] )
{
   int *neighs = xmalloc( sizeof(int)*cgraph_size(cgraph)*100  );
   {
      int i;
      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      {
         printf("[%d] ", i+1 );
         int n = cgraph_get_all_conflicting( cgraph, i, neighs, cgraph_size(cgraph)*100 );
         {
            int j;
            for ( j=0 ; (j<n) ; ++j )
               printf("%d ", neighs[j]+1 );
            printf("\n");
         }
      }
   }
   if (w)
   {
      int i;
      for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
         printf( "w[%d] %d\n", i+1, w[i] );
   }

   free( neighs );
}

void cgraph_update_min_max_degree( CGraph *cgraph )
{
   int i;

   cgraph->minDegree = INT_MAX;
   for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      cgraph->minDegree = MIN( cgraph->minDegree, cgraph->degree[i] );

   cgraph->maxDegree = 0;
   for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      cgraph->maxDegree = MAX( cgraph->maxDegree, cgraph->degree[i] );
}

int cgraph_degree( const CGraph *cgraph, const int node )
{
   return cgraph->degree[ node ];
}

int cgraph_min_degree( const CGraph *cgraph )
{
   return cgraph->minDegree;
}

int cgraph_max_degree( const CGraph *cgraph )
{
   return cgraph->maxDegree;
}

void cgraph_free( CGraph **cgraph )
{
   int i;
   for ( i=0 ; (i<(*cgraph)->nodeSize) ; ++i )
      vint_set_clean( (*cgraph)->nodeConflicts+i );

   for ( i=0 ; (i<(*cgraph)->nodeSize) ; ++i )
      vint_set_clean( (*cgraph)->nodeCliques+i );

   free( (*cgraph)->nodeConflicts );
   free( (*cgraph)->nodeCliques );
   free( (*cgraph)->nd );

   clq_set_free( &((*cgraph)->clqSet) );

   free( (*cgraph)->degree );
   free( (*cgraph)->iv );

   if ((*cgraph)->origIdx)
      free( (*cgraph)->origIdx );

   if ((*cgraph)->nodeNames)
   {
      free((*cgraph)->nodeNames[0]);
      free((*cgraph)->nodeNames);
      (*cgraph)->nodeNames = NULL;
   }

   if ((*cgraph)->w)
      free((*cgraph)->w);

   free( (*cgraph) );

   (*cgraph) = NULL;
}

void cgraph_load_dimensions( const char *fileName, int *nodes, int *edges )
{
#define LINE_SIZE 64
   FILE *f = fopen( fileName, "r" );
   if (!f)
   {
      fprintf( stderr, "Could not open file %s\n", &(fileName[0]) );
         exit( EXIT_FAILURE );
   }

   *nodes = -1;
   *edges = -1;

   char line[LINE_SIZE];
   while ( fgets( line, LINE_SIZE, f ) )
   {
      if ( (line[0] == 'p') || (line[0] == 'P') )
      {
         const char *p;
         if ( ((p = strstr( line, "edge" ))!=NULL) )
         {
            int shift = 7;
            if ( p[5] == 's' )
               ++shift;
            if ( sscanf( line+shift, "%d %d", nodes, edges ) != 2 )
            {
               fprintf( stderr, "Error reading dimensions.\n" );
               exit( EXIT_FAILURE );
            }
         }

      }
   }

   fclose( f );
#undef LINE_SIZE
}

CGraph *cgraph_create_induced_subgraph( const CGraph *cgraph, const int nindexes[] )
{
   const int n = cgraph_size( cgraph );

   int newN = 0;

   int i;
   for ( i=0 ; (i<n) ; ++i )
      if ( nindexes[i] != -1 )
         ++newN;


   CGraph *result = cgraph_create( newN );
   if ( result->origIdx )
      free( result->origIdx );
   result->origIdx = xmalloc( sizeof(int)*newN );

   {
      int last = 0;
      for ( i=0 ; (i<n) ; ++i )
         if ( nindexes[i] != -1 )
            result->origIdx[ last++ ] = i;
   }

   int *neighs = xmalloc( sizeof(int)*newN );


   if (cgraph->w)
      result->w = xmalloc( sizeof(int)*newN );
   if (cgraph->nodeNames)
   {
      result->nodeNames = xmalloc( sizeof(char*)*cgraph_size(result) );
      result->nodeNames[0] = xmalloc( sizeof(char)*MAX_NAME_SIZE*cgraph_size(result) );
      int i;
      for ( i=1 ; (i<cgraph_size(result)) ; ++i )
         result->nodeNames[i] = result->nodeNames[i-1] + MAX_NAME_SIZE;
   }

   /* filling new conflicts */
   for ( i=0 ; (i<n) ; ++i )
   {
      if ( nindexes[i] == -1 )
         continue;

      int nNeighs = 0;

      const IntSet *nodeConflicts = cgraph->nodeConflicts+i;
      const int *c = vint_set_get_elements( nodeConflicts );
      const int nc = vint_set_size( nodeConflicts );

      int j;
      for ( j=0 ; (j<nc) ; ++j )
      {
         if ( nindexes[ c[j] ] != -1 )
            neighs[ nNeighs++ ] = nindexes[ c[j] ];
      }

      cgraph_add_node_conflicts_no_sim( result, nindexes[i], neighs, nNeighs );

      if ((cgraph->w)&&(result->w))
         result->w[nindexes[i]] = cgraph->w[i];
      if ((cgraph->nodeNames)&&(result->nodeNames))
         strncpy( result->nodeNames[nindexes[i]], cgraph->nodeNames[i], MAX_NAME_SIZE );
   }

   /* filling new cliques */
   const int nCliques = clq_set_number_of_cliques( cgraph->clqSet );
   for ( i=0 ; (i<nCliques) ; ++i )
   {
      const int nEl = clq_set_clique_size( cgraph->clqSet, i );
      const int *el = clq_set_clique_elements( cgraph->clqSet, i );

      int newClqSize = 0;
      int j;
      for ( j=0 ; (j<nEl) ; ++j )
      {
         if ( nindexes[ el[j] ] != -1 )
            neighs[ newClqSize++ ] = nindexes[ el[j] ];
      }
      if ( newClqSize > 1 )
      {
         if ( newClqSize < 50 )
         {
            /* not adding as clique anymore */
            cgraph_add_clique_as_normal_conflicts( result, newClqSize, neighs );

            continue;
         }
         else
            cgraph_add_clique( result, neighs, newClqSize );
      }
   }

   free( neighs );

   cgraph_update_min_max_degree( result );

   return result;
}

int cgraph_get_node_weight( const CGraph *cgraph, int node )
{
   if (cgraph->w)
      return cgraph->w[node];

   return 0;
}

void cgraph_set_node_weight( CGraph *cgraph, int node, int weight )
{
   if (!cgraph->w)
      cgraph->w = xmalloc( sizeof(int)*cgraph_size(cgraph));

   cgraph->w[node] = weight;
}

const char *cgraph_get_node_name( const CGraph *cgraph, int node )
{
   if (cgraph->nodeNames)
      return cgraph->nodeNames[node];

   return NULL;
}

void cgraph_set_node_name( CGraph *cgraph, int node, const char *name )
{
   if (!cgraph->nodeNames)
   {
      cgraph->nodeNames = xmalloc( sizeof(char*)*cgraph_size(cgraph) );
      cgraph->nodeNames[0] = xmalloc( sizeof(char)*MAX_NAME_SIZE*cgraph_size(cgraph) );
      int i;
      for ( i=1 ; (i<cgraph_size(cgraph)) ; ++i )
         cgraph->nodeNames[i] = cgraph->nodeNames[i-1] + MAX_NAME_SIZE;
   }

   strncpy( cgraph->nodeNames[node], name, MAX_NAME_SIZE );
}

const int *cgraph_get_node_weights( const CGraph *cgraph )
{
   return cgraph->w;
}

int cgraph_weight( const double w )
{
   return (int)(w*1000.0);
}

int cgraph_get_original_node_index( const CGraph *cgraph, const int node )
{
   if (cgraph->origIdx)
      return cgraph->origIdx[node];

   return -1;
}


#ifdef DEBUG
void cgraph_check_node_cliques( const CGraph *cgraph )
{
   /* all nodes */
   int i,j;
   for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
   {
      const int nCliques = vint_set_size( cgraph->nodeCliques + i );
      const int *el = vint_set_get_elements( cgraph->nodeCliques + i );
      for ( j=0 ; ( j<nCliques ) ; ++j )
      {
         const int idx = el[j];
         if (!(clq_set_clique_has_element( cgraph->clqSet, idx, i )))
         {
            printf("\nnode %d should appear on clique %d but does not\n", i, idx );
            fflush(stdout);
         }
         assert( clq_set_clique_has_element( cgraph->clqSet, idx, i ) );
      }
   }

   /* printf("Information in node cliques indicate cliques which really have the node.\n"); */

   const int nc = clq_set_number_of_cliques( cgraph->clqSet );
   for ( i=0 ; (i<nc) ; ++i )
   {
      const int n   = clq_set_clique_size( cgraph->clqSet, i );
      const int *el = clq_set_clique_elements( cgraph->clqSet, i );

      for ( j=0 ; (j<n) ; ++j )
      {
         const int currNode = el[j];
         /* this must appear in node clique */
         int l;
         int nn = vint_set_size( cgraph->nodeCliques+currNode );
         const int *eel = vint_set_get_elements( cgraph->nodeCliques+currNode );

         for ( l=0 ; (l<nn) ; ++l )
            if (eel[l] == i)
               break;

         if ( l == nn )
         {
            fprintf( stderr, "ERROR: in clique %d node %d appears but this clique does not appears in its node list.\n", i, currNode );
            fflush( stdout );
            fflush( stderr );
         }

         assert( l<nn );
      }
   }
}

void cgraph_check_neighs( const CGraph *cgraph )
{
   const int n = cgraph_size( cgraph );
   int nNeighs;
   int *neighs = xmalloc( sizeof(int)*n*100 );

   int i, j;
   for ( i=0 ; (i<n) ; ++i )
   {
      nNeighs = cgraph_get_all_conflicting( cgraph, i, neighs, n*100 );

      /* computed number of neighs */
      int cn = 0;
      for ( j=0 ; ( j<n ) ; ++j )
         if (cgraph_conflicting_nodes( cgraph, i, j ))
            cn++;

      assert( cn == nNeighs );

      for ( j=0 ; ( j<n ) ; ++j )
         if (cgraph_conflicting_nodes( cgraph, i, j ))
            assert ( bsearch( &j, neighs, nNeighs, sizeof(int), vint_set_cmp_int ) );
         else
            assert ( !bsearch( &j, neighs, nNeighs, sizeof(int), vint_set_cmp_int ) );
   }

   free( neighs );
}

void cgraph_check_preproc( const CGraph *ppgraph, const CGraph *cgraph )
{
   assert( ppgraph->origIdx );
   assert( cgraph_size(cgraph) >= cgraph_size(ppgraph) );
   const int n = cgraph_size( ppgraph );
   int i,j;

   int *ppIdx = xmalloc( sizeof(int)*cgraph_size(cgraph) );
   for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
      ppIdx[i] = -1;

   for ( i=0 ; (i<cgraph_size(ppgraph)) ; ++i )
      ppIdx[ ppgraph->origIdx[i] ] = i;

   for ( i=0 ; (i<n) ; ++i )
   {
      for ( j=i+1 ; (i<n) ; ++i )
      {
         const int cij = cgraph_conflicting_nodes( ppgraph, i, j );
         const int cji = cgraph_conflicting_nodes( ppgraph, j, i );
         assert( cij == cji );

         const int oi = ppgraph->origIdx[i];
         assert( oi != -1 );
         const int oj = ppgraph->origIdx[j];
         assert( oj != -1 );

         const int cijo = cgraph_conflicting_nodes( cgraph, oi, oj );
         const int cjio = cgraph_conflicting_nodes( cgraph, oj, oi );

         assert( cijo == cjio );

         if ( cijo != cij )
         {
            fprintf( stderr, "Error: Conflicts different in original and ppgraph.\n" );
         }
         assert( cijo == cij );
         assert( cjio == cji );
      }
   }

   for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
   {
      if ( ppIdx[i] == -1 )
         continue;
      const int pIdxI = ppIdx[i];
      for ( j=i+1 ; (j<cgraph_size(cgraph)) ; ++j )
      {
         if ( ppIdx[j] == -1 )
            continue;
         const int pIdxJ = ppIdx[j];
         const int co = cgraph_conflicting_nodes( cgraph, i, j );
         const int ci = cgraph_conflicting_nodes( ppgraph, pIdxI, pIdxJ );
         if (co!=ci)
            printf("bra");
         assert( co==ci );
      }
   }


   free( ppIdx );
}

#endif

int cmp_by_node_degree( const void *v1, const void *v2 )
{
   const NodeDegree *nd1 = v1;
   const NodeDegree *nd2 = v2;

   if ( nd2->degree != nd1->degree )
      return (nd1->degree - nd2->degree);

   return (nd1->node - nd2->node);
}

const int *cgraph_get_original_node_indexes( const CGraph *cgraph )
{
   return cgraph->origIdx;
}

void cgraph_print_summary( CGraph *cgraph, const char *name )
{
   cgraph_update_min_max_degree( cgraph );
   int nodes = cgraph_size( cgraph );
   int nSmallDegree = 0;
   int edges = 0;
   {
      int i;
      for ( i=0 ; (i<nodes) ; ++i )
      {
         edges += cgraph_degree( cgraph,  i );
         if ( cgraph_degree( cgraph, i ) <= cgraph->lowDegree )
            nSmallDegree++;
      }
   }
   printf("%s : n %d e %d - minD %d maxD %d smallD %d\n", name, cgraph_size(cgraph), edges, cgraph_min_degree(cgraph), cgraph_max_degree(cgraph), nSmallDegree );
}

void cgraph_set_low_degree( CGraph *cgraph, const int lowDegree )
{
   cgraph->lowDegree = lowDegree;
}

void cgraph_add_clique_as_normal_conflicts( CGraph *cgraph, const int nEl, const int elClique[] )
{
   /* using iv as temporary area */
   int i,j,count;
   const int nconflicts=nEl-1;
   int *iv = cgraph->iv;

   for ( i=0 ; (i<nEl) ; ++i )
   {
      const int node = elClique[i];
      count = 0;
      for ( j=0 ; (j<nEl) ; ++j )
      {
         if ( elClique[j] == node )
            continue;

         iv[count++] = elClique[j];
      }
      cgraph_add_node_conflicts_no_sim( cgraph, node, iv, nconflicts );
   }
}

int cgraph_get_n_conflicting( const CGraph *cgraph, int node, int neighs[], int n,
      int v[], const int vcap )
{
   return n;
}


/* NeighIterator functions */

NeighIterator *nit_create( )
{
   NeighIterator *result = xmalloc( sizeof(NeighIterator) );

   result->bhv    = NULL;
   result->bhvCap = 0;

   result->previousNode = -1;

   return result;
}

void nit_check_capacity( NeighIterator *nit, int required )
{
   required++;

   if ( nit->bhvCap < required )
   {
      if ( nit->bhvCap )
      {
         nit->bhvCap = MAX( required, nit->bhvCap*2 );
         nit->bhv = xrealloc( nit->bhv, sizeof(NodeCost)*nit->bhvCap );
      }
      else
      {
         nit->bhvCap = MAX( required, 1024 );
         nit->bhv = xmalloc( sizeof(NodeCost)*nit->bhvCap );
      }
   }
}

/* Heap macros */
/* position of root node in vector */
#define rootPosBHV( node )  ( (((node+1)/2)-1) )
/* position of the first child node in vector */
#define childPosBHV( node ) ( (node*2)+1 )

#define bhv_swap( bhv,  p1, p2 ) { NodeCost tempBHV=bhv[p1]; bhv[p1]=bhv[p2]; bhv[p2]=tempBHV; }

void nit_bhv_up( NodeCost *bhv, int index )
{
   int root, child = index;

   while ( (root=rootPosBHV(child)) >=0 )
   {
      if ( bhv[root].cost > bhv[child].cost )
      {
         bhv_swap( bhv, child, root );
         child = root;
      }
      else
         return;
   }
}

void nit_bhv_down( NodeCost *bhv, int size, int index )
{
   int root = index, child;

#ifdef DEBUG
   assert( index >= 0 );
#endif

   while ( (child=childPosBHV(root)) < size )
   {
#ifdef DEBUG
      assert( child >= 0 );
      assert( child < size );
#endif
      /* child with the smallest cost */
      if ( (child+1<size) && (bhv[child].cost>bhv[child+1].cost) )
         child++;

#ifdef DEBUG
      assert( child >= 0 );
      assert( child < size );
#endif
      if ( bhv[root].cost > bhv[child].cost )
      {
         bhv_swap( bhv, root, child );
         root = child;
      }
      else
         break;
   }
}

void nit_bhv_add( NeighIterator *nit, int node, int cost )
{
#ifdef DEBUG
   assert( node >= 0 );
   assert( (nit->bhvSize>=0) && (nit->bhvSize<nit->bhvCap) );
#endif


   nit->bhv[nit->bhvSize].node = node;
   nit->bhv[nit->bhvSize].cost = cost;

   nit_bhv_up( nit->bhv, nit->bhvSize );

   nit->bhvSize++;
}

void nit_fill_bhv( NeighIterator *nit, const CGraph *cgraph, int node, const int costs[] )
{
#ifdef DEBUG
   assert( node >= 0 );
   assert( node < cgraph_size(cgraph) );
#endif
   {
      const IntSet *nodeConflicts = cgraph->nodeConflicts + node;
      const int nConf = vint_set_size( nodeConflicts );
      const int *conf = vint_set_get_elements( nodeConflicts );
      int i;
      for ( i=0 ; (i<nConf) ; ++i )
      {
#ifdef DEBUG
         assert( (conf[i]>=0) );
         assert( (conf[i]<cgraph_size(cgraph)) );
#endif
         nit_bhv_add( nit, conf[i], costs[conf[i]] );
      }
   }

   {
      /* conflicts stored in cliques */
      const int *cliques = vint_set_get_elements(cgraph->nodeCliques+node);
      const int nCliques = vint_set_size(cgraph->nodeCliques+node);
      const CliqueSet *clqSet = cgraph->clqSet;
      int i;
      for ( i=0 ; (i<nCliques) ; ++i )
      {
         const IntSet *clique = clq_set_get_clique( clqSet, cliques[i] );
         const int size = vint_set_size(clique);
         const int *el = vint_set_get_elements(clique);
         int j;
         for ( j=0 ; (j<size) ; ++j )
            nit_bhv_add( nit, el[j], costs[el[j]] );
      }
   }
}

void nit_start( NeighIterator *nit, const CGraph *cgraph, int node, const int costs[] )
{
   /* nodespace : needed for store neighs and all cliques */
   int nodeSpace      = vint_set_size(cgraph->nodeConflicts+node);
   int nCliques       = vint_set_size(cgraph->nodeCliques+node);
   const int *cliques = vint_set_get_elements(cgraph->nodeCliques+node);
   const CliqueSet *clqSet = cgraph->clqSet;
   {
      int i;
      for ( i=0 ; (i<nCliques) ; ++i )
         nodeSpace += vint_set_size(clq_set_get_clique(clqSet, cliques[i]));
   }

   nit_check_capacity( nit, nodeSpace );

   nit->bhvSize      =  0;
   nit->previousNode = -1;

   nit_fill_bhv( nit, cgraph, node, costs );

   nit->node = node;
}

int nit_next( NeighIterator *nit )
{
   int result = -1;

NEXT:

   if (nit->bhvSize==0)
      return INT_MAX;

   NodeCost *bhv = nit->bhv;

   result = bhv[0].node;

   --nit->bhvSize;
   bhv[0] = bhv[nit->bhvSize];

   nit_bhv_down( bhv, nit->bhvSize, 0 );

   if ((result==nit->previousNode)||(result==nit->node))
      goto NEXT;

   nit->previousNode = result;
   return result;
}

void cgraph_recompute_degree( CGraph *cgraph )
{
   const size_t size = cgraph->nodeSize;
   memset( cgraph->degree, 0, sizeof(unsigned long int)*size );
   int i, j;
   for ( i=0 ; (i<size) ; ++i )
      for ( j=i+1 ; (j<size) ; ++j )
      {
         if (cgraph_conflicting_nodes( cgraph, i, j ))
         {
            cgraph->degree[i]++;
            cgraph->degree[j]++;
         }
      }

   cgraph_update_min_max_degree( cgraph );
}

void nit_free( NeighIterator **nit )
{
   free( (*nit)->bhv );
   free( (*nit) );
   (*nit) = NULL;
}

