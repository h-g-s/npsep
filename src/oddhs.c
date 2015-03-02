#include <float.h>
#include <limits.h>
#include <alloca.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include "spaths.h"
#include "oddhs.h"
#include "memory.h"
#include "macros.h"
#include "node_heap.h"

#define EPS 1e-5             // minimum to a variable to be considered active
#define MIN_FRAC 0.001

#define DIST_INF 33554432    // not INT_MAX to not cause overflow (2^25)

#define DIST_MAX 16384       // maximum  distance for nodes with conflicts
// (< infinity) in the distance graph

#define MAX_WHEEL_CENTERS 256

#include "vectormgm.h"

struct _OddHoleSep
{
   const CGraph *cgraph;

   // from integer interesting columns
   // those which are active in solution
   int icaCount;
   int icaCapacity;
   int *icaIdx;      // original index
   int *icaActivity; // mapping of the fractional solution value to an
   // integer value to made further computations easier

   // shortest path data
   int spNodes;
   int spNodeCap;
   int spArcs;
   int spArcsCap;
   int *spArcStart; // start index for arcs of each node
   int *spArcTo;    // destination of each arc
   int *spArcDist;  // distance for each arc

   // all these vector have a size compatible with ICA
   int *spPathIdx;      // to store discovered shortest path
   int *spOrigPathIdx;      // to store discovered shortest path

   // incidence vector to check for repeated entries
   int *ivreIdx;
   int ivreCap;

   ShortestPathsFinder *spf;

   // discovered odd holes
   // indexes are stored related to icaIdx
   // (only considering active variables)
   int dohCapOH;       // capacity for odd holes
   int dohCapIdx;      // capacity for odd hole elements
   int dohCapWCIdx;    // capacity for odd hole wheel centers
   int dohCount;       // number of stored odd holes
   int *dohIdx;        // indexes of all odd holes
   int *dohWCIdx;      // indexes of all wheel centers
   int *dohStart;      // start indexes for the i-th odd hole
   int *dohWCStart;    // wheel center starts

   char *dohIV;         // incidence vector for variables considered here
   int dohIVCap;
   // icaSize related

   NodeHeap *nheap;

   CliqueSet *clqSet;
};

double oddhs_lhs( const int size, const int idx[], const double x[] )
{
   double r = 0.0;
   int i;
   for ( i=0; (i<size) ; i++ )
      r += x[idx[i]];
   return r;
}

double oddhs_rhs( const int size )
{
   return ((double)(size/2));
}

double oddhs_viol( const int size, const int idx[], const double x[] )
{
   return oddhs_lhs( size, idx, x) - oddhs_rhs( size );
}

double oddhs_frac_part( const double x )
{
   return MIN( fabs(x-floor(x)), fabs(x-ceil(x)) );
}

OddHoleSep *oddhs_create( )
{
   OddHoleSepPtr oddhs = (OddHoleSepPtr) xmalloc( sizeof(OddHoleSep) );

   oddhs->icaCount     = 0;
   oddhs->icaCapacity  = 1024;
   ALLOCATE_INT_VECTOR( oddhs->icaIdx, oddhs->icaCapacity );
   ALLOCATE_INT_VECTOR( oddhs->icaActivity, oddhs->icaCapacity );
   ALLOCATE_INT_VECTOR( oddhs->spPathIdx, oddhs->icaCapacity );
   ALLOCATE_INT_VECTOR( oddhs->spOrigPathIdx, oddhs->icaCapacity );

   oddhs->dohIVCap = 32768;
   oddhs->dohIV = xmalloc( sizeof(char)*oddhs->dohIVCap );

   oddhs->ivreCap    = 2048;
   ALLOCATE_INT_VECTOR( oddhs->ivreIdx, oddhs->ivreCap );

   oddhs->spf          = spf_create();

   oddhs->spNodes    = 0;
   oddhs->spArcs     = 0;
   oddhs->spNodeCap  = 1024;
   oddhs->spArcsCap  = 2048;
   ALLOCATE_INT_VECTOR( oddhs->spArcStart, (oddhs->spNodeCap+1) );
   ALLOCATE_INT_VECTOR( oddhs->spArcTo, oddhs->spArcsCap );
   ALLOCATE_INT_VECTOR( oddhs->spArcDist, oddhs->spArcsCap );

   oddhs->dohCount    = 0;
   oddhs->dohCapOH    = 1024;
   oddhs->dohCapIdx   = 4096;
   ALLOCATE_INT_VECTOR( oddhs->dohIdx, oddhs->dohCapIdx );

   ALLOCATE_INT_VECTOR( oddhs->dohStart, (oddhs->dohCapOH+1) );
   ALLOCATE_INT_VECTOR( oddhs->dohWCStart, (oddhs->dohCapOH+1) );
   oddhs->dohStart[0] = 0;
   oddhs->dohWCStart[0] = 0;

   oddhs->dohCapWCIdx = 8192;
   ALLOCATE_INT_VECTOR( oddhs->dohWCIdx, oddhs->dohCapWCIdx );

   oddhs->nheap = NULL;

   oddhs->clqSet = clq_set_create();

   oddhs->cgraph = NULL;

   return oddhs;
}

// if necessary expands ica to hold a new element
void oddhs_check_ica_capacity( OddHoleSepPtr oddhs ) __attribute__ ((visibility("hidden")));
// methods to work with current fractional solution
int oddhs_fill_active_intcols( OddHoleSepPtr oddhs, const int cols, const double x[], const CGraph * conflicts ) __attribute__ ((visibility("hidden")));
void oddhs_prepare_dist_graph( OddHoleSepPtr oddhs, const double x[], const CGraph * conflicts ) __attribute__ ((visibility("hidden")));

/**
 * Odd Hole store management functions.
 * Odd-holes are processed using ica indexes.
 * Only after all processing that they are
 * translated to original variable indexes.
 **/

/** Clears odd holes set. **/
void oddhs_clear_doh( OddHoleSep *oddhs ) __attribute__ ((visibility("hidden")));;

/** Tries to add a newly discovered odd hole (not necessarily violated).
 *  If it is a repeated entry, ignores it and returns zero, if suceed in
 *  inserting, returns one.
 **/
int oddhs_add_doh( OddHoleSep *oddhs, const int nz, int _idx[], const CGraph *cgraph )  __attribute__ ((visibility("hidden")));

/** tests if a newly discovered odd hole is not a repeated entry **/
int oddhs_doh_already_exists( OddHoleSep *oddhs, const int nz, const int idxs[] )  __attribute__ ((visibility("hidden")));

/** translate all DOHs to original variables indexes **/
/*void oddhs_translate_dohs( OddHoleSep *oddhs );*/

/**
 * Tries to find an Odd Hole using the conflict graph
 * this odd hole may correspond to a violated cut.
 * the search may be aggressive or not
 * returns how many new odd holes were found
 **/
int oddhs_find_odd_holes_with_node( OddHoleSep *oddhs, const int node, const double x[], const CGraph *cgraph ) __attribute__ ((visibility("hidden")));
int oddhs_try_path( OddHoleSep *oddhs, const int i, const int nz, int _idx[], const CGraph *cgraph ) __attribute__ ((visibility("hidden")));
int oddhs_vector_has_repeated_entries( OddHoleSep *oddhs, const int maxIndex, const int nz, const int idx[] ) __attribute__ ((visibility("hidden")));

void oddhs_make_maximal_cliques( OddHoleSep *oddhs, const CGraph *cgraph, const double x[] ) __attribute__ ((visibility("hidden")));

/**
 * tries to find wheel centers so that Odd Hole inequelities may be lifted
 */
int oddhs_search_wheel_centers_all_dohs( OddHoleSep *oddhs, const CGraph *cgraph, const double x[], const double rc[] ) __attribute__ ((visibility("hidden")));
/**
 * finds at most "maxCenters" wheel centers, for one odd hole
 * priority is:
 *   1. most fractional one
 *   2. those with smallest reduced cost
 **/
int oddhs_find_wheel_centers( const int cols, const double *x, const double rc[],
        const CGraph *cgraph, const int *oh, const int ohSize, int centers[],
        const int maxCenters ) __attribute__ ((visibility("hidden")));


void oddhs_check_space_sp_nodes( OddHoleSep *oddhs, const int nodes ) __attribute__ ((visibility("hidden")));
void oddhs_check_space_sp_arcs( OddHoleSep *oddhs, const int arcs ) __attribute__ ((visibility("hidden")));

// transform a double [0 , 1] into an integer
// ready to use in the distance graph
int oddhs_get_activity( const double val ) __attribute__ ((const,visibility("hidden")));
double oddhs_get_activityD( const long int val ) __attribute__ ((const,visibility("hidden")));

double oddhs_get_x_from_ica( OddHoleSep *oddhs, const long int idx )  __attribute__ ((const,visibility("hidden")));

int oddhs_search_odd_holes( OddHoleSep *oddhs, const int cols, const double x[],
        const double rc[], const CGraph *conflicts )
{
   oddhs->cgraph = conflicts;
   oddhs->nheap = nh_create( cols, DIST_INF );

   int icCount = oddhs_fill_active_intcols( oddhs, cols, x, conflicts );
   if ( icCount == 0 )
      return 0;

   int result = 0;

   oddhs_prepare_dist_graph( oddhs, x, conflicts );

   oddhs_clear_doh( oddhs );

   {
      int i;
      for ( i=0 ; (i<oddhs->icaCount) ; i++ )
         result += oddhs_find_odd_holes_with_node( oddhs, i, x, conflicts );
   }

   /* cliques are already translated */
   //oddhs_translate_dohs( oddhs );

   oddhs_make_maximal_cliques( oddhs, conflicts, x );

   oddhs_search_wheel_centers_all_dohs( oddhs, conflicts, x, rc);

   nh_free( &(oddhs->nheap) );

   return result;
}

int oddhs_search_wheel_centers_all_dohs( OddHoleSep *oddhs, const CGraph *cgraph, const double x[], const double rc[] )
{
   const int cols = cgraph_size( cgraph );
   /* trying to extend odd holes by including wheel centers */
   int *wheelc = alloca( MIN(MAX_WHEEL_CENTERS,cols)*sizeof(int) );
   int i, insertedWC = 0;
   const int nOH = oddhs->dohCount;
   for ( i=0 ; (i<nOH) ; ++i )
   {
      oddhs->dohWCStart[i] = insertedWC;

      const int *ohs = oddhs_get_odd_hole( oddhs, i );
      const int sizeOH = oddhs_get_odd_hole( oddhs, i+1 ) - ohs;
      const int nwc =
         oddhs_find_wheel_centers( cols, x, rc, cgraph, ohs, sizeOH, wheelc, MAX_WHEEL_CENTERS );

      if (!nwc)
         continue;

      /* adding more indexes */
      if (insertedWC+nwc+1 > oddhs->dohCapWCIdx)
      {
         /* realocating indexes and updating
            * pointers in dohWCStart */
         oddhs->dohCapWCIdx = MAX( oddhs->dohCapWCIdx*2, insertedWC+nwc+1 );
         oddhs->dohWCIdx = xrealloc( oddhs->dohWCIdx, sizeof(int)*oddhs->dohCapWCIdx );
      }

      memcpy( oddhs->dohWCIdx + oddhs->dohWCStart[i], wheelc, sizeof(int)*nwc );

      insertedWC += nwc;
   }
   oddhs->dohWCStart[i] = insertedWC;

   return insertedWC;
}

void oddhs_check_ica_capacity( OddHoleSepPtr oddhs )
{
   if (oddhs->icaCount+1>=oddhs->icaCapacity)
   {
      oddhs->icaCapacity *= 2;
      oddhs->icaIdx = xrealloc( oddhs->icaIdx, sizeof(int)*oddhs->icaCapacity );
      oddhs->icaActivity = xrealloc( oddhs->icaActivity, sizeof(int)*oddhs->icaCapacity );
      oddhs->spPathIdx = xrealloc( oddhs->spPathIdx, sizeof(int)*oddhs->icaCapacity );
      oddhs->spOrigPathIdx = xrealloc( oddhs->spOrigPathIdx, sizeof(int)*oddhs->icaCapacity );
   }
}

int oddhs_fill_active_intcols( OddHoleSepPtr oddhs, const int cols, const double x[], const CGraph *conflicts )
{
   // conflict nodes
   int cnCount = cgraph_size( conflicts );

   if ( cnCount<3 )
      return 0;
   //   const int *cn     = cgraph_size( conflicts );
   //int posIca = 0;
   oddhs->icaCount  = 0;
   {
      int j;
      for ( j=0 ; (j<cnCount)  ; j++ )
      {
         oddhs_check_ica_capacity( oddhs );
         if ((x[j]<=EPS)||(oddhs_get_activity(x[j])==0))
            continue;
         if ( cgraph_degree( conflicts, j ) < 3 )
            continue;

         //printf("j %d   x %g   activity %d posIca %d\n", j, x[j], oddhs_get_activity( x[j] ), posIca++ );
         oddhs->icaIdx[oddhs->icaCount]       = j;
         oddhs->icaActivity[oddhs->icaCount]  = oddhs_get_activity( x[j] );
         oddhs->icaCount++;
      }
   }

   return oddhs->icaCount;
}

void oddhs_prepare_dist_graph( OddHoleSepPtr oddhs, const double x[], const CGraph * conflicts )
{
   const int nodes = oddhs->icaCount*2;
   oddhs_check_space_sp_nodes( oddhs, nodes+1 );

   //printf("icaCount: %d\n", oddhs->icaCount );

   int idxArc = 0;

   // process arcs from G to G' and from G' to G
   {
      int  i;
      for ( i=0 ; i<2 ; ++i )
      {
         int shift1=0, shift2=0;
         switch (i)
         {
            case 0:
               shift1 = 0;
               shift2 = oddhs->icaCount;
               break;
            case 1:
               shift1 = oddhs->icaCount;
               shift2 = 0;
               break;
         }

         {
            int i1;
            for ( i1=0 ; (i1<oddhs->icaCount) ; ++i1 )
            {
               //printf("AAA: %d --- %d === %d\n",shift1+i1, oddhs->icaCount, oddhs->spNodeCap );
               oddhs->spArcStart[ shift1+i1 ] = idxArc;
               const int idx1 = oddhs->icaIdx[i1];
               {
                  int i2;
                  for ( i2=0 ; (i2<oddhs->icaCount) ; ++i2 )
                  {
                     const int idx2 = oddhs->icaIdx[i2];
                     if ( cgraph_conflicting_nodes( conflicts, idx1, idx2 ) )
                     {
                        //printf("arc  %d_%d (%d-%d) - %d\n", i1, i2, shift1+i1, shift2+i2, oddhs->icaActivity[i2] );
                        oddhs_check_space_sp_arcs( oddhs, idxArc );
                        oddhs->spArcTo[idxArc]   = shift2+i2;
                        oddhs->spArcDist[idxArc] = oddhs->icaActivity[i2];
                        ++idxArc;
                     } // conflict found
                  } // i2
               }
            } // i1
         }
      }  // for G and G'
   }
   oddhs->spArcStart[nodes] = idxArc;

   //printf("n: %d a: %d\n", nodes, idxArc ); fflush(stdout);
   spf_update_graph( oddhs->spf, nodes, idxArc, oddhs->spArcStart, oddhs->spArcTo, oddhs->spArcDist );
}

int oddhs_get_activity( const double val )
{
   int intval = (int) ((((double)DIST_MAX)*val) + 0.5);

   if ( intval>DIST_MAX )
      intval = DIST_MAX;
   else
      if ( intval<0 )
         intval = 0;

   return DIST_MAX-intval;
}

double oddhs_get_activityD( const long int val )
{
   const double dmax = DIST_MAX;
   const double dval = dmax-val;
   double result = dval / dmax;
   return result;
}

int oddhs_try_path( OddHoleSep *oddhs, const int i, const int nz, int _idx[], const CGraph *cgraph )
{
   if ( nz <=2  )
      return 0;

      // translating indexes
   {
      int j;
      for ( j=0 ; j<nz ; ++j )
         _idx[j] %= oddhs->icaCount;
   }

   // checking for repeated entries
   if ( oddhs_vector_has_repeated_entries( oddhs, oddhs->icaCount, nz-1, _idx ) )
      return 0;

   return oddhs_add_doh( oddhs, nz-1, _idx, cgraph );
}

int oddhs_find_odd_holes_with_node( OddHoleSep *oddhs, const int node, const double x[], const CGraph *cgraph )
{
   //printf("node %d.\n", node );

   const int icaCount = oddhs->icaCount;
   const int dest = icaCount + node;
   //const int distBlocked = DIST_MAX - 1;
   ShortestPathsFinder *spf = oddhs->spf;
   int origNZ;
   int *origPathIdx = oddhs->spOrigPathIdx;
   int result = 0;

   spf_find( spf, node );
   origNZ = spf_get_path( spf, dest, origPathIdx );
   if (!origNZ)
      return 0;

   result += oddhs_try_path( oddhs, node, origNZ, origPathIdx, cgraph );

   if (!result)
      return 0;

   /* checking if it is violated */
   {
      int i;
      double lhs = 0;
      for ( i=0; (i<origNZ-1) ; i++ )
      {
/*
         printf("1 %d\n", i);
         printf("2 %d %d\n", origPathIdx[i], oddhs->icaCount );
         printf("3 %d\n", oddhs->icaIdx[origPathIdx[i]]);
         printf("4 %d\n", cgraph_size(cgraph) );
*/
         assert( origPathIdx[i]<oddhs->icaCount );
         lhs += x[oddhs->icaIdx[origPathIdx[i]]];
      }
      const double viol = lhs - ((double)( (origNZ-1) / 2 ));
      if ( viol <= EPS )
         return result;
   }

   return result;
   /*  if it is, modifying every arc in the path and running again to try to find other ones */
   /* THIS IDEA DOES NOT WORKED, DEAD CODE FOR NOW */
   /*{
      int *otherPath = alloca( sizeof(int)*(oddhs->icaCount+2) );
      int i;
      for ( i=0 ; (i<origNZ-1) ; i++)
      {
         const int head = origPathIdx[i];
         const int tail = origPathIdx[i+1];
         const int origDist = spf_get_arc( spf, head, oddhs->icaCount+tail );

         spf_update_arc( spf, head, oddhs->icaCount+tail, distBlocked );
         spf_update_arc( spf, oddhs->icaCount+tail, head, distBlocked );

         int newNZ = spf_get_path( spf, dest, otherPath );
         if (newNZ)
         {
            const int newOH = oddhs_try_path( oddhs, node, newNZ, otherPath );
            result += newOH;
            if (newOH)
            {
               printf("====== Modified odd hole found!.\n");
            }
         }

         spf_update_arc( spf, head, oddhs->icaCount+tail, origDist );
         spf_update_arc( spf, oddhs->icaCount+tail, head, origDist );
      }
   }*/

}

int oddhs_vector_has_repeated_entries( OddHoleSep *oddhs, const int maxIndex, const int nz, const int idx[] )
{
   if (nz<=1)
      return 0;

   const int ivSpace = maxIndex + 1;
   ADJUST_INT_VECTOR_CAPACITY( oddhs->ivreIdx, oddhs->ivreCap, ivSpace );
   memset( oddhs->ivreIdx, 0, sizeof(int)*ivSpace );

   register int i;
   for ( i=0 ; (i<nz) ; i++ )
   {
//      printf(">>>%d %d\n", idx[i], oddhs->icaCount );
      assert( idx[i]<ivSpace );
      if ( oddhs->ivreIdx[ idx[i] ]++ )
         return 1;
   }

   return 0;
}

void oddhs_check_space_sp_nodes( OddHoleSep *oddhs, const int nodes )
{
   if ( nodes+1 >= oddhs->spNodeCap )
   {
      oddhs->spNodeCap = MAX(2*oddhs->spNodeCap, nodes );
      oddhs->spArcStart = xrealloc( oddhs->spArcStart, oddhs->spNodeCap*sizeof(int) );
      return;
   }
}

void oddhs_check_space_sp_arcs( OddHoleSep *oddhs, const int arcs )
{
   if ( arcs+1 >= oddhs->spArcsCap )
   {
      oddhs->spArcsCap *= 2;
      oddhs->spArcTo    = xrealloc( oddhs->spArcTo, sizeof(int)*oddhs->spArcsCap  );
      oddhs->spArcDist  = xrealloc( oddhs->spArcDist, sizeof(int)*oddhs->spArcsCap  );
      return;
   }
}

double oddhs_get_x_from_ica( OddHoleSep *oddhs, const long int idx )
{
   int i;
   for ( i=0 ; (i<oddhs->icaCount) ; i++ )
      if ( oddhs->icaIdx[i] == idx )
         return oddhs_get_activityD( oddhs->icaActivity[i] );

   fprintf( stderr, "Error in oddhs_get_x_from_ica: idx %d not found.\n", (int)idx );
   exit(1);
}

void oddhs_clear_doh( OddHoleSep *oddhs )
{
   oddhs->dohCount    = 0;
   oddhs->dohStart[0] = 0;
}

int oddhs_add_doh( OddHoleSep *oddhs, const int nz, int _idx[], const CGraph *cgraph )
{

   int *idx = alloca( sizeof(int)*nz );
   {
      int i;
      for ( i=0 ; (i<nz) ; ++i )
         idx[i] = oddhs->icaIdx[_idx[i]];
   }

   // checking for repeated entries
   if (oddhs_doh_already_exists( oddhs, nz, idx ))
      return 0;


   if ( nz == 3 )
   {
      int w = 0;
      {
         int i;
         for ( i=0; (i<nz) ; i++ ) /* computing weight */
            w += cgraph_get_node_weight( cgraph, idx[i] );
      }
      qsort_int( idx, nz );
      return clq_set_add( oddhs->clqSet, nz, idx, w );
   }

   // checking space
   if ( oddhs->dohCount+2 > oddhs->dohCapOH )
   {
      oddhs->dohCapOH *= 2;

      oddhs->dohStart = xrealloc( oddhs->dohStart, sizeof(int)*oddhs->dohCapOH );
      oddhs->dohWCStart = xrealloc( oddhs->dohWCStart, sizeof(int)*oddhs->dohCapOH );
   }

   const int reqSizeIdx = oddhs->dohStart[oddhs->dohCount] + nz + 1;
   ADJUST_INT_VECTOR_CAPACITY( oddhs->dohIdx, oddhs->dohCapIdx , ((int)(reqSizeIdx)) );

   // inserting
   const int start = oddhs->dohStart[oddhs->dohCount];
   memcpy( oddhs->dohIdx+start, idx, (sizeof(int)*nz) );

   oddhs->dohCount++;
   oddhs->dohStart[ oddhs->dohCount ] = start + nz;

   return 1;
}

int oddhs_doh_already_exists ( OddHoleSep *oddhs, const int nz, const int idxs[] )
{
   const CGraph *cgraph = oddhs->cgraph;
   const int cols = cgraph_size( cgraph );
   if ( cols+1 > oddhs->dohIVCap )
   {
      oddhs->dohIVCap = MAX( oddhs->dohIVCap*2, cols+1 );
      oddhs->dohIV = xrealloc( oddhs->dohIV, sizeof(char)*oddhs->dohIVCap );
   }

   memset( oddhs->dohIV, 0, sizeof(char)*(cols+1) );

   {
      int i;
      for ( i=0 ; (i<nz) ; ++i )
         oddhs->dohIV[ idxs[i] ]++;
   }

   {
      int idxDOH;
      for ( idxDOH=0 ; (idxDOH<oddhs->dohCount) ; ++idxDOH )
      {
         // checking size
         const int otherSize = oddhs->dohStart[idxDOH+1]-oddhs->dohStart[idxDOH];
         if (nz!=otherSize)
            continue;

         const int *dohIdx = oddhs->dohIdx + oddhs->dohStart[idxDOH];
         int j = 0;
         for (  ; (j<nz) ; ++j )
         {
            if (!(oddhs->dohIV[dohIdx[j]]))
               break;
         }
         if ( j == nz )
            return 1;
      }
   }

   return 0;
}

int *oddhs_get_odd_hole( OddHoleSep *oddhs, const int idx )
{
   return oddhs->dohIdx + oddhs->dohStart[idx];
}

int oddhs_get_odd_hole_count( OddHoleSep *oddhs )
{
   return oddhs->dohCount;
}

/*
void oddhs_translate_dohs( OddHoleSep *oddhs )
{
   int *p = oddhs->dohIdx;
   const int *pend = oddhs->dohIdx + oddhs->dohStart[ oddhs->dohCount ];
   for ( ; (p<pend) ; ++p )
      *p = oddhs->icaIdx[ *p ];
}
*/

int oddhs_find_wheel_centers( const int cols, const double *x, const double rc[],
        const CGraph *cgraph, const int *oh, const int ohSize, int centers[],
        const int maxCenters )
{
   int *soh = alloca( sizeof(int)*ohSize );
   memcpy( soh, oh, sizeof(int)*ohSize );
   qsort_int( soh, ohSize );

   int nCenters = 0;
   NodeHeap *nh = nh_create( cols, INT_MAX );

   /* filling priority queue */
   int i;
   for ( i=0 ; (i<cgraph_size(cgraph)) ; ++i )
   {
      if ( cgraph_degree( cgraph, i ) < ohSize )
         continue;

      /* checking if it conflicts with all nodes and is not included in oh */
      if ( bsearch_int( soh, ohSize, i ) )
         continue;

      int j;
      for ( j=0 ; (j<ohSize) ; ++j )
      {
         const int el = oh[j];

         if (!cgraph_conflicting_nodes( cgraph, i, el ))
            break;
      }

      if ( j != ohSize )
         continue;

      int priority = INT_MAX / 2;
      if ( x[i] >= MIN_FRAC )
         priority = (int)(x[i] * 1000.0);
      else /* avoiding negative numbers */
         priority = (int)(1000000.0 + rc[i]);

      nh_update( nh, i, priority );
   }

   int p = INT_MAX, n = INT_MAX;
   while ( ((p=nh_remove_first(nh, &n))!=INT_MAX) &&
           (nCenters<MAX_WHEEL_CENTERS) )
   {
      /* must have conflict with all other centers */
      for ( i=0 ; (i<nCenters) ; ++i )
          if (!cgraph_conflicting_nodes( cgraph, n, centers[i]))
              continue;

      centers[nCenters] = n;
      
      nCenters++;
   }

   nh_free( &nh );

   return nCenters;
}

const CliqueSet *oddhs_get_cliques( OddHoleSep *oddhs )
{
   return oddhs->clqSet;
}

void oddhs_make_maximal_cliques( OddHoleSep *oddhs, const CGraph *cgraph, const double x[] )
{
   NodeHeap *nh = oddhs->nheap;
   CliqueSet *clqSet = oddhs->clqSet;

   CliqueSet *newClqSet = clq_set_create();

   int *newClique = xmalloc( sizeof(int)*oddhs->icaCount );
   int *neighs = xmalloc( sizeof(int)*(cgraph_size(cgraph)*2) ); /* temp area needed (*2) */

   int ic;
   for ( ic=0 ; (ic<clq_set_number_of_cliques(clqSet)) ; ++ic )
   {
      const IntSet *clq = clq_set_get_clique( clqSet, ic );

      const int *el = vint_set_get_elements( clq );
      const int n = vint_set_size( clq );
      assert( (n)&&(el) );

      /* checking item with less neighbors */
      int minNeig = INT_MAX;
      int nodeLN = -1;
      int i;
      for ( i=0 ; (i<n) ; ++i )
      {
         if (cgraph_degree( cgraph, el[i] ) < minNeig )
         {
            minNeig = cgraph_degree( cgraph, el[i] );
            nodeLN = el[i];
         }
      }
      assert( nodeLN!=-1 );

      /* getting these neighbors and adding fractional ones
         (which are not in the clique) to nodeheap */
      nh_reset( nh );
      int nn = cgraph_get_all_conflicting( cgraph, nodeLN, neighs, cgraph_size(cgraph)*2 );
      for ( i=0 ; (i<nn) ; ++i )
      {
         const int idxn = neighs[i];
         const double fp = oddhs_frac_part( x[idxn] );
         if ( (fp<MIN_FRAC) || (vint_set_find( clq, idxn )) || (cgraph_degree(cgraph, idxn)<3) )
            continue;
         const int ifp = 501 - (int)(fp*1000.0);
         nh_update( nh, idxn, ifp );
      } /* filling heap */

      /* new clique */
      memcpy( newClique, el, sizeof(int)*n );
      int newCliqueSize = n;

      /* removing elements from heap and trying to extend clique */
      int fnh, cost;  /* first node heap */
      while ( (cost=nh_remove_first(nh,&fnh))!=DIST_INF )
      {
         /* checking if it is neighbor of all others */
         int j;
         for ( j=0 ; (j<newCliqueSize) ; ++j )
         {
            if (!cgraph_conflicting_nodes(cgraph, newClique[j], fnh ))
               break;
         } /* checking conflicts */
         if (j==newCliqueSize)
            newClique[newCliqueSize++] = fnh;
      } /* trying to extend */
      qsort_int( newClique, newCliqueSize );
      clq_set_add( newClqSet, newCliqueSize, newClique, 0 );
   } /* clique */

   clq_set_clear( clqSet );
   clq_set_add_cliques( clqSet, newClqSet );

   clq_set_free( &newClqSet );
   free( newClique );
   free( neighs );
}

int oddhs_get_nwc_doh( OddHoleSep *oddhs, const int doh )
{
   assert( doh<oddhs->dohCount );
   return oddhs->dohWCStart[doh+1] - oddhs->dohWCStart[doh];
}


const int *oddhs_get_wc_doh( OddHoleSep *oddhs, const int doh )
{
   assert( doh<oddhs->dohCount );
   return oddhs->dohWCIdx + oddhs->dohWCStart[doh];
}


void oddhs_free( OddHoleSepPtr *oddhs )
{
   free( (*oddhs)->icaIdx );
   free( (*oddhs)->icaActivity );
   free( (*oddhs)->spArcStart );
   free( (*oddhs)->spArcTo );
   free( (*oddhs)->spArcDist );
   free( (*oddhs)->dohIdx );
   free( (*oddhs)->dohStart );
   free( (*oddhs)->dohWCIdx );
   free( (*oddhs)->dohWCStart );
   free( (*oddhs)->dohIV );
   free( (*oddhs)->spPathIdx );
   free( (*oddhs)->spOrigPathIdx );
   free( (*oddhs)->ivreIdx );
   clq_set_free( &((*oddhs)->clqSet) );

   spf_free( &((*oddhs)->spf) );

   free((*oddhs));
   *oddhs = NULL;
}
