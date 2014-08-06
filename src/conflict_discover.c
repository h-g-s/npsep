#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include "memory.h"
#include "conflict_discover.h"
#include "vectormgm.h"

#define EPS 1e-6
#define IS_GREATER_THAN( lb, ub ) ( (ub == DBL_MAX) && (lb < DBL_MAX) )
/* exchanges values of 2 variables */
#define SWAP( a,b ) { typeof(a) c=a; a=b; b=c; }
/* mininum size for a row to be considered a clique row */
#define MIN_CLIQUE_ROW 250

#define IS_NEGATIVE_DBL( v ) ( v < -1e-5 )

// conflict vector dimensions
#define CVEC_CAP   1800000
#define CVEC_FLUSH 1700000

typedef struct
{
   int i;
   int j;
} Conflict;

struct _ConflictDiscover
{
   CGraph *cgraph;

   /* lp data, which will be checked */
   int lpCols;
   int lpRows;

   /* Callbacks : */
   LPCallbacks lpc;

   double *colLB;
   double *colUB;
   char   *isInteger;

   /* to pick all coeficients
      of the inequality */
   int *idx;
   double *coef;

   /* indexes of conflicting variables */
   int *confIdx;

   /* conflict vector */
   Conflict *cvec;
   int cvecSize;
   int cvecCap;

   int colCap;

   /* discovered conflicts per row */
   int *rowDC;
   int capRowDC;

   char *sense;
   double *rhs;
};

void cd_add_conflicts_cgraph( ConflictDiscover *cd, Conflict *conflicts, int nConflicts );

int cd_conflict_order_i( const void *e1, const void *e2 )
{
   const Conflict *c1 = (const Conflict *) e1;
   const Conflict *c2 = (const Conflict *) e2;

   if ( c1->i != c2->i )
      return ( c1->i - c2->i );

   return ( c1->j - c2->j );
}

int cd_conflict_order_j( const void *e1, const void *e2 )
{
   const Conflict *c1 = (const Conflict *) e1;
   const Conflict *c2 = (const Conflict *) e2;

   if ( c1->j != c2->j )
      return ( c1->j - c2->j );

   return ( c1->i - c2->i );
}


ConflictDiscover *cd_create(
   CGraph *_cgraph,
   const void *_lp,
   LPCallbacks *lpc
   )
{
   printf("\ncreating conflict graph ... "); fflush( stdout );

   clock_t start = clock();

   ConflictDiscover *cd = xmalloc( sizeof(ConflictDiscover) );

   cd->cgraph = _cgraph;

   cd->lpc = *lpc;

   assert( lpc->cb_get_col_info != NULL );
   assert( lpc->cb_get_row_info != NULL );
   assert( lpc->cb_get_num_cols != NULL );
   assert( lpc->cb_get_num_rows != NULL );
   assert( lpc->cb_get_row_name != NULL );
   assert( lpc->cb_get_col_name != NULL );
   assert( lpc->cb_get_row_sense != NULL );
   assert( lpc->cb_get_rhs != NULL );

   cd->lpCols = cd->lpc.cb_get_num_cols( (void*) _lp );
   cd->lpRows = cd->lpc.cb_get_num_rows( (void*) _lp );

   cd->colLB     = xmalloc( sizeof(double)*cd->lpCols );
   cd->colUB     = xmalloc( sizeof(double)*cd->lpCols );
   cd->isInteger = xmalloc( sizeof(char)*cd->lpCols );
   cd->idx       = xmalloc( sizeof(int)*cd->lpCols );
   cd->coef      = xmalloc( sizeof(double)*cd->lpCols );
   cd->confIdx   = xmalloc( sizeof(int)*cd->lpCols );
   cd->colCap    = cd->lpCols;
   cd->rhs       = xmalloc( sizeof(double)*cd->lpRows );
   cd->sense     = xmalloc( sizeof(char)*cd->lpRows );

   /* Conflict vector */
   cd->cvecCap   = CVEC_CAP;
   cd->cvecSize  = 0;
   cd->cvec      = xmalloc( sizeof(Conflict)*cd->cvecCap );

   cd->capRowDC = cd->lpRows;
   cd->rowDC = xmalloc( sizeof(int)*cd->capRowDC );

   /* first, getting rhs and sense to see candidates for GUB constraints */
   lpc->cb_get_row_sense( (void*)_lp, cd->sense );
   lpc->cb_get_rhs( (void*)_lp, cd->rhs );
   
   /* filling information from initial LP */
   cd_update_col_info( cd, (void*)_lp, 0, cd->lpCols );
   int i;
   cd->cvecSize = 0;
   for ( i=0 ; (i<cd->lpRows) ; ++i )
   {
      if ( ( cd->sense[i] == 'E' || cd->sense[i] == 'L' ) && ( (fabs(cd->rhs[i]-1.0)<2.0001) ) )
      {
         cd_process_row( cd, (void*)_lp, i );
         if ( cd->cvecSize > CVEC_FLUSH )
         {
            cd_add_conflicts_cgraph( cd, cd->cvec, cd->cvecSize );
            cd->cvecSize = 0;
         }
      }
   }

   if (cd->cvecSize)
   {
      cd_add_conflicts_cgraph( cd, cd->cvec, cd->cvecSize );
      cd->cvecSize = 0;
   }
  

   clock_t end = clock();
   
   const double seconds = ((double)(end-start)) / ((double)CLOCKS_PER_SEC);
   printf("done in %.4f seconds.\n", seconds); fflush( stdout );

   return cd;
}

void cd_update_col_info( ConflictDiscover *cd, void *lp, int start, int end )
{
   int i;
   for ( i=start ; (i<end) ; ++i )
   {
      cd->lpc.cb_get_col_info( lp, i, cd->colLB+i, cd->colUB+i, cd->isInteger+i );
      if ( cd->lpc.cb_get_col_name != NULL )
         cgraph_set_node_name( cd->cgraph, i, cd->lpc.cb_get_col_name(lp,i) );
   }
}

void cd_update_graph( ConflictDiscover *cd, const void *lp )
{
   int changed = 0;

   /* checking possible new dimensions */
   int cols = cd->lpc.cb_get_num_cols( (void*)lp );
   int rows = cd->lpc.cb_get_num_rows( (void*)lp );

   if ( cols > cd->lpCols ) /* more information needed */
   {
      if ( cols > cd->colCap )  /* more capacity needed */
      {
         cd->colCap = (int) ( 1.5*((double)cd->colCap) );   /* growing at */
         if ( cd->colCap < cols )                           /* at least 50% */
            cd->colCap = cols;

         cd->colLB = xrealloc( cd->colLB, sizeof(double)*cd->colCap );
         cd->colUB = xrealloc( cd->colUB, sizeof(double)*cd->colCap );
         cd->isInteger = xrealloc( cd->isInteger, sizeof(char)*cd->colCap );
         cd->idx = xrealloc( cd->idx, sizeof(int)*cd->colCap );
         cd->confIdx = xrealloc( cd->confIdx, sizeof(int)*cd->colCap );
         cd->coef = xrealloc( cd->coef, sizeof(double)*cd->colCap );
      }

      cd_update_col_info( cd, (void*)lp, cd->lpCols, cols );
      cd->lpCols = cols;
      changed = 1;
   }

   if ( rows > cd->lpRows )
   {
      int i;
      for ( i=cd->lpRows ; (i<rows) ; ++i )
         cd_process_row( cd, (void*)lp, i );

      cd->lpRows = rows;
      changed = 1;
   }

   if (changed)
      cgraph_update_min_max_degree( cd->cgraph );
}

void fill_row_descr( ConflictDiscover *cd, void *lp, int row, char *descr )
{
   sprintf( descr, "%s;%c;%g", cd->lpc.cb_get_row_name( lp, row ), cd->sense[row], cd->rhs[row] );
}

int cd_process_row( ConflictDiscover *cd, void *lp, int row )
{
   double lb, ub;
   int size;

   cd->lpc.cb_get_row_info ( lp, row, &lb, &ub, &size, cd->idx, cd->coef );

   if (size<2)
      return 0;

   if ( row+1>cd->capRowDC )
   {
      cd->capRowDC = ((double)cd->capRowDC)*1.5;
      if ( row+1>cd->capRowDC )
         cd->capRowDC = row+1;
      cd->rowDC = xrealloc( cd->rowDC, sizeof(int)*cd->capRowDC );
   }
   cd->rowDC[row] = 0;

   /* for convenience if it is an inequality,
      lets put in the cx <= -cx + ub   format,  */
   if ( IS_GREATER_THAN( lb,ub ) )
   {
      int j;
      for ( j=0 ; (j<size) ; ++j )
         cd->coef[j] *= -1.0;
      const double u = lb;
      lb = ub;
      ub = u;

      lb *= -1.0;
      ub *= -1.0;
   }

   const int *idx = cd->idx;
   const double *coef = cd->coef;

   const double *colLB = cd->colLB;
   const double *colUB = cd->colUB;
   const char *isInteger = cd->isInteger;

   if ( ub == DBL_MAX )
      return 0;

   /* negative variables in the let side will change the bound */
   int j;
   for ( j=0 ; ( j<size ) ; ++j )
   {
      const int cIdx = idx[j];
      const double cCoef = coef[j];
      const double cLB = colLB[cIdx];
      const double cUB = colUB[cIdx];

      if ( ( fabs( cLB ) < EPS) && ( fabs( cUB ) < EPS ) )   /* fixed in zero, */
         continue;                                           /* does not matters */

      if ( ( cLB < -EPS ) && ( cUB > EPS ) )  /* not a positive only (or negative only) variable */
         return 0;

      /* coefficient in inequality is negative */
      if ( IS_NEGATIVE_DBL( cCoef ) )    /* there should not exist an positive unbounded variable with negative coefficient */
      {
         if ( cUB == DBL_MAX )
            return 0;
         else
         {
            if ( cUB > EPS )
            {
               double oldUB = ub;
               ub += cUB*cCoef*-1.0;
               if ( ub < oldUB )
                  return 0;
            }
         }
      }
      else
      {
         if ( cLB == -DBL_MAX )
            return 0;
         else
         {
            if ( cLB < -EPS )
            {
               double oldUB = ub;
               ub += cUB*cCoef*-1.0;
               if ( ub < oldUB )
                  return 0;
            }
         }
      }  /* positive coefficient */

   }

   /* now checking if it is a constraint with a large clique */
   int nColsTightUB = 0;
   int i;
   for ( i=0 ; (i<size) ; ++i )
   {
      const int cIdx = idx[i];

      const double cCoef = coef[ i ];

      if ( ( fabs( cCoef-ub ) <= EPS ) && ( isInteger[cIdx] ) )
      {
         cd->confIdx[ nColsTightUB ] = cIdx;
         ++nColsTightUB;
      }
   }

   if ( nColsTightUB >= MIN_CLIQUE_ROW )
   {
      cd->rowDC[row] += (nColsTightUB*(nColsTightUB-1))/2;
#ifdef VERBOSE_CGRAPH
      {
         char rName[256];
         fill_row_descr( cd, lp, row, rName );
         printf("cd-> added clique of size %d using row %s\n", nColsTightUB, rName );
      }
#endif
      cgraph_add_clique( cd->cgraph, cd->confIdx, nColsTightUB );

      return 0;
   }

   /* checking for pairwise conflicts */

   const int lastJ1 = size - 1;
   int j1;
   for ( j1=0 ; (j1<lastJ1) ; ++j1 )
   {
      const int cIdx1 = idx[j1];
      if ( !isInteger[cIdx1] )
         continue;

      const double partialBound = coef[ j1 ];

      int j2;
      for ( j2=j1+1 ; (j2<size) ; ++j2 )
      {
         const int cIdx2 = idx[ j2 ];

         if ( (( partialBound + coef[ j2 ] - ub) > EPS) && ( isInteger[ cIdx2 ] ) )
         {
            Conflict c;
            c.i = cIdx1;
            c.j = cIdx2;
            cd_check_cv_cap( cd, cd->cvecSize+1 );
            cd->cvec[cd->cvecSize] = c;
            cd->cvecSize++;
            cd->rowDC[row]++;
         }
      }
   }

#ifdef VERBOSE_CGRAPH
   if (cd->rowDC[row])
   {
      char rName[256];
      fill_row_descr( cd, lp, row, rName );
      printf( "row %s generated %d conflicts.\n", rName,  cd->rowDC[row] );
   }
#endif

   return cd->cvecSize;
}

void cd_add_conflicts_cgraph( ConflictDiscover *cd, Conflict *conflicts, int nConflicts )
{
   qsort( conflicts, nConflicts, sizeof(Conflict), cd_conflict_order_i );

#ifdef DEBUG
   int i;
   for ( i=0 ; (i<nConflicts-1) ; ++i )
   {
      assert( conflicts[i].i <= conflicts[i+1].i );
      if ( conflicts[i].i == conflicts[i+1].i )
         assert( conflicts[i].j <= conflicts[i+1].j );
   }
#endif

   int lastNode = -1;
   int nConfNode = 0;

   int *idx = cd->idx;
   CGraph *cgraph = cd->cgraph;

   {
      int i;
      for ( i=0 ; (i<nConflicts) ; ++i )
      {
         const Conflict *c = conflicts + i;
         if ( c->i != lastNode )
         {
            if ( nConfNode )
               cgraph_add_node_conflicts_no_sim( cgraph, lastNode, idx, nConfNode );

            lastNode = c->i;
            idx[ 0 ] = c->j;
            nConfNode = 1;
         }
         else
         {
            if ((nConfNode) && (c->j==idx[nConfNode-1]))
               continue;
            idx[ nConfNode ] = c->j;
            nConfNode++;
         }
      }
   }

   if ( nConfNode )
      cgraph_add_node_conflicts_no_sim( cgraph, lastNode, idx, nConfNode );

   qsort( conflicts, nConflicts, sizeof(Conflict), cd_conflict_order_j );

   lastNode = -1;
   nConfNode = 0;

   {
      int i;
      for ( i=0 ; (i<nConflicts) ; ++i )
      {
         const Conflict *c = conflicts + i;
         if ( c->j != lastNode )
         {
            if ( nConfNode )
               cgraph_add_node_conflicts_no_sim( cgraph, lastNode, idx, nConfNode );

            lastNode = c->j;
            idx[ 0 ] = c->i;
            nConfNode = 1;

         }
         else
         {
            if ((nConfNode) && (c->i==idx[nConfNode-1]))
               continue;
            idx[ nConfNode ] = c->i;
            nConfNode++;
         }
      }
   }

   if ( nConfNode )
      cgraph_add_node_conflicts_no_sim( cgraph, lastNode, idx, nConfNode );
}

void cd_check_cv_cap( ConflictDiscover *cd, int size )
{
   if ( size > cd->cvecCap )
   {
      cd->cvecCap = ((double)cd->cvecCap)*1.5;
      if ( size > cd->cvecCap )
         cd->cvecCap = size;

      cd->cvec = xrealloc( cd->cvec, sizeof(Conflict)*cd->cvecCap );
   }
}

const int *cd_get_row_dc( ConflictDiscover *cd )
{
   return cd->rowDC;
}

char cd_col_integer( ConflictDiscover *cd, int col )
{
   return cd->isInteger[col];
}

double cd_col_lb( ConflictDiscover *cd, int col )
{
   return cd->colLB[col];
}

double cd_col_ub( ConflictDiscover *cd, int col )
{
   return cd->colUB[col];
}

void cd_free( ConflictDiscover **cd )
{
   free( (*cd)->colLB );
   free( (*cd)->colUB );
   free( (*cd)->isInteger );
   free( (*cd)->idx );
   free( (*cd)->confIdx );
   free( (*cd)->coef );
   free( (*cd)->cvec );
   free( (*cd)->rowDC );
   free( (*cd)->sense );
   free( (*cd)->rhs );

   free( *cd );
   *cd = NULL;
}
