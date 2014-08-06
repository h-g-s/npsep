#include <glpk.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "glp_hooks.h"

int cb_glp_get_num_cols( void *lp )
{
   return glp_get_num_cols( lp );
}

int cb_glp_get_num_rows( void *lp )
{
   return glp_get_num_rows( lp );
}

void cb_glp_get_col_info( void *lp, int col, double *lb, double *ub, char *isInteger )
{
   const int cIdx = col+1;
   const int kind = glp_get_col_kind( lp, cIdx );
   *isInteger =  ( ( kind == GLP_IV ) || ( kind == GLP_BV ) );
   *lb = glp_get_col_lb( lp, cIdx );
   *ub = glp_get_col_ub( lp, cIdx );
}

void cb_glp_get_row_info( void *lp, int row, double *lb, double *ub, int *size, int *idx, double *coef )
{
   const int rIdx = row + 1;

   const int rType = glp_get_row_type( lp, rIdx );

   switch ( rType )
   {
   case GLP_FX :
      *lb = *ub = glp_get_row_lb( lp , rIdx );
      break;
   case GLP_DB:
      *lb = glp_get_row_lb( lp , rIdx );
      *ub = glp_get_row_ub( lp , rIdx );
      break;
   case GLP_LO :
      *lb = glp_get_row_lb( lp , rIdx );
      *ub = DBL_MAX;
      break;
   case GLP_UP :
      *lb = -DBL_MAX;
      *ub = glp_get_row_ub( lp , rIdx );
      break;
   default:
      fprintf( stderr, "Not prepared to handle." );
      exit( EXIT_FAILURE );
   }

   *size = glp_get_mat_row( lp, rIdx, idx-1, coef-1 );
   int i;
   for ( i=0 ; (i<*size) ; ++i )
      (idx[i])--;
}

const char *cb_glp_get_col_name( void *lp, int col )
{
   return glp_get_col_name( lp, col+1 );
}

const char *cb_glp_get_row_name( void *lp, int row )
{
   return glp_get_row_name( lp, row+1 );
}

double cb_glp_get_col_dual( void *lp, int col )
{
   return glp_get_col_dual( lp, col+1 );
}

/* in glpk this function receives a tree instead of lp */
void cb_glp_add_cut( void *tree, int nz, const int *idx, const double *coef, double rhs )
{
   glp_ios_add_row( tree, NULL, 0, 0, nz, idx, coef, GLP_UP, rhs );
}

void cd_ini_lp_callbacks( LPCallbacks *lpc )
{
   lpc->cb_get_col_info = cb_glp_get_col_info;
   lpc->cb_get_row_info = cb_glp_get_row_info;
   lpc->cb_get_num_rows = cb_glp_get_num_rows;
   lpc->cb_get_num_cols = cb_glp_get_num_cols;
   lpc->cb_get_col_name = cb_glp_get_col_name;
   lpc->cb_get_row_name = cb_glp_get_row_name;
   lpc->cb_get_col_dual = cb_glp_get_col_dual;
   lpc->cb_add_cut      = cb_glp_add_cut;
}


