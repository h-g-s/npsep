#ifndef GLP_HOOKS_H_INCLUDED
#define GLP_HOOKS_H_INCLUDED

/**
 * set of functions needed to discover
 * LP properties for a glpk problem
 **/

#include "lp_callbacks.h"

int cb_glp_get_num_cols( void *lp );
int cb_glp_get_num_rows( void *lp );
void cb_glp_get_col_info( void *lp, int col, double *lb, double *ub, char *isInteger );
void cb_glp_get_row_info( void *lp, int row, double *lb, double *ub, int *size, int *idx, double *coef );
const char *cb_glp_get_col_name( void *lp, int col );
const char *cb_glp_get_row_name( void *lp, int row );
double cb_glp_get_col_dual( void *lp, int col );
void cb_glp_add_cut( void *tree, int nz, const int *idx, const double *coef, double rhs );

void cd_ini_lp_callbacks( LPCallbacks *lpc );

#endif

