#ifndef CPX_HOOKS_HPP_INCLUDED
#define CPX_HOOKS_HPP_INCLUDED

/// all this functions receive an OsiClpSolverInterface
int cpx_get_num_cols( void *lp );
int cpx_get_num_rows( void *lp );
void cpx_get_col_info (void *lp, int col, double *lb, double *ub, char *isInteger);
void cpx_get_row_info (void *lp, int row, double *lb, double *ub, int *size, int *idx, double *coef);
double cpx_get_col_dual( void *lp, int col );
const char *cpx_get_row_name( void *lp, int row );
const char *cpx_get_col_name( void *lp, int col );
void cpx_get_row_sense( void *lp, char sense[] );
void cpx_get_rhs( void *lp, double rhs[] );

void cpx_ini_lp_callbacks( LPCallbacks *lpc );


#endif // CPX_HOOKS_HPP_INCLUDED

