#ifndef CBC_HOOKS_HPP_INCLUDED
#define CBC_HOOKS_HPP_INCLUDED

/// all this functions receive an OsiClpSolverInterface
int cbc_get_num_cols( void *lp );
int cbc_get_num_rows( void *lp );
void cbc_get_col_info (void *lp, int col, double *lb, double *ub, char *isInteger);
void cbc_get_row_info (void *lp, int row, double *lb, double *ub, int *size, int *idx, double *coef);
double cbc_get_col_dual( void *lp, int col );
const char *cbc_get_row_name( void *lp, int row );
const char *cbc_get_col_name( void *lp, int col );
void cbc_get_row_sense( void *lp, char sense[] );
void cbc_get_rhs( void *lp, double rhs[] );

void cbc_ini_lp_callbacks( LPCallbacks *lpc );


#endif // CBC_HOOKS_HPP_INCLUDED

