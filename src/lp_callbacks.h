#ifndef LP_CALLBACKS
#define LP_CALLBACKS

typedef struct
{
   int  (*cb_get_num_cols) (void *lp);
   int  (*cb_get_num_rows) (void *lp);
   void (*cb_get_col_info) (void *lp, int col, double *lb, double *ub, char *isInteger);
   void (*cb_get_row_info) (void *lp, int row, double *lb, double *ub, int *size, int *idx, double *coef);
   const char *(*cb_get_col_name) (void *lp, int col);
   const char *(*cb_get_row_name) (void *lp, int row);
   double (*cb_get_col_dual) (void *lp, int col);
   void (*cb_add_cut) (void *lp, int nz, const int *idx, const double *coef, double rhs );
   void (*cb_get_row_sense) ( void *lp, char sense[] );
   void (*cb_get_rhs) ( void *lp, double *rhs );
} LPCallbacks;

void lp_callbacks_init( LPCallbacks *lpc );

const double *lp_get_ones_vector( int size );

/* increments by one every index */
void lp_update_idxs_to_glpk( int idx[], int size );

#endif

