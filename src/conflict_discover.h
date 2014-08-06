#ifndef CONFLICT_DISCOVER_H_INCLUDED
#define CONFLICT_DISCOVER_H_INCLUDED

#include "cgraph.h"
#include "lp_callbacks.h"

typedef struct _ConflictDiscover ConflictDiscover;

ConflictDiscover *cd_create(
   CGraph *_cgraph,
   const void *_lp,
   LPCallbacks *lpc
   );

void cd_update_graph( ConflictDiscover *cd, const void *lp );

void cd_free( ConflictDiscover **cd );


/* "Private functions" */

void cd_update_col_info( ConflictDiscover *cd, void *lp, int start, int end );

int cd_process_row( ConflictDiscover *cd, void *lp, int row );

int cd_conflict_order_i( const void *e1, const void *e2 );

/* per row discoverd conflicts */
const int *cd_get_row_dc( ConflictDiscover *cd );

/* col properties functions */
char cd_col_integer( ConflictDiscover *cd, int col );
double cd_col_lb( ConflictDiscover *cd, int col );
double cd_col_ub( ConflictDiscover *cd, int col );

/* make sure that conflict vector can hold "size" elements */
void cd_check_cv_cap( ConflictDiscover *cd, int size );


#endif
