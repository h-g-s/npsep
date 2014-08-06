#include <stdlib.h>
#include "lp_callbacks.h"
#include "memory.h"

void lp_callbacks_init( LPCallbacks *lpc )
{
   lpc->cb_get_num_rows   = NULL;
   lpc->cb_get_num_cols   = NULL;
   lpc->cb_get_col_info   = NULL;
   lpc->cb_get_row_info   = NULL;
   lpc->cb_get_col_name   = NULL;
   lpc->cb_get_row_name   = NULL;
   lpc->cb_get_col_dual   = NULL;
   lpc->cb_add_cut        = NULL;
   lpc->cb_get_row_sense  = NULL;
   lpc->cb_get_rhs        = NULL;
}

void lp_update_idxs_to_glpk( int idx[], int size )
{
   int i;
   for ( i=0 ; (i<size) ; ++i )
      idx[i]++;
}

const double *lp_get_ones_vector( int size )
{
   static double *vOnes = NULL;
   static int vOnesCap = 0;
   if ( vOnesCap < size )
   {
      if ( vOnesCap == 0 )
      {
         vOnesCap  = ((size+1)*2)|(1024);
         vOnes = xmalloc( sizeof(double)*vOnesCap );
      }
      else
      {
         vOnesCap = (size+1) * 2;
         vOnes = xrealloc( vOnes, sizeof(double)*vOnesCap );
      }

      int i;
      for ( i=0 ; (i<vOnesCap) ; ++i )
         vOnes[i] = 1.0;
   }

   return vOnes;
}

