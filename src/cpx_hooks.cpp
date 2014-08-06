#include <cplex.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cfloat>

extern "C"
{
#include "lp_callbacks.h"
}

static CPXENVptr LPcpxDefaultEnv = NULL;

void cpx_check_error(CPXENVptr env, int errorCode, const char *sourceFile, int sourceLine);

void cpx_check_cpx_env()
{
   if (!LPcpxDefaultEnv)
   {
      int cpxError = 0;
      LPcpxDefaultEnv = CPXopenCPLEX( &cpxError );
      if (!LPcpxDefaultEnv)
      {
         fprintf( stderr, "Error opening CPLEX environment.\n" );
         exit( EXIT_FAILURE );
      }
   }
}

int cpx_get_num_cols( void *lp )
{
   assert( lp!=NULL ); cpx_check_cpx_env(); CPXLPptr cpxLP = (CPXLPptr) lp;
   return CPXgetnumcols(LPcpxDefaultEnv, cpxLP );
}

int cpx_get_num_rows( void *lp )
{
   assert( lp!=NULL ); cpx_check_cpx_env(); CPXLPptr cpxLP = (CPXLPptr) lp;
   return CPXgetnumrows(LPcpxDefaultEnv, cpxLP );
}

void cpx_get_col_info (void *lp, int col, double *lb, double *ub, char *isInteger)
{
   assert( lp!=NULL ); cpx_check_cpx_env(); CPXLPptr cpxLP = (CPXLPptr) lp;
   assert( (col>=0) && (col<cpx_get_num_cols(lp)) );
   char cType = 0; int cpxError = 0;
   
   cpxError = CPXgetlb( LPcpxDefaultEnv, cpxLP, lb, col, col );
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );
   cpxError = CPXgetub( LPcpxDefaultEnv, cpxLP, ub, col, col );
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );
   cpxError = CPXgetctype( LPcpxDefaultEnv, cpxLP, &cType, col, col );
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );
   *isInteger = ((cType==CPX_BINARY)||(cType==CPX_INTEGER));
}

void cpx_get_row_info (void *lp, int row, double *lb, double *ub, int *size, int *idx, double *coef)
{
   assert( lp!=NULL ); assert( (row>=0) && (row<cpx_get_num_rows(lp)) );

   int cpxError = 0; char sense = 0; double rhs = 0.0;
   *lb = -DBL_MAX; *ub = DBL_MAX; 
   
   cpx_check_cpx_env();
   CPXLPptr cpxLP = (CPXLPptr) lp;
   
   cpxError = CPXgetrhs( LPcpxDefaultEnv, cpxLP, &rhs, row, row );
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );
   
   cpxError = CPXgetsense( LPcpxDefaultEnv, cpxLP, &sense, row, row );
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );
   switch (sense)
   {
      case 'E':
         *lb = *ub = rhs;
         break;
      case 'G':
         *lb = rhs;
         break;
      case 'L':
         *ub = rhs;
         break;
      default:
         fprintf( stderr, "Error: cpx_get_row_info %s:%d - not ready to handle ranged constraints.  fix it.\n", __FILE__, __LINE__ );
         abort();
   }

   int matbeg;
   int surplus = 0;
   cpxError = CPXgetrows( LPcpxDefaultEnv, cpxLP, size, &matbeg, idx, coef, cpx_get_num_cols(lp) + 1, &surplus, row, row);
   assert( surplus>=0 );
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );
}

const char *cpx_get_row_name( void *lp, int row )
{
   assert( lp!=NULL ); assert( (row>=0) && (row<cpx_get_num_rows(lp)) );

   cpx_check_cpx_env();
   CPXLPptr cpxLP = (CPXLPptr) lp;
   static char cstore[512];
   static char *namev[] = { NULL };
   int sp = 0;
   
   int cpxError = CPXgetrowname( LPcpxDefaultEnv, cpxLP, namev, cstore, 512, &sp, row, row );
   if (sp<0)
   {
      fprintf( stderr, "Not enough space to store row name.\n" );
      abort();
   }
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );

   return &(cstore[0]);
}

const char *cpx_get_col_name( void *lp, int col )
{
   assert( lp!=NULL ); assert( (col>=0) && (col<cpx_get_num_cols(lp)) );

   cpx_check_cpx_env();
   CPXLPptr cpxLP = (CPXLPptr) lp;
   static char cstore[512];
   static char *namev[] = { NULL };
   int sp = 0;
   
   int cpxError = CPXgetcolname( LPcpxDefaultEnv, cpxLP, namev, cstore, 512, &sp, col, col );
   if (sp<0)
   {
      fprintf( stderr, "Not enough space to store col name.\n" );
      abort();
   }
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );

   return &(cstore[0]);
}

double cpx_get_col_dual( void *lp, int col )
{
   double res = 0.0;
   assert( lp!=NULL ); assert( (col>=0) && (col<cpx_get_num_cols(lp)) );

   cpx_check_cpx_env();
   CPXLPptr cpxLP = (CPXLPptr) lp;
   
   int cpxError = CPXgetdj( LPcpxDefaultEnv, cpxLP, &res, col, col );
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );

   return res;
}

void cpx_get_row_sense( void *lp, char sense[] )
{
   assert( lp!=NULL );
   cpx_check_cpx_env();
   CPXLPptr cpxLP = (CPXLPptr) lp;
   
   int cpxError = CPXgetsense(LPcpxDefaultEnv, cpxLP, sense, 0, cpx_get_num_rows(lp)-1 );
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );
}

void cpx_get_rhs( void *lp, double rhs[] )
{
   assert( lp!=NULL );
   cpx_check_cpx_env();
   CPXLPptr cpxLP = (CPXLPptr) lp;
   
   int cpxError = CPXgetrhs(LPcpxDefaultEnv, cpxLP, rhs, 0, cpx_get_num_rows(lp)-1 );
   cpx_check_error( LPcpxDefaultEnv, cpxError, __FILE__, __LINE__ );
}

void cpx_ini_lp_callbacks( LPCallbacks *lpc )
{
   lpc->cb_get_num_cols = cpx_get_num_cols;
   lpc->cb_get_num_rows = cpx_get_num_rows;

   lpc->cb_get_col_info = cpx_get_col_info;
   lpc->cb_get_row_info = cpx_get_row_info;

   lpc->cb_get_col_dual = cpx_get_col_dual;

   lpc->cb_get_col_name = cpx_get_col_name;
   lpc->cb_get_row_name = cpx_get_row_name;
   lpc->cb_get_row_sense = cpx_get_row_sense;
   lpc->cb_get_rhs = cpx_get_rhs;
}

void cpx_check_error(CPXENVptr env, int errorCode, const char *sourceFile, int sourceLine)
{
    if (errorCode) {
        char errorStr[256];
        CPXgeterrorstring(LPcpxDefaultEnv, errorCode, errorStr);
        fprintf(stderr, "CPLEX Error: %s\n", errorStr);
        fprintf(stderr, "Inside LP Library - %s:%d\n\n", sourceFile, sourceLine);
        abort();
        exit(EXIT_FAILURE);
    }
}

