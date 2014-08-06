#include <OsiSolverInterface.hpp>
#include <CoinPackedMatrix.hpp>

extern "C"
{
#include "lp_callbacks.h"
}

int cbc_get_num_cols( void *lp )
{
   OsiSolverInterface *solver = (OsiSolverInterface *) lp;
   return solver->getNumCols();
}

int cbc_get_num_rows( void *lp )
{
   OsiSolverInterface *solver = (OsiSolverInterface *) lp;
   return solver->getNumRows();
}

void cbc_get_col_info (void *lp, int col, double *lb, double *ub, char *isInteger)
{
   OsiSolverInterface *solver = (OsiSolverInterface *) lp;
   *lb = solver->getColLower()[col];
   *ub = solver->getColUpper()[col];
   if ( solver->isInteger( col ) )
      *isInteger = 1;
   else
      *isInteger = 0;
}

void cbc_get_row_info (void *lp, int row, double *lb, double *ub, int *size, int *idx, double *coef)
{
   OsiSolverInterface *solver = (OsiSolverInterface *) lp;
   *lb = solver->getRowLower()[row];
   *ub = solver->getRowUpper()[row];
   const CoinPackedMatrix *cpm = solver->getMatrixByRow();
   *size = cpm->getVectorSize( row );
   for ( CoinBigIndex i = cpm->getVectorFirst( row ), p = 0 ; (i<cpm->getVectorLast( row )) ; ++i,++p  )
      idx[p] = cpm->getIndices( )[i];
   for ( CoinBigIndex i = cpm->getVectorFirst( row ), p = 0 ; (i<cpm->getVectorLast( row )) ; ++i,++p  )
      coef[p] = cpm->getElements( )[i];
}

const char *cbc_get_row_name( void *lp, int row )
{
   OsiSolverInterface *solver = (OsiSolverInterface *) lp;

   return solver->getRowName( row ).c_str();
}

const char *cbc_get_col_name( void *lp, int col )
{
   OsiSolverInterface *solver = (OsiSolverInterface *) lp;

   return solver->getColName( col ).c_str();
}

double cbc_get_col_dual( void *lp, int col )
{
   OsiSolverInterface *solver = (OsiSolverInterface *) lp;

   return solver->getReducedCost()[col];
}

void cbc_get_row_sense( void *lp, char sense[] )
{
   OsiSolverInterface *solver = (OsiSolverInterface *) lp;
   memcpy( sense, solver->getRowSense(), sizeof(char)*solver->getNumRows() );
}

void cbc_get_rhs( void *lp, double rhs[] )
{
   OsiSolverInterface *solver = (OsiSolverInterface *) lp;
   memcpy( rhs, solver->getRightHandSide() , sizeof(double)*solver->getNumRows() );
}

void cbc_ini_lp_callbacks( LPCallbacks *lpc )
{
   lpc->cb_get_num_cols = cbc_get_num_cols;
   lpc->cb_get_num_rows = cbc_get_num_rows;

   lpc->cb_get_col_info = cbc_get_col_info;
   lpc->cb_get_row_info = cbc_get_row_info;

   lpc->cb_get_col_dual = cbc_get_col_dual;

   lpc->cb_get_col_name = cbc_get_col_name;
   lpc->cb_get_row_name = cbc_get_row_name;
   lpc->cb_get_row_sense = cbc_get_row_sense;
   lpc->cb_get_rhs = cbc_get_rhs;
}

