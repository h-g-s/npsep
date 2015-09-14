/*
	Structure Problem is used to store MIPs
	without dependence on solvers.
*/

#ifndef PROBLEM_H
#define PROBLEM_H

#include <OsiSolverInterface.hpp>

#define CONTINUOUS 0
#define BINARY 1
#define INTEGER 2

typedef struct _Problem Problem;

Problem* problem_create(int numCols, int numRows, double infty);
Problem* problem_create_using_osi(const OsiSolverInterface *solver); /* Fills the structure Problem using an OsiSolverInterface pointer */
void problem_free(Problem **p);

/* general information */
int problem_num_cols(const Problem *p);
int problem_num_rows(const Problem *p);
int problem_num_elements(const Problem *p);
double problem_get_infinity(const Problem *p);

/* information about constraints */
int problem_row_size(const Problem *p, int idxRow);
char problem_row_sense(const Problem *p, int idxRow);
const int* problem_row_idxs(const Problem *p, int idxRow);
const double* problem_row_coefs(const Problem *p, int idxRow);
double problem_row_rhs(const Problem *p, int idxRow);
const char* problem_row_name(const Problem *p, int idxRow);

/* information about variables */
char problem_var_is_binary(const Problem *p, int idxVar);
char problem_var_type(const Problem *p, int idxVar);
const char* problem_var_name(const Problem *p, int idxVar);
double problem_var_lower_bound(const Problem *p, int idxVar);
double problem_var_upper_bound(const Problem *p, int idxVar);
double problem_var_obj_coef(const Problem *p, int idxVar);
int problem_var_n_rows(const Problem *p, int idxVar);
const int* problem_var_rows_idxs(const Problem *p, int idxVar);
const double* problem_var_rows_coefs(const Problem *p, int idxVar);
const double* problem_vars_lower_bound(const Problem *p);
const double* problem_vars_upper_bound(const Problem *p);
const double* problem_vars_obj_coefs(const Problem *p);

/* general set methods */
void problem_set_num_cols(Problem *p, double numCols);
void problem_set_num_rows(Problem *p, double numRows);
void problem_set_num_elements(Problem *p, int nElements);
void problem_set_infinity(Problem *p, double infty);

/* set methods for variables */
void problem_var_set_lower_bound(Problem *p, int idxVar, double value);
void problem_var_set_upper_bound(Problem *p, int idxVar, double value);
void problem_var_set_obj_coef(Problem *p, int idxVar, double value);
void problem_var_set_type(Problem *p, int idxVar, char type);
void problem_var_set_name(Problem *p, int idxVar, const char *name);

/* set methods for constraints */
void problem_set_row(Problem *p, int idxRow, const int *idxs, const double *coefs, int nElements, double rhs, char sense);
void problem_row_set_name(Problem *p, int idxRow, const char *name);
/* used after call problem_set_row (after filling the matrices) */
void problem_update_matrices_by_col(Problem *p);

/* converts Problem structure to OsiSolverInterface */
OsiSolverInterface* problem_convert_to_osi(Problem *p);

#endif