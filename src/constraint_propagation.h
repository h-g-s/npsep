#ifndef CONSTRAINT_PROPAGATION_H
#define CONSTRAINT_PROPAGATION_H

extern "C"
{
	#include "cgraph.h"
}

#include "problem.h"
#include <vector>

using namespace std;

#define DEACTIVATE 0
#define ACTIVATE 1
#define UNFIXED 2
#define CONFLICT 3
#define FIXATION 4
#define NOIMPLICATION 5
#define EPS 1e-8

typedef struct Fixation
{
	int idxVar;
	double valueToFix;

    Fixation() { }
    Fixation(int _idxVar, double _valueToFix) : idxVar(_idxVar), valueToFix(_valueToFix) { }
} Fixation;

inline bool operator<(const Fixation &f1, const Fixation &f2)
{
	if(f1.idxVar != f2.idxVar)
		return f1.idxVar < f2.idxVar;
	return f1.valueToFix < f2.valueToFix;
}

/* stores constraints which have only binary variables */
typedef struct CPropagation
{
	int numCols, numRows;
	char *varIsBinary; /* for each variable, stores 1 if variable is binary and 0 otherwise */

	vector<vector<pair<int, double> > > matrixByRow; /* matrix of indexes and coefficients by row */
    vector<vector<pair<int, double> > > matrixByCol; /* matrix of indexes and coefficients by column */
	vector<double> rhs; /* right-hand side for each constraint */
	int binaryVars; /* number of binary variables and number of unfixed binary variables */

	/* Used for preprocess module. For each variable:
		0 - if variable must be fixed to 0
        1 - if variable must be fixed to 1
        2 - if variable is unfixed */
	char *isToFix;
    int varsToFix; /* number of variables to fix */

    Problem *problem;
} CPropagation;

CPropagation* cpropagation_create(const Problem *problem);
void cpropagation_free(CPropagation *cp);
void cpropagation_get_vars_to_fix(CPropagation *cp, const CGraph* cgraph);
int cpropagation_get_num_vars_to_fix(CPropagation *cp);
char cpropagation_var_is_to_fix(CPropagation *cp, int idxVar);
OsiSolverInterface* cpropagation_preprocess(CPropagation *cp, int nindexes[]);

char constraintPropagation(CPropagation *cp, int var, double value, vector<Fixation> &fixations,
                            int &unfixedVars, double *colLb, double *colUb, double *constrBound, int *unfixedVarsByRow);
void unfixVariable(CPropagation *cp, int idx, int &unfixedVars, double *colLb, double *colUb, double *constrBound,
                    int *unfixedVarsByRow);
void fixVariable(CPropagation *cp, int idx, double value, int &unfixedVars, double *colLb, double *colUb, double *constrBound,
                    int *unfixedVarsByRow);
char evaluateBound(const CPropagation *cp, const pair<int, double> &var, int constraint, double *constrBound);

#endif