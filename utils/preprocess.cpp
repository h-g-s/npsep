#include "preprocess.h"
#include <cmath>
#include <vector>
#include <queue>
#include "constraint_propagation.h"

using namespace std;

extern "C"
{
    #include "memory.h"
    #include "cgraph.h"
}

#define EPS 1e-8

struct _Preprocess
{
    const Problem *problem;

    double **coefficients, *rhs;
    double *colLb, *colUb;
    int *colNConstraints; /* count how many constraints each variable appears  */
    char *rowSense;

    double *lhsMin, *lhsMax; /* minimum and maximum values for LHS of each constraint. */
    int *countInftyPosBounds, *countInftyNegBounds; /* counts how many variables at each constraint have infinity as upper bound */
    char *negBounds;  /* marks rows containing variables with negative bounds. (used in preprocess step) */
    int *nindexes; /* new indexes of variables (-1 if variable does not appear at preprocessed problem) */
    char *removeRow; /* redundant rows will be deleted. */
    int rowsRemoved; /* redundant rows removed at preprocessing */
    int varsRemoved; /* number of variables fixed at preprocessing */
    int coefsImproved; /* number of coefficients improved at preprocessing */
    int boundsImproved; /* number of bounds of variables improved at preprocessing */

    /*----------------------------------CONSTRAINT PROPAGATION AND GUB STUFF---------------------------------------------*/
    char *isGubConstraint; /* for each constraint, stores 1 if it is a GUB constraint. 0 otherwise */
    CPropagation *cp;
    double *cpColLb, *cpColUb, *constrBound;
    int unfixedVars, *unfixedVarsByRow;
    int numGub, transformed;
};

void gub_constraint(Preprocess *pp, int idxRow);
void handling_inequality_constraints(Preprocess *pp, int idxRow);
void handling_equality_constraints(Preprocess *pp, int idxRow);
void one_element_constraint(Preprocess *pp, int idxRow);
void fix_var(Preprocess *pp, int idxVar, double valueToFix);
void improve_lower_bound(Preprocess *pp, int idxVar, double newBound);
void improve_upper_bound(Preprocess *pp, int idxVar, double newBound);
void execute_basic_preprocessing(Preprocess *pp);

Preprocess* preprocess_create(const Problem *problem)
{
    int i, j, nCols = problem_num_cols(problem), nRows = problem_num_rows(problem);
    Preprocess* pp = new Preprocess;
    pp->cp = cpropagation_create(problem);

    pp->problem = problem;
    pp->coefficients = (double**) xmalloc(sizeof(double*) * nRows);
    pp->rhs = (double*) xmalloc(sizeof(double) * nRows);
    pp->colLb = (double*) xmalloc(sizeof(double) * nCols);
    pp->colUb = (double*) xmalloc(sizeof(double) * nCols);
    pp->colNConstraints = (int*) xmalloc(sizeof(int) * nCols);
    pp->rowSense = (char*) xmalloc(sizeof(char) * nRows);
    pp->lhsMin = (double*) xmalloc(sizeof(double) * nRows);
    pp->lhsMax = (double*) xmalloc(sizeof(double) * nRows);
    pp->countInftyPosBounds = (int*) xmalloc(sizeof(int) * nRows);
    pp->countInftyNegBounds = (int*) xmalloc(sizeof(int) * nRows);
    pp->negBounds = (char*) xmalloc(sizeof(char) * nRows);
    pp->nindexes = (int*) xmalloc(sizeof(int) * nCols);
    pp->removeRow = (char*) xmalloc(sizeof(char) * nRows);
    pp->rowsRemoved = pp->varsRemoved = pp->coefsImproved = pp->boundsImproved = 0;
    pp->isGubConstraint = (char*) xmalloc(sizeof(char) * nRows);
    pp->numGub = 0;

    for(i = 0; i < nCols; i++)
    {
        pp->colLb[i] = problem_var_lower_bound(pp->problem, i);
        pp->colUb[i] = problem_var_upper_bound(pp->problem, i);
        pp->nindexes[i] = -1;
        pp->colNConstraints[i] = problem_var_n_rows(pp->problem, i);

        if(pp->colLb[i] == pp->colUb[i])
            pp->varsRemoved++;
    }

    const double infty = problem_get_infinity(pp->problem);

    int count = 0;

    for(i = 0; i < nRows; i++)
    {
        const int *idxs = problem_row_idxs(pp->problem, i);
        const double *coefs = problem_row_coefs(pp->problem, i);
        const char rowSense = problem_row_sense(pp->problem, i);
        const int rowSize = problem_row_size(pp->problem, i);
        const double mult = (rowSense == 'G') ? -1.0 : 1.0; /* used to convert >= rows to <= */
        double minCoef = infty;
        double maxCoef = -infty;
        int countBin = 0;

        pp->rowSense[i] = rowSense;

        pp->removeRow[i] = 0;
        pp->rhs[i] = mult * problem_row_rhs(pp->problem, i);
        pp->coefficients[i] = (double*) xmalloc(sizeof(double) * rowSize);
        pp->lhsMin[i] = pp->lhsMax[i] = 0.0;
        pp->countInftyPosBounds[i] = pp->countInftyNegBounds[i] = 0;
        pp->negBounds[i] = 0;

        for(j = 0; j < rowSize; j++)
        {
        	if(problem_var_is_binary(pp->problem, idxs[j]))
        		countBin++;

        	minCoef = min(minCoef, coefs[j]);
        	maxCoef = max(maxCoef, coefs[j]);

            pp->coefficients[i][j] = mult * coefs[j];

            if(pp->negBounds[i]) continue;

            if(pp->colLb[idxs[j]] <= -EPS) /* variable with negative bounds */
            {
                pp->negBounds[i] = 1;
                pp->lhsMin[i] = pp->lhsMax[i] = 0.0;
                pp->countInftyPosBounds[i] = pp->countInftyNegBounds[i] = 0;
                continue;
            }

            if(pp->colUb[idxs[j]] >= infty)
                (pp->coefficients[i][j] <= -EPS) ? pp->countInftyNegBounds[i]++ : pp->countInftyPosBounds[i]++;
            else
            {
                if(pp->coefficients[i][j] <= -EPS)
                {
                    pp->lhsMin[i] += (pp->colUb[idxs[j]] * pp->coefficients[i][j]);
                    pp->lhsMax[i] += (pp->colLb[idxs[j]] * pp->coefficients[i][j]);
                }
                else
                {
                    pp->lhsMin[i] += (pp->colLb[idxs[j]] * pp->coefficients[i][j]);
                    pp->lhsMax[i] += (pp->colUb[idxs[j]] * pp->coefficients[i][j]);
                }
            }
        }

        //Detects GUB constraints:
        if( (countBin == rowSize) && (fabs(maxCoef - minCoef) < EPS) && (fabs(maxCoef - pp->rhs[i]) < EPS) &&
           ( ( (pp->rowSense[i] == 'L') && (maxCoef > EPS) ) || ( (pp->rowSense[i] == 'G') && (maxCoef < EPS) ) ) )
        {
        	pp->isGubConstraint[i] = 1;
        	pp->numGub++;
        }
        else pp->isGubConstraint[i] = 0;
    }

/*----------------------------------CONSTRAINT PROPAGATION AND GUB STUFF---------------------------------------------*/
    int cpRows = pp->cp->numRows;

    pp->cpColLb = (double*) xmalloc(sizeof(double) * nCols);
    pp->cpColUb = (double*) xmalloc(sizeof(double) * nCols);
    pp->constrBound = (double*) xmalloc(sizeof(double) * cpRows);
    pp->unfixedVars = problem_num_binaries(pp->problem);
    pp->unfixedVarsByRow = (int*) xmalloc(sizeof(int) * cpRows);
    pp->transformed = 0;

    for(int i = 0; i < nCols; i++)
    {
        pp->cpColLb[i] = pp->colLb[i];
        pp->cpColUb[i] = pp->colUb[i];
    }
    for(int i = 0; i < pp->cp->numRows; i++)
    {
        pp->constrBound[i] = pp->cp->rhs[i];
        pp->unfixedVarsByRow[i] = (int)pp->cp->matrixByRow[i].size();

        for(int j = 0; j < (int)pp->cp->matrixByRow[i].size(); j++)
        {
            const int idx = pp->cp->matrixByRow[i][j].first;
            const double coef = pp->cp->matrixByRow[i][j].second;

            if(coef > 0.0) pp->constrBound[i] -= (coef * pp->cpColLb[idx]);
            else pp->constrBound[i] -= (coef * pp->cpColUb[idx]);
        }
    }

    return pp;
}

void preprocess_free(Preprocess **pp)
{
    int i, nRows = problem_num_rows((*pp)->problem);
    for(i = 0; i < nRows; i++) free((*pp)->coefficients[i]);
    free((*pp)->coefficients);
    free((*pp)->rhs);
    free((*pp)->colLb);
    free((*pp)->colUb);
    free((*pp)->colNConstraints);
    free((*pp)->rowSense);
    free((*pp)->lhsMin);
    free((*pp)->lhsMax);
    free((*pp)->countInftyPosBounds);
    free((*pp)->countInftyNegBounds);
    free((*pp)->negBounds);
    free((*pp)->nindexes);
    free((*pp)->removeRow);

    free((*pp)->isGubConstraint);
    free((*pp)->cpColLb);
    free((*pp)->cpColUb);
    free((*pp)->constrBound);
    free((*pp)->unfixedVarsByRow);
    cpropagation_free((*pp)->cp);
    

    free(*pp);
    (*pp) = NULL;
}

void execute_basic_preprocessing(Preprocess *pp)
{
    int i, j;
    int nCols = problem_num_cols(pp->problem), nRows = problem_num_rows(pp->problem);
    double infty = problem_get_infinity(pp->problem);

    for(i = 0; i < nRows; i++)
    {
        int rowSize = problem_row_size(pp->problem, i);
        const int *idxs = problem_row_idxs(pp->problem, i);
        char sense = pp->rowSense[i];

        if(sense != 'E' && sense != 'L' && sense != 'G')
        {
            fprintf(stderr, "Error: invalid type of constraint!\n");
            exit(EXIT_FAILURE);
        }

        if(rowSize == 1)
        {
            one_element_constraint(pp, i);
            continue;
        }

        if(pp->negBounds[i]) continue; /* ignoring rows containing variables with negative bounds */
        if(pp->removeRow[i]) continue; /* ignoring rows already removed */

        if(pp->isGubConstraint[i])
        	gub_constraint(pp, i);

        if(sense == 'L' || sense == 'G')
        	handling_inequality_constraints(pp, i);
        else
        	handling_equality_constraints(pp, i);
    }
}

void gub_constraint(Preprocess *pp, int idxRow)
{
	char status;
	int rowSize = problem_row_size(pp->problem, idxRow);
    const int *idxs = problem_row_idxs(pp->problem, idxRow);
    vector<Fixation> fixations;

    for(int i = 0; i < rowSize; i++)
    {
    	const int idx = idxs[i];
    	fixations.push_back(Fixation(idx, 0.0));
    	status = constraintPropagation(pp->cp, idx, 0.0, fixations, pp->unfixedVars, pp->cpColLb, pp->cpColUb,
    								   pp->constrBound, pp->unfixedVarsByRow);
    	if(status == CONFLICT)
    		break;
    }

    if(status == CONFLICT)
    	pp->transformed++;

	for(int j = 0; j < (int)fixations.size(); j++)
		unfixVariable(pp->cp, fixations[j].idxVar, pp->unfixedVars, pp->cpColLb, pp->cpColUb, pp->constrBound, pp->unfixedVarsByRow);
}

void one_element_constraint(Preprocess *pp, int idxRow)
{
    char sense = pp->rowSense[idxRow];
    int idx = (problem_row_idxs(pp->problem, idxRow))[0];
    double coef = (problem_row_coefs(pp->problem, idxRow))[0], rhs = problem_row_rhs(pp->problem, idxRow);

    /* constraints with one element are removed */
    pp->removeRow[idxRow] = 1;
    pp->colNConstraints[idx]--;
    pp->rowsRemoved++;

    if(coef <= -EPS)
    {
        if(sense == 'L') sense = 'G';
        else if(sense == 'G') sense = 'L';

        coef = -coef;
        rhs = -rhs;
    }

    const double newRhs = (rhs / coef);

    switch(sense)
    {
        case 'E':
            if(pp->colLb[idx] == pp->colUb[idx])
            {
                if(fabs(newRhs - pp->colLb[idx]) > EPS)
                {
                    fprintf(stderr, "Error: variable %s has already been fixed!\n", problem_var_name(pp->problem, idx));
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                fix_var(pp, idx, newRhs);
                pp->varsRemoved++;
            }
        break;

        case 'L':
            if(pp->colUb[idx] > newRhs + EPS)
            {
                if(problem_var_is_binary(pp->problem, idx))
                {
                    fix_var(pp, idx, 0.0);
                    pp->varsRemoved++;
                }
                else
                {
                    improve_upper_bound(pp, idx, newRhs);
                    pp->boundsImproved++;

                    if(pp->colLb[idx] == pp->colUb[idx])
                        pp->varsRemoved++;
                }
            }
        break;

        case 'G':
            if(pp->colLb[idx] + EPS < newRhs)
            {
                if(problem_var_is_binary(pp->problem, idx))
                {
                    fix_var(pp, idx, 1.0);
                    pp->varsRemoved++;
                }
                else
                {
                    improve_lower_bound(pp, idx, newRhs);
                    pp->boundsImproved++;

                    if(pp->colLb[idx] == pp->colUb[idx])
                        pp->varsRemoved++;
                }
            }
        break;

        default:
            fprintf(stderr, "Error: invalid type of constraint!\n");
            exit(EXIT_FAILURE);
    }
}

void handling_inequality_constraints(Preprocess *pp, int idxRow)
{
    int j, rowSize = problem_row_size(pp->problem, idxRow);
    char sense = pp->rowSense[idxRow];
    const double infty = problem_get_infinity(pp->problem);
    const int *idxs = problem_row_idxs(pp->problem, idxRow);

    if(pp->countInftyNegBounds[idxRow] == 0 && pp->lhsMin[idxRow] > pp->rhs[idxRow] + EPS)
    {
        fprintf(stderr, "Preprocessing says problem is infeasible!\n");
        exit(EXIT_SUCCESS);
    }

    if(pp->countInftyPosBounds[idxRow] == 0 && pp->lhsMax[idxRow] <= pp->rhs[idxRow])
    {
        pp->removeRow[idxRow] = 1; /* redundant row has been discovered. */
        pp->rowsRemoved++;

        for(j = 0; j < rowSize; j++)
        	pp->colNConstraints[idxs[j]]--;

        return;
    }

    for(j = 0; j < rowSize; j++)
    {
        const int idx = idxs[j];
        const double coef = pp->coefficients[idxRow][j];

        if(pp->colLb[idx] == pp->colUb[idx]) continue; /* ignoring fixed variables */
        
        if(problem_var_is_binary(pp->problem, idx))
        {
            double newLhsMin = (coef <= -EPS) ? (pp->lhsMin[idxRow] - (pp->colUb[idx] * coef))
                                              : (pp->lhsMin[idxRow] + (pp->colUb[idx] * coef));
            double newLhsMax = (coef <= -EPS) ? (pp->lhsMax[idxRow] + (pp->colUb[idx] * coef))
                                              : (pp->lhsMax[idxRow] - (pp->colUb[idx] * coef));

            if(pp->countInftyNegBounds[idxRow] == 0 && newLhsMin > pp->rhs[idxRow] + EPS)
            {
                double valueToFix = (coef <= -EPS) ? 1.0 : 0.0;
                fix_var(pp, idx, valueToFix);
                pp->varsRemoved++;
            }
            else if(pp->countInftyPosBounds[idxRow] == 0 && newLhsMax + EPS < pp->rhs[idxRow])
            {
                double delta = pp->rhs[idxRow] - newLhsMax;

                if(coef <= -EPS)
                {
                	double newCoef = coef + delta;
                	assert(newCoef > coef + EPS);
                	pp->coefficients[idxRow][j] = newCoef;
                	pp->coefsImproved++;
                	/* updating lhsMin e lhsMax */
                	pp->lhsMin[idxRow] += ((newCoef - coef) * pp->colUb[idx]);
                	pp->lhsMax[idxRow] += ((newCoef - coef) * pp->colLb[idx]);
                }
                else
                {
                	double newCoef = coef - delta;
                	assert(newCoef + EPS < coef);
                	pp->coefficients[idxRow][j] = newCoef;
                	pp->rhs[idxRow] = pp->rhs[idxRow] - delta;
                	pp->coefsImproved++;
                	/* updating lhsMin e lhsMax */
                	pp->lhsMin[idxRow] -= ((coef - newCoef) * pp->colLb[idx]);
                	pp->lhsMax[idxRow] -= ((coef - newCoef) * pp->colUb[idx]);
                }
            }
        }
        
        else
        {
            int countInftyNegBounds = (coef <= -EPS && pp->colUb[idx] >= infty) ? pp->countInftyNegBounds[idxRow] - 1
                                                                                : pp->countInftyNegBounds[idxRow];

            if(countInftyNegBounds == 0)
            {
                if(coef <= -EPS)
                {
                    double newLhsMin = (pp->colUb[idx] >= infty) ? pp->lhsMin[idxRow]
                     											 : (pp->lhsMin[idxRow] - (coef * pp->colUb[idx]));
                    double newBound = ((newLhsMin - pp->rhs[idxRow]) / -coef);

                    if(newBound > pp->colLb[idx] + EPS)
                    {
                        improve_lower_bound(pp, idx, newBound);
                        pp->boundsImproved++;

                        assert(pp->colLb[idx] <= pp->colUb[idx]);

                        if(pp->colLb[idx] == pp->colUb[idx])
                            pp->varsRemoved++;
                    }
                }
                else
                {
                    double newLhsMin = pp->lhsMin[idxRow] - (coef * pp->colLb[idx]);
                    double newBound = ((pp->rhs[idxRow] - newLhsMin) / coef);

                    if(pp->colUb[idx] > newBound + EPS)
                    {
                        improve_upper_bound(pp, idx, newBound);
                        pp->boundsImproved++;

                        if(pp->colLb[idx] == pp->colUb[idx])
                            pp->varsRemoved++;
                    }
                }
            }
    	}
    }
}

void handling_equality_constraints(Preprocess *pp, int idxRow)
{
	int j, rowSize = problem_row_size(pp->problem, idxRow);
    char sense = pp->rowSense[idxRow];
    const double infty = problem_get_infinity(pp->problem);
    const int *idxs = problem_row_idxs(pp->problem, idxRow);

	if(pp->countInftyNegBounds[idxRow] == 0 && pp->lhsMin[idxRow] > pp->rhs[idxRow] + EPS)
    {
        fprintf(stderr, "Preprocessing says problem is infeasible!\n");
        exit(EXIT_SUCCESS);
    }

    if(pp->countInftyPosBounds[idxRow] == 0 && -pp->lhsMax[idxRow] > -pp->rhs[idxRow] + EPS)
    {
        fprintf(stderr, "Preprocessing says problem is infeasible!\n");
        exit(EXIT_SUCCESS);
    }

    if( (pp->countInftyPosBounds[idxRow] == 0 && pp->lhsMax[idxRow] <= pp->rhs[idxRow]) &&
    	(pp->countInftyNegBounds[idxRow] == 0 && -pp->lhsMin[idxRow] <= -pp->rhs[idxRow]) )
    {
        pp->removeRow[idxRow] = 1; /* redundant row has been discovered. */
        pp->rowsRemoved++;

        for(j = 0; j < rowSize; j++)
        	pp->colNConstraints[idxs[j]]--;

        return;
    }

    for(j = 0; j < rowSize; j++)
    {
        const int idx = idxs[j];
        const double coef = pp->coefficients[idxRow][j];

        if(pp->colLb[idx] == pp->colUb[idx]) continue; /* ignoring fixed variables */

        if(pp->countInftyPosBounds[idxRow] == 0 && pp->lhsMax[idxRow] > pp->rhs[idxRow] + EPS)
        {
	        if(problem_var_is_binary(pp->problem, idx))
	        {
	            double newLhsMin = (coef <= -EPS) ? (pp->lhsMin[idxRow] - (pp->colUb[idx] * coef))
	                                              : (pp->lhsMin[idxRow] + (pp->colUb[idx] * coef));
	            double newLhsMax = (coef <= -EPS) ? (pp->lhsMax[idxRow] + (pp->colUb[idx] * coef))
	                                              : (pp->lhsMax[idxRow] - (pp->colUb[idx] * coef));

	            if(pp->countInftyNegBounds[idxRow] == 0 && newLhsMin > pp->rhs[idxRow] + EPS)
	            {
	                double valueToFix = (coef <= -EPS) ? 1.0 : 0.0;
	                fix_var(pp, idx, valueToFix);
	                pp->varsRemoved++;
	            }
	            else if(pp->countInftyPosBounds[idxRow] == 0 && newLhsMax + EPS < pp->rhs[idxRow])
	            {
	                double delta = pp->rhs[idxRow] - newLhsMax;

	                if(coef <= -EPS)
	                {
	                	double newCoef = coef + delta;
	                	assert(newCoef > coef + EPS);
	                	pp->coefficients[idxRow][j] = newCoef;
	                	pp->coefsImproved++;
	                	/* updating lhsMin e lhsMax */
	                	pp->lhsMin[idxRow] += ((newCoef - coef) * pp->colUb[idx]);
	                	pp->lhsMax[idxRow] += ((newCoef - coef) * pp->colLb[idx]);
	                }
	                else
	                {
	                	double newCoef = coef - delta;
	                	assert(newCoef + EPS < coef);
	                	pp->coefficients[idxRow][j] = newCoef;
	                	pp->rhs[idxRow] = pp->rhs[idxRow] - delta;
	                	pp->coefsImproved++;
	                	/* updating lhsMin e lhsMax */
	                	pp->lhsMin[idxRow] -= ((coef - newCoef) * pp->colLb[idx]);
	                	pp->lhsMax[idxRow] -= ((coef - newCoef) * pp->colUb[idx]);
	                }
	            }
	        }
	        
	        else
	        {
	            int countInftyNegBounds = (coef <= -EPS && pp->colUb[idx] >= infty) ? pp->countInftyNegBounds[idxRow] - 1
	                                                                                : pp->countInftyNegBounds[idxRow];

	            if(countInftyNegBounds == 0)
	            {
	                if(coef <= -EPS)
	                {
	                    double newLhsMin = (pp->colUb[idx] >= infty) ? pp->lhsMin[idxRow]
	                                                                 : (pp->lhsMin[idxRow] - (coef * pp->colUb[idx]));
	                    double newBound = ((newLhsMin - pp->rhs[idxRow]) / -coef);

	                    if(newBound > pp->colLb[idx] + EPS)
	                    {
	                        improve_lower_bound(pp, idx, newBound);
	                        pp->boundsImproved++;

	                        assert(pp->colLb[idx] <= pp->colUb[idx]);

	                        if(pp->colLb[idx] == pp->colUb[idx])
	                            pp->varsRemoved++;
	                    }
	                }
	                else
	                {
	                    double newLhsMin = pp->lhsMin[idxRow] - (coef * pp->colLb[idx]);
	                    double newBound = ((pp->rhs[idxRow] - newLhsMin) / coef);

	                    if(pp->colUb[idx] > newBound + EPS)
	                    {
	                        improve_upper_bound(pp, idx, newBound);
	                        pp->boundsImproved++;

	                        if(pp->colLb[idx] == pp->colUb[idx])
	                            pp->varsRemoved++;
	                    }
	                }
	            }
	    	}
	    }

	    if(pp->countInftyNegBounds[idxRow] == 0 && -pp->lhsMin[idxRow] > -pp->rhs[idxRow] + EPS)
        {
	        if(problem_var_is_binary(pp->problem, idx))
	        {
	            double newLhsMin = (-coef <= -EPS) ? (-pp->lhsMax[idxRow] + (pp->colUb[idx] * coef))
	                                               : (-pp->lhsMax[idxRow] - (pp->colUb[idx] * coef));
	            double newLhsMax = (-coef <= -EPS) ? (-pp->lhsMin[idxRow] - (pp->colUb[idx] * coef))
	                                               : (-pp->lhsMin[idxRow] + (pp->colUb[idx] * coef));

	            if(pp->countInftyPosBounds[idxRow] == 0 && newLhsMin > -pp->rhs[idxRow] + EPS)
	            {
	                double valueToFix = (-coef <= -EPS) ? 1.0 : 0.0;

	                if(pp->colLb[idx] == pp->colUb[idx])
                    {
                        if(valueToFix != pp->colLb[idx])
                        {
                            fprintf(stderr, "Error: variable %s has already been fixed!\n", problem_var_name(pp->problem, idx));
                            exit(EXIT_FAILURE);
                        }
                    }
                    else
                    {
                        fix_var(pp, idx, valueToFix);
                        pp->varsRemoved++;
                    }
	            }
	            else if(pp->countInftyNegBounds[idxRow] == 0 && newLhsMax + EPS < -pp->rhs[idxRow])
	            {
	                double delta = -pp->rhs[idxRow] - newLhsMax;

	                if(-coef <= -EPS)
	                {
	                	double newCoef = -coef + delta;
	                	assert(newCoef > -coef + EPS);
	                	pp->coefficients[idxRow][j] = -newCoef;
	                	pp->coefsImproved++;
	                	/* updating lhsMin e lhsMax */
	                	pp->lhsMin[idxRow] += ((-newCoef - coef) * pp->colUb[idx]);
	                	pp->lhsMax[idxRow] += ((-newCoef - coef) * pp->colLb[idx]);
	                }
	                else
	                {
	                	double newCoef = -coef - delta;
	                	assert(newCoef + EPS < -coef);
	                	pp->coefficients[idxRow][j] = -newCoef;
	                	pp->rhs[idxRow] = -pp->rhs[idxRow] + delta;
	                	pp->coefsImproved++;
	                	/* updating lhsMin e lhsMax */
	                	pp->lhsMin[idxRow] -= ((coef + newCoef) * pp->colLb[idx]);
	                	pp->lhsMax[idxRow] -= ((coef + newCoef) * pp->colUb[idx]);
	                }
	            }
	        }
	        
	        else
	        {
	            int countInftyNegBounds = (-coef <= -EPS && pp->colUb[idx] >= infty) ? pp->countInftyPosBounds[idxRow] - 1
	                                                                                 : pp->countInftyPosBounds[idxRow];

	            if(countInftyNegBounds == 0)
	            {
	                if(-coef <= -EPS)
	                {
	                    double newLhsMin = (pp->colUb[idx] >= infty) ? -pp->lhsMax[idxRow]
	                                                                 : (-pp->lhsMax[idxRow] + (coef * pp->colUb[idx]));
	                    double newBound = ((newLhsMin + pp->rhs[idxRow]) / coef);

	                    if(newBound > pp->colLb[idx] + EPS)
	                    {
	                        improve_lower_bound(pp, idx, newBound);
	                        pp->boundsImproved++;

	                        assert(pp->colLb[idx] <= pp->colUb[idx]);

	                        if(pp->colLb[idx] == pp->colUb[idx])
	                            pp->varsRemoved++;
	                    }
	                }
	                else
	                {
	                    double newLhsMin = -pp->lhsMax[idxRow] + (coef * pp->colLb[idx]);
	                    double newBound = ((pp->rhs[idxRow] + newLhsMin) / coef);

	                    if(pp->colUb[idx] > newBound + EPS)
	                    {
	                        improve_upper_bound(pp, idx, newBound);
	                        pp->boundsImproved++;

	                        if(pp->colLb[idx] == pp->colUb[idx])
	                            pp->varsRemoved++;
	                    }
	                }
	            }
	    	}
	    }
    }
}

void fix_var(Preprocess *pp, int idxVar, double valueToFix)
{
    assert(idxVar >= 0 && idxVar < problem_num_cols(pp->problem));
    assert(pp->colLb[idxVar] != pp->colUb[idxVar]);
    assert(valueToFix >= pp->colLb[idxVar] && valueToFix <= pp->colUb[idxVar]);

    int i;
    int nElements = problem_var_n_rows(pp->problem, idxVar);
    const int* idxs = problem_var_rows_idxs(pp->problem, idxVar);
    const double* coefs = problem_var_rows_coefs(pp->problem, idxVar);
    const double prevLb = pp->colLb[idxVar], prevUb = pp->colUb[idxVar];

    /* fixing bounds */
    pp->colLb[idxVar] = pp->colUb[idxVar] = valueToFix;
    
    for(i = 0; i < nElements; i++)
    {
        const int idxRow = idxs[i];
        const double mult = (pp->rowSense[idxRow] == 'G') ? -1.0 : 1.0; /* used to convert >= rows to <= */
        const double coefRow = mult * coefs[i];

        if(coefRow <= -EPS)
        {
            pp->lhsMin[idxRow] += (coefRow * (valueToFix - prevUb));
            pp->lhsMax[idxRow] += (coefRow * (valueToFix - prevLb));
        }
        else
        {
            pp->lhsMin[idxRow] += (coefRow * (valueToFix - prevLb));
            pp->lhsMax[idxRow] += (coefRow * (valueToFix - prevUb));
        }
    }
}

void improve_lower_bound(Preprocess *pp, int idxVar, double newBound)
{
    assert(idxVar >= 0 && idxVar < problem_num_cols(pp->problem));
    assert(newBound > pp->colLb[idxVar] + EPS);

    int i;
    int nElements = problem_var_n_rows(pp->problem, idxVar);
    const int* idxs = problem_var_rows_idxs(pp->problem, idxVar);
    const double* coefs = problem_var_rows_coefs(pp->problem, idxVar);
    const double prevLb = pp->colLb[idxVar];

    /* improving lower bound */
    char vType = problem_var_type(pp->problem, idxVar);
    assert(vType == INTEGER || vType == CONTINUOUS);
    if(vType == INTEGER)
    	newBound = ceil(newBound);
    pp->colLb[idxVar] = newBound;
    
    for(i = 0; i < nElements; i++)
    {
        const int idxRow = idxs[i];
        const double mult = (pp->rowSense[idxRow] == 'G') ? -1.0 : 1.0; /* used to convert >= rows to <= */
        const double coefRow = mult * coefs[i];

        if(coefRow <= -EPS)
            pp->lhsMax[idxRow] += (coefRow * (newBound - prevLb));
        else
            pp->lhsMin[idxRow] += (coefRow * (newBound - prevLb));
    }
}

void improve_upper_bound(Preprocess *pp, int idxVar, double newBound)
{
    assert(idxVar >= 0 && idxVar < problem_num_cols(pp->problem));
    assert(newBound + EPS < pp->colUb[idxVar]);

    int i;
    int nElements = problem_var_n_rows(pp->problem, idxVar);
    const int* idxs = problem_var_rows_idxs(pp->problem, idxVar);
    const double* coefs = problem_var_rows_coefs(pp->problem, idxVar);
    const double prevUb = pp->colUb[idxVar], infty = problem_get_infinity(pp->problem);

    /* improving upper bound */
    char vType = problem_var_type(pp->problem, idxVar);
    assert(vType == INTEGER || vType == CONTINUOUS);
    if(vType == INTEGER)
    	newBound = floor(newBound);
    pp->colUb[idxVar] = newBound;

    if(prevUb >= infty)
        for(i = 0; i < nElements; i++)
        {
            const int idxRow = idxs[i];
            const double mult = (pp->rowSense[idxRow] == 'G') ? -1.0 : 1.0; /* used to convert >= rows to <= */
            const double coefRow = mult * coefs[i];

            if(coefRow <= -EPS)
            {
                pp->countInftyNegBounds[idxRow]--;
                pp->lhsMin[idxRow] += (coefRow * newBound);
            }
            else
            {
                pp->countInftyPosBounds[idxRow]--;
                pp->lhsMax[idxRow] += (coefRow * newBound);
            }
        }

    else
        for(i = 0; i < nElements; i++)
        {
            const int idxRow = idxs[i];
            const double mult = (pp->rowSense[idxRow] == 'G') ? -1.0 : 1.0; /* used to convert >= rows to <= */
            const double coefRow = mult * coefs[i];

            if(coefRow <= -EPS)
                pp->lhsMin[idxRow] += (coefRow * (newBound - prevUb));
            else
                pp->lhsMax[idxRow] += (coefRow * (newBound - prevUb));
        }
}

Problem* preprocess_basic_preprocessing(Preprocess *pp)
{
    execute_basic_preprocessing(pp);

    int i = 0, j = 0;
    int initialNCols = problem_num_cols(pp->problem) - pp->varsRemoved + 1;
    int initialNRows = problem_num_rows(pp->problem) - pp->rowsRemoved + 1;
    int countCols = 0, countRows = 0;
    double sumFixedObj = 0.0;
    Problem *preProc = problem_create(initialNCols, initialNRows, problem_get_infinity(pp->problem));

    for(i = 0; i < problem_num_cols(pp->problem); i++)
    {
        if(pp->colLb[i] == pp->colUb[i])
        {
            sumFixedObj += (pp->colLb[i] * problem_var_obj_coef(pp->problem, i));
            continue;
        }

        assert(pp->colNConstraints[i] >= 0 && pp->colNConstraints[i] <= initialNRows);

        if(pp->colNConstraints[i] == 0 && fabs(problem_var_obj_coef(pp->problem, i)) < EPS)
        	continue;

        pp->nindexes[i] = countCols++;
        problem_var_set_lower_bound(preProc, pp->nindexes[i], pp->colLb[i]);
        problem_var_set_upper_bound(preProc, pp->nindexes[i], pp->colUb[i]);
        problem_var_set_obj_coef(preProc, pp->nindexes[i], problem_var_obj_coef(pp->problem, i));
        problem_var_set_type(preProc, pp->nindexes[i], problem_var_type(pp->problem, i));
        problem_var_set_name(preProc, pp->nindexes[i], problem_var_name(pp->problem, i));
    }

    for(i = 0; i < problem_num_rows(pp->problem); i++)
    {
        if(pp->removeRow[i]) continue;

        const int *idxs = problem_row_idxs(pp->problem, i), maxElements = problem_row_size(pp->problem, i);
        const char sense = pp->rowSense[i];
        double mult = ((sense == 'G') ? -1.0 : 1.0);

        int newIdxs[maxElements], newNElements = 0;
        double newCoefs[maxElements], newRhs = mult * pp->rhs[i];
        double discount = 0.0;

        for(j = 0; j < maxElements; j++)
        {
            const int idx = idxs[j];
            const double coef = pp->coefficients[i][j] * mult;

            if(pp->colLb[idx] == pp->colUb[idx])
            {
                discount += (coef * pp->colLb[idx]);
                continue;
            }

            if(fabs(coef) < EPS) continue;

            newIdxs[newNElements] = pp->nindexes[idx];
            newCoefs[newNElements] = coef;
            newNElements++;
        }

        if(newNElements > 0)
        {
            problem_set_row(preProc, countRows, newIdxs, newCoefs, newNElements, newRhs-discount, sense);
            problem_row_set_name(preProc, countRows, problem_row_name(pp->problem, i));
            countRows++;
        }
    }

    if(fabs(sumFixedObj) > EPS)
    {
        problem_var_set_lower_bound(preProc, countCols, 0.0);
        problem_var_set_upper_bound(preProc, countCols, 1.0);
        problem_var_set_obj_coef(preProc, countCols, sumFixedObj);
        problem_var_set_type(preProc, countCols, BINARY);
        problem_var_set_name(preProc, countCols, "sum_fixed_obj");

        int newIdxs[] = { countCols };
        double newCoefs[] = { 1.0 };
        problem_set_row(preProc, countRows, newIdxs, newCoefs, 1, 1.0, 'E');
        problem_row_set_name(preProc, countRows, "fix_sum_obj_coef");

        problem_set_num_cols(preProc, countCols + 1);
        problem_set_num_rows(preProc, countRows + 1);
    }
    else
    {
        problem_set_num_cols(preProc, countCols);
        problem_set_num_rows(preProc, countRows);
    }

    problem_update_matrices_by_col(preProc);

    //printf("%d %d %d %d", pp->rowsRemoved, pp->varsRemoved, pp->coefsImproved, pp->boundsImproved);
    printf("%d %d", pp->numGub, pp->transformed);

    return preProc;
}
