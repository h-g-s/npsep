#include "constraint_propagation.h"
#include <limits>
#include <cfloat>
#include <vector>
#include <cmath>

extern "C"
{
	#include "memory.h"
}

using namespace std;

#define DEACTIVATE 0
#define ACTIVATE 1
#define UNFIXED 2

#define CONFLICT 3
#define FIXATION 4
#define NOIMPLICATION 5

/* stores constraints which have only binary variables */
struct _CPropagation
{
	int numCols, numRows;
	vector<vector<pair<int, double> > > matrixByRow; /* matrix of indexes and coefficients by row */
    vector<vector<pair<int, double> > > matrixByCol; /* matrix of indexes and coefficients by column */
    double *colLb, *colUb; /* bounds of variables */
    char *varIsBinary; /* for each variable, stores 1 if variable is binary and 0 otherwise */
	vector<double> rhs; /* right-hand side for each constraint */
	vector<double> constrBound; /* lower bound for LHS of each constraint */
	vector<int> unfixedVarsByRow; /* number of unfixed variables for each constraint */
	int fixedVars; /* number of fixed variables */
	int maxDepth;

    /* nogoods discovered in backtracking */
	vector<vector<int> > nogoods;

	/* variables to remove of the problem - discovered in function discoverConflicts
	   index of variable >= numCols indicates that the variable can not be assigned with 0
	   index of variable < numCols indicates that the variable can not be assigned with 1 */
	vector<int> varsToRemove;

    /* conflicts between two variables - discovered in function discoverConflicts */
	vector<pair<int, int> > newConflicts;
};

void fillMatrices(const OsiSolverInterface *solver, CPropagation *cp);

char evaluateBound(const CPropagation *cp, const pair<int, double> &var, int constraint);

char constraintPropagation(CPropagation *cp, int idxVar, char *varBounds);

void fixVariable(CPropagation *cp, int idx, double value);

void unfixVariable(CPropagation *cp, int idx);

void backtrack(CPropagation *cp, int idxVar, double valueToFix, int currDepth, vector<int> confs);

CPropagation *cpropagation_create(const OsiSolverInterface *solver)
{
	CPropagation *cp = new CPropagation;
	const double *colLb = solver->getColLower();
    const double *colUb = solver->getColUpper();

    fillMatrices(solver, cp);
	assert(cp->matrixByRow.size() == cp->rhs.size());

	cp->fixedVars = 0;
	cp->maxDepth = 2;
    cp->numCols = solver->getNumCols();
	cp->numRows = (int)cp->matrixByRow.size();
    cp->constrBound.resize(cp->numRows);
    cp->unfixedVarsByRow.resize(cp->numRows);
    cp->colLb = new double[cp->numCols];
    cp->colUb = new double[cp->numCols];
    cp->varIsBinary = new char[cp->numCols];

    for(int i = 0; i < cp->numCols; i++)
    {
    	cp->colLb[i] = colLb[i];
    	cp->colUb[i] = colUb[i];

        if(cp->colLb[i] == 0.0 && cp->colUb[i] == 1.0)
            cp->varIsBinary[i] = 1;
        else cp->varIsBinary[i] = 0;
    }

    /* calculating the initial lower bound for each constraint */
    for(int i = 0; i < cp->numRows; i++)
    {
        cp->constrBound[i] = cp->rhs[i];
        cp->unfixedVarsByRow[i] = 0;

        for(int j = 0; j < (int)cp->matrixByRow[i].size(); j++)
        {
            const int idx = cp->matrixByRow[i][j].first;
            const double coef = cp->matrixByRow[i][j].second;

            if(coef > 0.0) cp->constrBound[i] -= (coef * cp->colLb[idx]);
            else cp->constrBound[i] -= (coef * cp->colUb[idx]);
        }
    }

	return cp;
}

void cpropagation_free(CPropagation *cp)
{
	delete[] cp->colLb;
	delete[] cp->colUb;
    delete[] cp->varIsBinary;
	delete cp;
}

void fillMatrices(const OsiSolverInterface *solver, CPropagation *cp)
{
    const CoinPackedMatrix *M = solver->getMatrixByRow();
    const double *rhs = solver->getRightHandSide();
    const char *sense = solver->getRowSense();
    const char *ctype = solver->getColType();

    for(int idxRow = 0; idxRow < solver->getNumRows(); idxRow++)
    {
        if(sense[idxRow] == 'R') /* ignoring ranged constarints */
            continue;

        const CoinShallowPackedVector &row = M->getVector(idxRow);
        const int nElements = row.getNumElements();
        const int *idxs = row.getIndices();
        const double *coefs = row.getElements();
        double mult = (sense[idxRow] == 'G') ? -1.0 : 1.0;
        vector<pair<int, double> > constraint(nElements);
        bool allBinaries = true;

        for(int i = 0; i < nElements; i++)
        {
            constraint[i] = pair<int, double> (idxs[i], mult * coefs[i]);
            if(ctype[idxs[i]] != 1)
            {
                allBinaries = false;
                break;
            }
        }

        if(allBinaries)
        {
            cp->matrixByRow.push_back(constraint);
            cp->rhs.push_back(rhs[idxRow] * mult);
            if(sense[idxRow] == 'E')
            {
                for(int j = 0; j < nElements; j++)
                    constraint[j].second = -1.0 * constraint[j].second;
                cp->matrixByRow.push_back(constraint);
                cp->rhs.push_back(-1.0 * rhs[idxRow]);
            }
        }
    }

    cp->matrixByCol.resize(solver->getNumCols());
    for(int i = 0; i < (int)cp->matrixByRow.size(); i++)
        for(int j = 0; j < (int)cp->matrixByRow[i].size(); j++)
        {
            const int var = cp->matrixByRow[i][j].first;
            const double coef = cp->matrixByRow[i][j].second;
            cp->matrixByCol[var].push_back(pair<int, double>(i, coef));
        }
}

char evaluateBound(const CPropagation *cp, const pair<int, double> &var, int constraint)
{
    const double bound = cp->constrBound[constraint];

    if(var.second > 0.0) /* positive coefficient */
    {
        if(bound < 0.0)
            return CONFLICT;
        else if(bound - var.second < 0.0)
            return DEACTIVATE;
        else return NOIMPLICATION;
    }

    else /* negative coefficient */
    {
        if(bound < 0.0)
            return CONFLICT;
        else if(bound + var.second < 0.0)
            return ACTIVATE;
        else return NOIMPLICATION;
    }

    return NOIMPLICATION; /* just to avoid errors */
}

char constraintPropagation(CPropagation *cp, int idxVar, char *varBounds)
{
	char status = NOIMPLICATION;
	vector<int> C;

	for(int i = 0; i < cp->numCols; i++)
		varBounds[i] = UNFIXED;

	for(int i = 0; i < (int)cp->matrixByCol[idxVar].size(); i++)
	{
		const pair<int, double> constraint = cp->matrixByCol[idxVar][i];
		const int idx = constraint.first;
		C.push_back(idx);
	}

	do {
		char newImplication = 0;
		for(int i = 0; i < (int)C.size(); i++)
		{
			const int row = C[i];

			for(int j = 0; j < (int)cp->matrixByRow[row].size(); j++)
			{
				const pair<int, double> var = cp->matrixByRow[row][j];
				const int idx = var.first;
	            const double coef = var.second;

	            if(cp->colLb[idx] == cp->colUb[idx]) continue;

	            char eval = evaluateBound(cp, var, row);

	            if(eval == DEACTIVATE)
	            {
	            	if(varBounds[idx] == ACTIVATE)
	            		return CONFLICT;
	            	newImplication = 1;
	            	varBounds[idx] = DEACTIVATE;
	            	status = FIXATION;
	            }
	            else if(eval == ACTIVATE)
	            {
	            	if(varBounds[idx] == DEACTIVATE)
		            	return CONFLICT;
	            	newImplication = 1;
	            	varBounds[idx] = ACTIVATE;
	            	status = FIXATION;
	            }
	            else if(eval == CONFLICT)
	            	return CONFLICT;
			}
		}

		C.clear();

		if(!newImplication)
			return status;

		/* fixing values */
		for(int i = 0; i < cp->numCols; i++)
		{
            if(cp->colLb[i] == cp->colUb[i]) continue;

            if(varBounds[i] == ACTIVATE || varBounds[i] == DEACTIVATE)
            {
            	cp->fixedVars++;
                cp->colLb[i] = varBounds[i];
                cp->colUb[i] = varBounds[i];

                for(int j = 0; j < (int)cp->matrixByCol[i].size(); j++)
                {
                    const pair<int, double> constraint = cp->matrixByCol[i][j];
                    const int idxRow = constraint.first;
                    const double coef = constraint.second;

                    cp->unfixedVarsByRow[idxRow]--;

                    if(cp->unfixedVarsByRow[idxRow] > 0)
                        C.push_back(idxRow);

                    if(varBounds[i] == ACTIVATE && coef > 0.0)
                        cp->constrBound[idxRow] -= coef;
                    else if(varBounds[i] == DEACTIVATE && coef < 0.0)
                        cp->constrBound[idxRow] += coef;
                }
            }
		}
	} while(!C.empty());

	return status;
}

void fixVariable(CPropagation *cp, int idx, double value)
{
	assert(value == 0.0 || value == 1.0);
	assert(idx >= 0 && idx < cp->numCols);

	cp->fixedVars++;
	cp->colLb[idx] = value;
	cp->colUb[idx] = value;

	for(int j = 0; j < (int)cp->matrixByCol[idx].size(); j++)
	{
		const pair<int, double> constraint = cp->matrixByCol[idx][j];
		const int idxRow = constraint.first;
		const double coef = constraint.second;

		cp->unfixedVarsByRow[idxRow]--;

		if(value == 1.0 && coef > 0.0)
			cp->constrBound[idxRow] -= coef;
		else if(value == 0.0 && coef < 0.0)
			cp->constrBound[idxRow] += coef;
	}
}

void unfixVariable(CPropagation *cp, int idx)
{
    if(cp->colLb[idx] == 0.0 && cp->colUb[idx] == 1.0)
        return;

    for(int j = 0; j < (int)cp->matrixByCol[idx].size(); j++)
    {
        const pair<int, double> constraint = cp->matrixByCol[idx][j];
        const int idxRow = constraint.first;
        const double coef = constraint.second;

        cp->unfixedVarsByRow[idxRow]++;

        if(cp->colLb[idx] == 1.0 && coef > 0.0) //fixed in 1
            cp->constrBound[idxRow] += coef;
        else if(cp->colUb[idx] == 0.0 && coef < 0.0) //fixex in 0
            cp->constrBound[idxRow] -= coef;
    }

    cp->fixedVars--;
    cp->colLb[idx] = 0.0;
    cp->colUb[idx] = 1.0;
}

void backtrack(CPropagation *cp, int idxVar, double valueToFix, int currDepth, vector<int> confs)
{
	if(currDepth > cp->maxDepth) //explores untill maxDepth
		return;

	fixVariable(cp, idxVar, valueToFix);
	if(valueToFix == 1.0)
		confs.push_back(idxVar);
	else confs.push_back(idxVar+cp->numCols);

	char status, varBounds[cp->numCols];
	status = constraintPropagation(cp, idxVar, varBounds);

	if(status == CONFLICT)
	{
		cp->nogoods.push_back(confs);
		unfixVariable(cp, idxVar);
		for(int i = 0; i < cp->numCols; i++)
			if(varBounds[i] == ACTIVATE || varBounds[i] == DEACTIVATE)
				unfixVariable(cp, i);
		return;
	}

	if(cp->fixedVars == cp->numCols) //all variables were fixed
	{
		unfixVariable(cp, idxVar);
		for(int i = 0; i < cp->numCols; i++)
			if(varBounds[i] == ACTIVATE || varBounds[i] == DEACTIVATE)
				unfixVariable(cp, i);
		return;
	}

	for(int i = idxVar+1; i < cp->numCols; i++)
	{
		if(cp->colLb[i] != cp->colUb[i])
		{
			backtrack(cp, i, 0.0, currDepth+1, confs);
			unfixVariable(cp, i);
			backtrack(cp, i, 1.0, currDepth+1, confs);
			unfixVariable(cp, i);
		}
	}

	unfixVariable(cp, idxVar);
	for(int i = 0; i < cp->numCols; i++)
		if(varBounds[i] == ACTIVATE || varBounds[i] == DEACTIVATE)
			unfixVariable(cp, i);
}

void discoverConflicts(CPropagation *cp)
{
    char status, varBounds1[cp->numCols], varBounds2[cp->numCols];

    for(int i = 0; i < cp->numCols; i++)
    {
        if(!cp->varIsBinary[i]) continue;
        assert(cp->fixedVars == 0);
        fixVariable(cp, i, 0.0);
        status = constraintPropagation(cp, i, varBounds1);
        if(status == CONFLICT)
            cp->varsToRemove.push_back(i+cp->numCols);
        else
        {
            for(int j = i + 1; j < cp->numCols; j++)
            {
                if(!cp->varIsBinary[j] || (cp->colLb[j] == cp->colUb[j]))
                    continue;

                fixVariable(cp, j, 0.0);
                status = constraintPropagation(cp, j, varBounds2);
                if(status == CONFLICT)
                    cp->newConflicts.push_back(pair<int, int>(i+cp->numCols, j+cp->numCols));
                unfixVariable(cp, j);
                for(int k = 0; k < cp->numCols; k++)
                    if(varBounds2[k] == ACTIVATE || varBounds2[k] == DEACTIVATE)
                        unfixVariable(cp, k);

                fixVariable(cp, j, 1.0);
                status = constraintPropagation(cp, j, varBounds2);
                if(status == CONFLICT)
                    cp->newConflicts.push_back(pair<int, int>(i+cp->numCols, j));
                unfixVariable(cp, j);
                for(int k = 0; k < cp->numCols; k++)
                    if(varBounds2[k] == ACTIVATE || varBounds2[k] == DEACTIVATE)
                        unfixVariable(cp, k);
            }
        }
        unfixVariable(cp, i);
        for(int j = 0; j < cp->numCols; j++)
            if(varBounds1[j] == ACTIVATE || varBounds1[j] == DEACTIVATE)
                unfixVariable(cp, j);
        assert(cp->fixedVars == 0);

        fixVariable(cp, i, 1.0);
        status = constraintPropagation(cp, i, varBounds1);
        if(status == CONFLICT)
            cp->varsToRemove.push_back(i);
        else
        {
            for(int j = i + 1; j < cp->numCols; j++)
            {
                if(!cp->varIsBinary[j] || (cp->colLb[j] == cp->colUb[j]))
                    continue;

                fixVariable(cp, j, 0.0);
                status = constraintPropagation(cp, j, varBounds2);
                if(status == CONFLICT)
                    cp->newConflicts.push_back(pair<int, int>(i, j+cp->numCols));
                unfixVariable(cp, j);
                for(int k = 0; k < cp->numCols; k++)
                    if(varBounds2[k] == ACTIVATE || varBounds2[k] == DEACTIVATE)
                        unfixVariable(cp, k);

                fixVariable(cp, j, 1.0);
                status = constraintPropagation(cp, j, varBounds2);
                if(status == CONFLICT)
                    cp->newConflicts.push_back(pair<int, int>(i, j));
                unfixVariable(cp, j);
                for(int k = 0; k < cp->numCols; k++)
                    if(varBounds2[k] == ACTIVATE || varBounds2[k] == DEACTIVATE)
                        unfixVariable(cp, k);
            }
        }
        unfixVariable(cp, i);
        for(int j = 0; j < cp->numCols; j++)
            if(varBounds1[j] == ACTIVATE || varBounds1[j] == DEACTIVATE)
                unfixVariable(cp, j);
        assert(cp->fixedVars == 0);
    }
}

void test(CPropagation *cp)
{
	/*clock_t start = clock();
	cp->maxDepth = 1;
	for(int i = 0; i < cp->numCols; i++)
	{
		vector<int> confs;
		assert(cp->fixedVars == 0);
		backtrack(cp, i, 0.0, 1, confs);
		confs.clear();
		assert(cp->fixedVars == 0);
		backtrack(cp, i, 1.0, 1, confs);
		assert(cp->fixedVars == 0);
	}
	printf("Time elapsed in backtracking: %.2lf seconds.\n", (double(clock() - start))/((double)CLOCKS_PER_SEC));
	for(int i = 0; i < (int)cp->nogoods.size(); i++)
	{
		for(int j = 0; j < (int)cp->nogoods[i].size(); j++)
			printf("%d ", cp->nogoods[i][j]);
		printf("\n");
	}*/

    clock_t start = clock();
    discoverConflicts(cp);
    printf("Time elapsed: %.2lf seconds.\n", (double(clock() - start))/((double)CLOCKS_PER_SEC));

    printf("Variables to remove: ");
    for(int i = 0; i < (int)cp->varsToRemove.size(); i++)
        printf("%d ", cp->varsToRemove[i]);
    printf("\n\n");

    printf("New pairs of conflicts:\n");
    for(int i = 0; i < (int)cp->newConflicts.size(); i++)
        printf("(%d, %d) ", cp->newConflicts[i].first, cp->newConflicts[i].second);
    printf("\n");
}