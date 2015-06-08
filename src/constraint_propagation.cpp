#include "constraint_propagation.h"
#include <limits>
#include <cfloat>
#include <cmath>
#include <vector>
#include <set>
#include <queue>

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

typedef struct Fixation
{
	int idxVar;
	double valueToFix;
} Fixation;

inline bool operator<(const Fixation &f1, const Fixation &f2)
{
	if(f1.idxVar != f2.idxVar)
		return f1.idxVar < f2.idxVar;
	return f1.valueToFix < f2.valueToFix;
}

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
	int binaryVars, unfixedVars; /* number of binary variables and number of unfixed variables */
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

char constraintPropagation(CPropagation *cp, int var, double value, vector<Fixation> &fixations);

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

	cp->maxDepth = 2;
    cp->numCols = solver->getNumCols();
	cp->numRows = (int)cp->matrixByRow.size();
    cp->constrBound.resize(cp->numRows);
    cp->unfixedVarsByRow.resize(cp->numRows);
    cp->colLb = new double[cp->numCols];
    cp->colUb = new double[cp->numCols];
    cp->varIsBinary = new char[cp->numCols];
    cp->binaryVars = 0;

    for(int i = 0; i < cp->numCols; i++)
    {
    	cp->colLb[i] = colLb[i];
    	cp->colUb[i] = colUb[i];

        if(cp->colLb[i] == 0.0 && cp->colUb[i] == 1.0)
        {
            cp->varIsBinary[i] = 1;
            cp->binaryVars++;
        }
        else cp->varIsBinary[i] = 0;
    }
    cp->unfixedVars = cp->binaryVars;

    /* calculating the initial lower bound for each constraint */
    for(int i = 0; i < cp->numRows; i++)
    {
        cp->constrBound[i] = cp->rhs[i];
        cp->unfixedVarsByRow[i] = (int)cp->matrixByRow[i].size();

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

char constraintPropagation(CPropagation *cp, int var, double value, vector<Fixation> &fixations)
{
	assert(value == 0.0 || value == 1.0);
	assert(var >= 0 && var < cp->numCols);
	assert(cp->colLb[var] != cp->colUb[var]);

	char status = NOIMPLICATION;
	vector<int> C;
	Fixation tmp;
	vector<char> rowIsInQueue(cp->numRows, 0);

	C.reserve(cp->numRows);

	fixations.clear(); //cleaning old fixations
	fixations.reserve(cp->numCols);

	fixVariable(cp, var, value);

	char varBounds[cp->numCols];
	for(int i = 0; i < cp->numCols; i++) varBounds[i] = UNFIXED;

	for(int i = 0; i < (int)cp->matrixByCol[var].size(); i++)
	{
		const pair<int, double> constraint = cp->matrixByCol[var][i];
		const int idxRow = constraint.first;
		const double coef = constraint.second;

		if((cp->unfixedVarsByRow[idxRow] > 0) && (value == 1.0 && coef > 0.0) || (value == 0.0 && coef < 0.0))
		{
			C.push_back(idxRow);
			rowIsInQueue[idxRow] = 1;
		}
	}

	if(C.empty())
		return status;

	do {
		char newImplication = 0;
		vector<Fixation> currFixations;
		for(int i = 0; i < (int)C.size(); i++)
		{
			const int row = C[i];
			rowIsInQueue[row] = 0;

			for(int j = 0; j < (int)cp->matrixByRow[row].size(); j++)
			{
				const pair<int, double> var = cp->matrixByRow[row][j];
				const int idx = var.first;
	            const double coef = var.second;

	            if(cp->colLb[idx] == cp->colUb[idx]) continue;

	            char eval = evaluateBound(cp, var, row);
            	
	            if(eval == DEACTIVATE && varBounds[idx] != DEACTIVATE)
	            {
	            	if(varBounds[idx] == ACTIVATE)
	            		return CONFLICT;
	            	tmp.idxVar = idx;
	            	tmp.valueToFix = 0.0;
	            	varBounds[idx] = DEACTIVATE;
	            	currFixations.push_back(tmp);
	            	newImplication = 1;
	            	status = FIXATION;
	            }
	            else if(eval == ACTIVATE && varBounds[idx] != ACTIVATE)
	            {
	            	if(varBounds[idx] == DEACTIVATE)
	            		return CONFLICT;
	            	tmp.idxVar = idx;
	            	tmp.valueToFix = 1.0;
	            	varBounds[idx] = ACTIVATE;
	            	currFixations.push_back(tmp);
	            	newImplication = 1;
	            	status = FIXATION;
	            }
	            else if(eval == CONFLICT)
	            	return CONFLICT;
			}
		}

		if(!newImplication)
			return status;

		C.clear();

		for(int i = 0; i < (int)currFixations.size(); i++)
		{
			const int idxVar = currFixations[i].idxVar;
			const double valueToFix = currFixations[i].valueToFix;

			assert(cp->colLb[idxVar] != cp->colUb[idxVar]);

			cp->unfixedVars--;
            cp->colLb[idxVar] = valueToFix;
            cp->colUb[idxVar] = valueToFix;
            fixations.push_back(currFixations[i]);

            for(int j = 0; j < (int)cp->matrixByCol[idxVar].size(); j++)
            {
                const pair<int, double> constraint = cp->matrixByCol[idxVar][j];
                const int idxRow = constraint.first;
                const double coef = constraint.second;

                cp->unfixedVarsByRow[idxRow]--;

                if(valueToFix == 1.0 && coef > 0.0)
                {
                    cp->constrBound[idxRow] -= coef;
                    if(cp->unfixedVarsByRow[idxRow] > 0 && !rowIsInQueue[idxRow])
                    {
                    	C.push_back(idxRow);
                    	rowIsInQueue[idxRow] = 1;
                    }
                }
                else if(valueToFix == 0.0 && coef < 0.0)
                {
                    cp->constrBound[idxRow] += coef;
                    if(cp->unfixedVarsByRow[idxRow] > 0 && !rowIsInQueue[idxRow])
                    {
                    	C.push_back(idxRow);
                    	rowIsInQueue[idxRow] = 1;
                    }
                }
                
                if(cp->unfixedVarsByRow[idxRow] == 0 && cp->constrBound[idxRow] < 0.0)
                	status = CONFLICT; //nao posso retornar direto pq tenho q atualizar os bounds das retricoes
            }
		}
		if(status == CONFLICT)
			return CONFLICT;
	} while(!C.empty());

	return status;
}

/*char constraintPropagation(CPropagation *cp, int var, double value, vector<Fixation> &fixations)
{
	assert(value == 0.0 || value == 1.0);
	assert(var >= 0 && var < cp->numCols);
	assert(cp->colLb[var] != cp->colUb[var]);
	fixations.clear(); fixations.reserve(cp->numCols); //cleaning old fixations and reserving memory

	char status = NOIMPLICATION;
	queue<int> C;
	Fixation tmp;
	vector<char> rowIsInQueue(cp->numRows, 0);

	fixVariable(cp, var, value);

	for(int i = 0; i < (int)cp->matrixByCol[var].size(); i++)
	{
		const pair<int, double> constraint = cp->matrixByCol[var][i];
		const int idxRow = constraint.first;
		const double coef = constraint.second;

		if((cp->unfixedVarsByRow[idxRow] > 0) && (value == 1.0 && coef > 0.0) || (value == 0.0 && coef < 0.0))
		{
			C.push(idxRow);
			rowIsInQueue[idxRow] = 1;
		}
	}

	while(!C.empty())
	{
		const int idxRow = C.front();
		C.pop();
		rowIsInQueue[idxRow] = 0;

		for(int i = 0; i < (int)cp->matrixByRow[idxRow].size(); i++)
		{
			const pair<int, double> var = cp->matrixByRow[idxRow][i];
			const int idxVar = var.first;
			const double coef = var.second;

			if(cp->colLb[idxVar] == cp->colUb[idxVar])
				continue;

			char eval = evaluateBound(cp, var, idxRow);

            if(eval == ACTIVATE || eval == DEACTIVATE)
            {
            	assert(cp->colLb[idxVar] != cp->colUb[idxVar]);
            	tmp.idxVar = idxVar;
            	tmp.valueToFix = eval;
            	fixations.push_back(tmp);
            	status = FIXATION;
				cp->unfixedVars--;
	            cp->colLb[idxVar] = cp->colUb[idxVar] = eval;

	            for(int j = 0; j < (int)cp->matrixByCol[idxVar].size(); j++)
	            {
	                const pair<int, double> constraint = cp->matrixByCol[idxVar][j];
	                const int idxRow2 = constraint.first;
	                const double coef2 = constraint.second;

	                cp->unfixedVarsByRow[idxRow2]--;

	                if(eval == ACTIVATE && coef2 > 0.0)
	                {
	                    cp->constrBound[idxRow2] -= coef2;
	                    if(cp->unfixedVarsByRow[idxRow2] > 0 && !rowIsInQueue[idxRow2])
	                    {
	                    	C.push(idxRow2);
	                    	rowIsInQueue[idxRow2] = 1;
	                    }
	                }
	                else if(eval == DEACTIVATE && coef2 < 0.0)
	                {
	                    cp->constrBound[idxRow2] += coef2;
	                    if(cp->unfixedVarsByRow[idxRow2] > 0 && !rowIsInQueue[idxRow2])
	                    {
	                    	C.push(idxRow2);
	                    	rowIsInQueue[idxRow2] = 1;
	                    }
	                }
                    
                    if(cp->unfixedVarsByRow[idxRow2] == 0 && cp->constrBound[idxRow2] < 0.0)
                    	status = CONFLICT; //nao posso retornar direto pq tenho q atualizar os bounds das retricoes
	            }
	            if(status == CONFLICT)
	            	return CONFLICT;
        	}
            else if(eval == CONFLICT)
            	return CONFLICT;
		}
	}
	return status;
}*/

void fixVariable(CPropagation *cp, int idx, double value)
{
	assert(value == 0.0 || value == 1.0);
	assert(idx >= 0 && idx < cp->numCols);

	cp->unfixedVars--;
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
    assert(cp->colLb[idx] == cp->colUb[idx]);

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

    cp->unfixedVars++;
    cp->colLb[idx] = 0.0;
    cp->colUb[idx] = 1.0;
}

void backtrack(CPropagation *cp, int idxVar, double valueToFix, int currDepth, vector<int> confs)
{
	if(currDepth > cp->maxDepth) //explores untill maxDepth
		return;

	char status;
	vector<Fixation> fixations;
	status = constraintPropagation(cp, idxVar, valueToFix, fixations);
	if(valueToFix == 1.0)
		confs.push_back(idxVar);
	else confs.push_back(idxVar+cp->numCols);

	if(status == CONFLICT)
	{
		cp->nogoods.push_back(confs);
		unfixVariable(cp, idxVar);
		for(int i = 0; i < (int)fixations.size(); i++)
			unfixVariable(cp, fixations[i].idxVar);
		return;
	}

	if(cp->unfixedVars == 0) //all variables were fixed
	{
		unfixVariable(cp, idxVar);
		for(int i = 0; i < (int)fixations.size(); i++)
			unfixVariable(cp, fixations[i].idxVar);
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
	for(int i = 0; i < (int)fixations.size(); i++)
			unfixVariable(cp, fixations[i].idxVar);
}

void discoverConflicts(CPropagation *cp)
{
    char status;
    vector<Fixation> fixations1, fixations2;

    for(int i = 0; i < cp->numCols; i++)
    {
        if(!cp->varIsBinary[i]) continue;
        assert(cp->unfixedVars == cp->binaryVars);
        status = constraintPropagation(cp, i, 0.0, fixations1);
        if(status == CONFLICT)
            cp->varsToRemove.push_back(i+cp->numCols);
        /*else
        {
            for(int j = i + 1; j < cp->numCols; j++)
            {
                if(!cp->varIsBinary[j] || (cp->colLb[j] == cp->colUb[j]))
                    continue;

                status = constraintPropagation(cp, j, 0.0, fixations2);
                if(status == CONFLICT)
                    cp->newConflicts.push_back(pair<int, int>(i+cp->numCols, j+cp->numCols));
                unfixVariable(cp, j);
                for(int k = 0; k < (int)fixations2.size(); k++)
					unfixVariable(cp, fixations2[k].idxVar);

                status = constraintPropagation(cp, j, 1.0, fixations2);
                if(status == CONFLICT)
                    cp->newConflicts.push_back(pair<int, int>(i+cp->numCols, j));
                unfixVariable(cp, j);
                for(int k = 0; k < (int)fixations2.size(); k++)
					unfixVariable(cp, fixations2[k].idxVar);
            }
        }*/
        unfixVariable(cp, i);
        for(int j = 0; j < (int)fixations1.size(); j++)
			unfixVariable(cp, fixations1[j].idxVar);
        assert(cp->unfixedVars == cp->binaryVars);

        status = constraintPropagation(cp, i, 1.0, fixations1);
        if(status == CONFLICT)
            cp->varsToRemove.push_back(i);
        /*else
        {
            for(int j = i + 1; j < cp->numCols; j++)
            {
                if(!cp->varIsBinary[j] || (cp->colLb[j] == cp->colUb[j]))
                    continue;

                status = constraintPropagation(cp, j, 0.0, fixations2);
                if(status == CONFLICT)
                    cp->newConflicts.push_back(pair<int, int>(i, j+cp->numCols));
                unfixVariable(cp, j);
                for(int k = 0; k < (int)fixations2.size(); k++)
					unfixVariable(cp, fixations2[k].idxVar);

                status = constraintPropagation(cp, j, 1.0, fixations2);
                if(status == CONFLICT)
                    cp->newConflicts.push_back(pair<int, int>(i, j));
                unfixVariable(cp, j);
                for(int k = 0; k < (int)fixations2.size(); k++)
					unfixVariable(cp, fixations2[k].idxVar);
            }
        }*/
        unfixVariable(cp, i);
        for(int j = 0; j < (int)fixations1.size(); j++)
			unfixVariable(cp, fixations1[j].idxVar);
        assert(cp->unfixedVars == cp->binaryVars);
    }
}

void test(CPropagation *cp, OsiSolverInterface *solver)
{
    clock_t start = clock();
    discoverConflicts(cp);
    printf("Time elapsed: %.2lf seconds.\n", (double(clock() - start))/((double)CLOCKS_PER_SEC));
    printf("Variables to remove: ");
    for(int i = 0; i < (int)cp->varsToRemove.size(); i++)
        printf("%d ", cp->varsToRemove[i]);
    printf("\n\n");

	/*for(int i = 0; i < (int)cp->varsToRemove.size(); i++)
	{
		solver->setColBounds(cp->varsToRemove[i], 1.0, 1.0);
		solver->initialSolve();
		assert(!solver->isProvenOptimal());
	    solver->setColBounds(cp->varsToRemove[i], 0.0, 1.0);
	}*/
}