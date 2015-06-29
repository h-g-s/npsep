#include "constraint_propagation.h"
#include <limits>
#include <cfloat>
#include <cmath>
#include <vector>
#include <set>
#include <queue>
#include <OsiClpSolverInterface.hpp>
#include <CoinBuild.hpp>
#include <omp.h>

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
struct _CPropagation
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

	/* a pointer to original problem */
	OsiSolverInterface *solver;
};

void fillMatrices(const OsiSolverInterface *solver, CPropagation *cp);

void fixVariable(CPropagation *cp, int idx, double value, int &unfixedVars, double *colLb, double *colUb, double *constrBound,
                    int *unfixedVarsByRow);

void unfixVariable(CPropagation *cp, int idx, int &unfixedVars, double *colLb, double *colUb, double *constrBound,
                    int *unfixedVarsByRow);

CPropagation *cpropagation_create(const OsiSolverInterface *solver)
{
	CPropagation *cp = new CPropagation;
	const double *colLb = solver->getColLower();
    const double *colUb = solver->getColUpper();

    cp->solver = (OsiSolverInterface*)solver;
    fillMatrices(solver, cp);
	assert(cp->matrixByRow.size() == cp->rhs.size());

    cp->numCols = solver->getNumCols();
	cp->numRows = (int)cp->matrixByRow.size();
    cp->varIsBinary = new char[cp->numCols];
    cp->isToFix = new char[cp->numCols];
    cp->binaryVars = 0;
    cp->varsToFix = 0;

    for(int i = 0; i < cp->numCols; i++)
    {
        cp->isToFix[i] = UNFIXED;

        if(colLb[i] == 0.0 && colUb[i] == 1.0)
        {
            cp->varIsBinary[i] = 1;
            cp->binaryVars++;
        }
        else cp->varIsBinary[i] = 0;
    }

	return cp;
}

void cpropagation_free(CPropagation *cp)
{
    delete[] cp->varIsBinary;
    delete[] cp->isToFix;

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

char evaluateBound(const CPropagation *cp, const pair<int, double> &var, int constraint, double *constrBound)
{
    const double bound = constrBound[constraint];

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

char constraintPropagation(CPropagation *cp, int var, double value, vector<Fixation> &fixations,
                            int &unfixedVars, double *colLb, double *colUb, double *constrBound, int *unfixedVarsByRow)
{
	assert(value == 0.0 || value == 1.0);
	assert(var >= 0 && var < cp->numCols);
	assert(colLb[var] != colUb[var]);
	fixations.clear(); fixations.reserve(cp->numCols); //cleaning old fixations and reserving memory

	char status = NOIMPLICATION;
	queue<int> C;
	Fixation tmp;
	vector<char> rowIsInQueue(cp->numRows, 0);

	fixVariable(cp, var, value, unfixedVars, colLb, colUb, constrBound, unfixedVarsByRow);

	for(int i = 0; i < (int)cp->matrixByCol[var].size(); i++)
	{
		const pair<int, double> constraint = cp->matrixByCol[var][i];
		const int idxRow = constraint.first;
		const double coef = constraint.second;

		if((unfixedVarsByRow[idxRow] > 0) && (value == 1.0 && coef > 0.0) || (value == 0.0 && coef < 0.0))
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

			if(colLb[idxVar] == colUb[idxVar])
				continue;

			char eval = evaluateBound(cp, var, idxRow, constrBound);

            if(eval == ACTIVATE || eval == DEACTIVATE)
            {
            	assert(colLb[idxVar] != colUb[idxVar]);
            	tmp.idxVar = idxVar;
            	tmp.valueToFix = eval;
            	fixations.push_back(tmp);
            	status = FIXATION;
				unfixedVars--;
	            colLb[idxVar] = colUb[idxVar] = eval;

	            for(int j = 0; j < (int)cp->matrixByCol[idxVar].size(); j++)
	            {
	                const pair<int, double> constraint = cp->matrixByCol[idxVar][j];
	                const int idxRow2 = constraint.first;
	                const double coef2 = constraint.second;

	                unfixedVarsByRow[idxRow2]--;

	                if(eval == ACTIVATE && coef2 > 0.0)
	                {
	                    constrBound[idxRow2] -= coef2;
	                    if(unfixedVarsByRow[idxRow2] > 0 && !rowIsInQueue[idxRow2])
	                    {
	                    	C.push(idxRow2);
	                    	rowIsInQueue[idxRow2] = 1;
	                    }
	                }
	                else if(eval == DEACTIVATE && coef2 < 0.0)
	                {
	                    constrBound[idxRow2] += coef2;
	                    if(unfixedVarsByRow[idxRow2] > 0 && !rowIsInQueue[idxRow2])
	                    {
	                    	C.push(idxRow2);
	                    	rowIsInQueue[idxRow2] = 1;
	                    }
	                }
                    
                    if(unfixedVarsByRow[idxRow2] == 0 && constrBound[idxRow2] < 0.0)
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
}

void fixVariable(CPropagation *cp, int idx, double value, int &unfixedVars, double *colLb, double *colUb, double *constrBound,
                    int *unfixedVarsByRow)
{
	assert(value == 0.0 || value == 1.0);
	assert(idx >= 0 && idx < cp->numCols);
	assert(colLb[idx] == 0.0 && colUb[idx] == 1.0);

	unfixedVars--;
	colLb[idx] = colUb[idx] = value;

	for(int j = 0; j < (int)cp->matrixByCol[idx].size(); j++)
	{
		const pair<int, double> constraint = cp->matrixByCol[idx][j];
		const int idxRow = constraint.first;
		const double coef = constraint.second;

		unfixedVarsByRow[idxRow]--;

		if(value == 1.0 && coef > 0.0)
			constrBound[idxRow] -= coef;
		else if(value == 0.0 && coef < 0.0)
			constrBound[idxRow] += coef;
	}
}

void unfixVariable(CPropagation *cp, int idx, int &unfixedVars, double *colLb, double *colUb, double *constrBound,
                    int *unfixedVarsByRow)
{
    assert(colLb[idx] == colUb[idx]);

    for(int j = 0; j < (int)cp->matrixByCol[idx].size(); j++)
    {
        const pair<int, double> constraint = cp->matrixByCol[idx][j];
        const int idxRow = constraint.first;
        const double coef = constraint.second;

        unfixedVarsByRow[idxRow]++;

        if(colLb[idx] == 1.0 && coef > 0.0) //fixed in 1
            constrBound[idxRow] += coef;
        else if(colUb[idx] == 0.0 && coef < 0.0) //fixex in 0
            constrBound[idxRow] -= coef;
    }

    unfixedVars++;
    colLb[idx] = 0.0;
    colUb[idx] = 1.0;
}

void cpropagation_get_vars_to_fix(CPropagation *cp)
{
    char status;
    vector<Fixation> fixations;
    const double *lb = cp->solver->getColLower(), *ub = cp->solver->getColUpper();
    double colLb[cp->numCols], colUb[cp->numCols];
    double constrBound[cp->numRows];
    int unfixedVars = cp->binaryVars, unfixedVarsByRow[cp->numRows];

    for(int i = 0; i < cp->numCols; i++)
    {
        colLb[i] = lb[i];
        colUb[i] = ub[i];
    }

    /* calculating the initial lower bound for each constraint */
    for(int i = 0; i < cp->numRows; i++)
    {
        constrBound[i] = cp->rhs[i];
        unfixedVarsByRow[i] = (int)cp->matrixByRow[i].size();

        for(int j = 0; j < (int)cp->matrixByRow[i].size(); j++)
        {
            const int idx = cp->matrixByRow[i][j].first;
            const double coef = cp->matrixByRow[i][j].second;

            if(coef > 0.0) constrBound[i] -= (coef * colLb[idx]);
            else constrBound[i] -= (coef * colUb[idx]);
        }
    }

    #pragma omp parallel for shared(cp) \
    						 firstprivate(colLb, colUb, constrBound, unfixedVarsByRow, unfixedVars) \
    						 private(status, fixations)
        for(int i = 0; i < cp->numCols; i++)
        {
            if(!cp->varIsBinary[i]) continue;

            status = constraintPropagation(cp, i, 0.0, fixations, unfixedVars, colLb, colUb, constrBound, unfixedVarsByRow);
            if(status == CONFLICT)
            {
                cp->isToFix[i] = ACTIVATE;
                #pragma omp atomic
                	cp->varsToFix++;
            }
            unfixVariable(cp, i, unfixedVars, colLb, colUb, constrBound, unfixedVarsByRow);
            for(int j = 0; j < (int)fixations.size(); j++)
    			unfixVariable(cp, fixations[j].idxVar, unfixedVars, colLb, colUb, constrBound, unfixedVarsByRow);

            status = constraintPropagation(cp, i, 1.0, fixations, unfixedVars, colLb, colUb, constrBound, unfixedVarsByRow);
            if(status == CONFLICT)
            {
                cp->isToFix[i] = DEACTIVATE;
                #pragma omp atomic
                	cp->varsToFix++;
            }
            unfixVariable(cp, i, unfixedVars, colLb, colUb, constrBound, unfixedVarsByRow);
            for(int j = 0; j < (int)fixations.size(); j++)
    			unfixVariable(cp, fixations[j].idxVar, unfixedVars, colLb, colUb, constrBound, unfixedVarsByRow);
        }
}

OsiSolverInterface* cpropagation_preprocess(CPropagation *cp, int nindexes[])
{
    if(cp->varsToFix == 0)
    {
        /* printf("There are no variables to remove from the problem!\n"); */
        return cp->solver; /* returns a pointer to original solver */
    }

    const double *rhs = cp->solver->getRightHandSide();
    const double *colLb = cp->solver->getColLower(), *colUb = cp->solver->getColUpper();
    const double *rowLb = cp->solver->getRowLower(), *rowUb = cp->solver->getRowUpper();
    const double *objCoef = cp->solver->getObjCoefficients();
    const char *sense = cp->solver->getRowSense();
    const char *ctype = cp->solver->getColType();
    const CoinPackedMatrix *M = cp->solver->getMatrixByRow();

    double sumFixedObj = 0.0; /* stores the sum of objective coefficients of all variables fixed to 1 */

    OsiSolverInterface *preProcSolver = new OsiClpSolverInterface();
    preProcSolver->setIntParam(OsiNameDiscipline, 2);
    preProcSolver->messageHandler()->setLogLevel(0);
    preProcSolver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    preProcSolver->setObjName(cp->solver->getObjName());

    for(int i = 0, j = 0; i < cp->solver->getNumCols(); i++)
    {
        nindexes[i] = -1;
        if(cp->isToFix[i] == UNFIXED)
        {
            preProcSolver->addCol(0, NULL, NULL, colLb[i], colUb[i], objCoef[i]);
            preProcSolver->setColName(j, cp->solver->getColName(i));
            if(ctype[i] == 1 || ctype[i] == 2)
            	preProcSolver->setInteger(j);
            else preProcSolver->setContinuous(j);
            nindexes[i] = j++;
        }
        else if(cp->isToFix[i] == ACTIVATE)
            sumFixedObj += objCoef[i];
    }

    /* adding a variable with cost equals to the sum of all coefficients of variables fixed to 1 */
    preProcSolver->addCol(0, NULL, NULL, 1.0, 1.0, sumFixedObj);
    preProcSolver->setColName(preProcSolver->getNumCols()-1, "sumFixedObj");
    preProcSolver->setInteger(preProcSolver->getNumCols()-1);

    for(int idxRow = 0; idxRow < cp->solver->getNumRows(); idxRow++)
    {
        const CoinShallowPackedVector &row = M->getVector(idxRow);
        const int nElements = row.getNumElements();
        const int *idxs = row.getIndices();
        const double *coefs = row.getElements();
        vector< int > vidx; vidx.reserve(cp->solver->getNumCols());
        vector< double > vcoef; vcoef.reserve(cp->solver->getNumCols());
        double activeCoefs = 0.0;

        for(int i = 0; i < nElements; i++)
        {
            if(cp->isToFix[idxs[i]] == UNFIXED)
            {
                assert(nindexes[idxs[i]] >= 0 && nindexes[idxs[i]] < cp->solver->getNumCols());
                vidx.push_back(nindexes[idxs[i]]);
                vcoef.push_back(coefs[i]);
            }
            else if(cp->isToFix[idxs[i]] == ACTIVATE)
            	activeCoefs += objCoef[idxs[i]];
        }

        if(!vidx.empty())
        {
        	double rlb, rub;
        	
        	if(sense[idxRow] == 'E')
            {
                rlb = rowLb[idxRow] - activeCoefs;
                rub = rowUb[idxRow] - activeCoefs;
            }
        	else if(sense[idxRow] == 'L')
            {
                rlb = rowLb[idxRow];
                rub = rowUb[idxRow] - activeCoefs;
            }
        	else if(sense[idxRow] == 'G')
            {
                rlb = rowLb[idxRow] - activeCoefs;
                rub = rowUb[idxRow];
            }
        	else
        	{
        		fprintf(stderr, "Error: invalid type of constraint!\n");
        		exit(EXIT_FAILURE);
        	}

            preProcSolver->addRow((int)vcoef.size(), &vidx[0], &vcoef[0], rlb, rub);
            preProcSolver->setRowName(idxRow, cp->solver->getRowName(idxRow));
        }
	}

    return preProcSolver;
}

int cpropagation_get_num_vars_to_fix(CPropagation *cp) { return cp->varsToFix; }