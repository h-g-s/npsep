#include "constraint_propagation.h"
#include <limits>
#include <cfloat>
#include <cmath>
#include <vector>
#include <set>
#include <queue>
#include <CoinBuild.hpp>
#include <omp.h>

#include <OsiClpSolverInterface.hpp>

extern "C"
{
	#include "memory.h"
}

void fillMatrices(CPropagation *cp);

void fixVariable(CPropagation *cp, int idx, double value, int &unfixedVars, double *colLb, double *colUb, double *constrBound,
                    int *unfixedVarsByRow);

void unfixVariable(CPropagation *cp, int idx, int &unfixedVars, double *colLb, double *colUb, double *constrBound,
                    int *unfixedVarsByRow);

void calculateCj(CPropagation *cp);

CPropagation *cpropagation_create(const Problem *_problem)
{
	CPropagation *cp = new CPropagation;
	const double *colLb = problem_vars_lower_bound(_problem);
    const double *colUb = problem_vars_upper_bound(_problem);

    cp->problem = (Problem*)_problem;
    fillMatrices(cp);
	assert(cp->matrixByRow.size() == cp->rhs.size());

    cp->numCols = problem_num_cols(cp->problem);
	cp->numRows = (int)cp->matrixByRow.size();
    cp->varIsBinary = new char[cp->numCols];
    cp->isToFix = new char[cp->numCols];
    cp->binaryVars = 0;
    cp->varsToFix = 0;

    int numThreads = omp_get_max_threads();
    int chunk = max(1, cp->numCols/numThreads);
    #pragma omp parallel for shared(cp, colLb, colUb) num_threads(numThreads) schedule(static, chunk)
        for(int i = 0; i < cp->numCols; i++)
        {
            cp->isToFix[i] = UNFIXED;

            if(problem_var_is_binary(cp->problem, i))
            {
                cp->varIsBinary[i] = 1;
                #pragma omp atomic
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

void fillMatrices(CPropagation *cp)
{
	const char *ctype = problem_vars_type(cp->problem);

    for(int idxRow = 0; idxRow < problem_num_rows(cp->problem); idxRow++)
    {
    	char sense = problem_row_sense(cp->problem, idxRow);

        if(sense == 'R') /* ignoring ranged constarints */
            continue;

        const int nElements = problem_row_size(cp->problem, idxRow);
        const int *idxs = problem_row_idxs(cp->problem, idxRow);
        const double *coefs = problem_row_coefs(cp->problem, idxRow);
        const double rhs = problem_row_rhs(cp->problem, idxRow);
        double mult = (sense == 'G') ? -1.0 : 1.0;
        vector<pair<int, double> > constraint(nElements);
        bool allBinaries = true;

        for(int i = 0; i < nElements; i++)
        {
            constraint[i] = pair<int, double> (idxs[i], mult * coefs[i]);
            if(ctype[idxs[i]] != BINARY)
            {
                allBinaries = false;
                break;
            }
        }

        if(allBinaries)
        {
            cp->matrixByRow.push_back(constraint);
            cp->rhs.push_back(rhs * mult);
            if(sense == 'E')
            {
                for(int j = 0; j < nElements; j++)
                    constraint[j].second = -1.0 * constraint[j].second;
                cp->matrixByRow.push_back(constraint);
                cp->rhs.push_back(-1.0 * rhs);
            }
        }
    }

    cp->matrixByCol.resize(problem_num_cols(cp->problem));
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
	//assert(colLb[var] != colUb[var]);
    if(colLb[var] == colUb[var])
    {
        // assert(colLb[var] == value);
        return NOIMPLICATION;
    }
	//assert(fixations.empty());

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
	// assert(colLb[idx] == 0.0 && colUb[idx] == 1.0);
    if(colLb[idx] == colUb[idx])
    {
        //assert(colLb[idx] == value);
        return;
    }

	unfixedVars--;
	colLb[idx] = colUb[idx] = value;

    int numThreads = omp_get_max_threads();
    int chunk = max(1, (int)cp->matrixByCol[idx].size()/numThreads);
    #pragma omp parallel for shared(cp, colLb, colUb, unfixedVarsByRow, constrBound) num_threads(numThreads) \
                         schedule(static, chunk)
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
    if(colLb[idx] != colUb[idx])
        return;
    //assert(colLb[idx] == colUb[idx]);

    int numThreads = omp_get_max_threads();
    int chunk = max(1, (int)cp->matrixByCol[idx].size()/numThreads);
    #pragma omp parallel for shared(cp, colLb, colUb, unfixedVarsByRow, constrBound) num_threads(numThreads) \
                         schedule(static, chunk)
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

void cpropagation_get_vars_to_fix(CPropagation *cp, const CGraph* cgraph)
{
    const double *lb = problem_vars_lower_bound(cp->problem), *ub = problem_vars_upper_bound(cp->problem);
    double colLb[cp->numCols], colUb[cp->numCols], constrBound[cp->numRows];
    int unfixedVars = cp->binaryVars, unfixedVarsByRow[cp->numRows];

    int numThreads = omp_get_max_threads();
    int chunk = max(1, cp->numCols/numThreads);

    #pragma omp parallel for shared(colLb, colUb, lb, ub) num_threads(numThreads) schedule(static, chunk)
        for(int i = 0; i < cp->numCols; i++)
        {
            colLb[i] = lb[i];
            colUb[i] = ub[i];
        }

    /* calculating the initial lower bound for each constraint */
    #pragma omp parallel for shared(colLb, colUb, constrBound, cp, unfixedVarsByRow) num_threads(numThreads) \
                             schedule(dynamic)
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

    #pragma omp parallel for shared(cp) firstprivate(colLb, colUb, constrBound, unfixedVarsByRow, unfixedVars) \
    						 num_threads(numThreads) schedule(dynamic)
        for(int i = 0; i < cp->numCols; i++)
        {
            if(!cp->varIsBinary[i]) continue;

            if(colLb[i] == colUb[i])
            {
                // cp->isToFix[i] = (int)colLb[i];
                continue;
            }

            vector<Fixation> fixations; fixations.reserve(cp->numCols);
            char status;

            if(cgraph_degree(cgraph, cp->numCols + i) > 1)
            {
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
	    		fixations.clear(); fixations.reserve(cp->numCols);
	    	}

	    	if(cgraph_degree(cgraph, i) > 1)
	    	{
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
}

int cpropagation_get_num_vars_to_fix(CPropagation *cp) { return cp->varsToFix; }
char cpropagation_var_is_to_fix(CPropagation *cp, int idxVar)
{
	assert(idxVar >= 0 && idxVar < problem_num_cols(cp->problem));
	return cp->isToFix[idxVar];
}


OsiSolverInterface* cpropagation_preprocess(CPropagation *cp, int nindexes[])
{
    if(cp->varsToFix == 0)
    {
        /* printf("There are no variables to remove from the problem!\n"); */
        return NULL; /* returns a pointer to original solver */
    }

    const double *colLb = problem_vars_lower_bound(cp->problem), *colUb = problem_vars_upper_bound(cp->problem);
    const double *objCoef = problem_vars_obj_coefs(cp->problem);
    const char *ctype = problem_vars_type(cp->problem);

    double sumFixedObj = 0.0; /* stores the sum of objective coefficients of all variables fixed to 1 */

    OsiSolverInterface *preProcSolver = new OsiClpSolverInterface();
    preProcSolver->setIntParam(OsiNameDiscipline, 2);
    preProcSolver->messageHandler()->setLogLevel(0);
    preProcSolver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    //preProcSolver->setObjName(cp->solver->getObjName());

    for(int i = 0, j = 0; i < problem_num_cols(cp->problem); i++)
    {
        nindexes[i] = -1;
        if(cp->isToFix[i] == UNFIXED)
        {
            preProcSolver->addCol(0, NULL, NULL, colLb[i], colUb[i], objCoef[i]);
            preProcSolver->setColName(j, problem_var_name(cp->problem, i));
            if(problem_var_type(cp->problem, i) == CONTINUOUS)
                preProcSolver->setContinuous(j);
            else 
                preProcSolver->setInteger(j);
            nindexes[i] = j++;
        }
        else if(cp->isToFix[i] == ACTIVATE)
            sumFixedObj += objCoef[i];
    }

    if(fabs(sumFixedObj) > EPS)
    {
        /* adding a variable with cost equals to the sum of all coefficients of variables fixed to 1 */
        preProcSolver->addCol(0, NULL, NULL, 1.0, 1.0, sumFixedObj);
        preProcSolver->setColName(preProcSolver->getNumCols()-1, "sumFixedObj");
        preProcSolver->setInteger(preProcSolver->getNumCols()-1);
    }

    for(int idxRow = 0; idxRow < problem_num_rows(cp->problem); idxRow++)
    {
        const int nElements = problem_row_size(cp->problem, idxRow);
        const int *idxs = problem_row_idxs(cp->problem, idxRow);
        const double *coefs = problem_row_coefs(cp->problem, idxRow);
        vector< int > vidx; vidx.reserve(problem_num_cols(cp->problem));
        vector< double > vcoef; vcoef.reserve(problem_num_cols(cp->problem));
        double activeCoefs = 0.0;

        for(int i = 0; i < nElements; i++)
        {
            if(cp->isToFix[idxs[i]] == UNFIXED)
            {
                assert(nindexes[idxs[i]] >= 0 && nindexes[idxs[i]] < problem_num_cols(cp->problem));
                vidx.push_back(nindexes[idxs[i]]);
                vcoef.push_back(coefs[i]);
            }
            else if(cp->isToFix[idxs[i]] == ACTIVATE)
            	activeCoefs += coefs[i];
        }

        if(!vidx.empty())
        {
        	double rlb, rub;
        	const char sense = problem_row_sense(cp->problem, idxRow);
        	
        	if(sense == 'E')
            {
                rlb = problem_row_rhs(cp->problem, idxRow) - activeCoefs;
                rub = problem_row_rhs(cp->problem, idxRow) - activeCoefs;
            }
        	else if(sense == 'L')
            {
                rlb = preProcSolver->getInfinity();
                rub = problem_row_rhs(cp->problem, idxRow) - activeCoefs;
            }
        	else if(sense == 'G')
            {
                rlb = problem_row_rhs(cp->problem, idxRow) - activeCoefs;
                rub = preProcSolver->getInfinity();
            }
        	else
        	{
        		fprintf(stderr, "Error: invalid type of constraint!\n");
        		exit(EXIT_FAILURE);
        	}

        	preProcSolver->addRow((int)vcoef.size(), &vidx[0], &vcoef[0], rlb, rub);
            preProcSolver->setRowName(idxRow, problem_row_name(cp->problem, idxRow));
        }
	}

    return preProcSolver;
}
