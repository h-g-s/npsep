extern "C"
{
    #include "memory.h"
}

#include "problem.h"
#include <OsiClpSolverInterface.hpp>
#include <CoinBuild.hpp>

#define MAX_NAME_SIZE 64

struct _Problem
{
    int numCols, numRows, numElements;
    double infty; /* stores the solver's value for infinity */
    double *colLb, *colUb, *rhs, *objCoef;
    char *rowSense, *colType, **colName, **rowName;

    /* stores matrices of indexes and coefficients by row */
    int **idxsByRow, *rowNElements;
    double **coefsByRow;

    /* stores matrices of indexes and coefficients by column */
    int **idxsByCol, *colNElements;
    double **coefsByCol;
};

Problem* problem_create(int numCols, int numRows, double infty)
{
    Problem *p = (Problem*) xmalloc (sizeof(Problem));

    p->numCols = numCols;
    p->numRows = numRows;
    p->numElements = 0;
    p->infty = infty;
    p->colLb = (double*) xmalloc(sizeof(double) * p->numCols);
    p->colUb = (double*) xmalloc(sizeof(double) * p->numCols);
    p->rhs = (double*) xmalloc(sizeof(double) * p->numRows);
    p->objCoef = (double*) xmalloc(sizeof(double) * p->numCols);
    p->rowSense = (char*) xmalloc(sizeof(char) * p->numRows);
    p->colType = (char*) xmalloc(sizeof(char) * p->numCols);
    p->colName = (char**) xmalloc(sizeof(char*) * p->numCols);
    p->rowName = (char**) xmalloc(sizeof(char*) * p->numRows);
    
    p->idxsByRow = (int**) xmalloc(sizeof(int*) * p->numRows);
    p->rowNElements = (int*) xmalloc(sizeof(int) * p->numRows);
    p->coefsByRow = (double**) xmalloc(sizeof(double*) * p->numRows);

    p->idxsByCol = (int**) xmalloc(sizeof(int*) * p->numCols);
    p->colNElements = (int*) xmalloc(sizeof(int) * p->numCols);
    p->coefsByCol = (double**) xmalloc(sizeof(double*) * p->numCols);

    int i;
    for(i = 0; i < p->numCols; i++)
    {
        p->colName[i] = (char*) xmalloc(sizeof(char) * MAX_NAME_SIZE);
        p->colNElements[i] = 0;
    }
    for(i = 0; i < p->numRows; i++)
    {
        p->rowName[i] = (char*) xmalloc(sizeof(char) * MAX_NAME_SIZE);
        p->rowNElements[i] = 0;
    }

    return p;
}

Problem* problem_create_using_osi(const OsiSolverInterface *solver)
{
    Problem *p = (Problem*) xmalloc (sizeof(Problem));

    p->numCols = solver->getNumCols();
    p->numRows = solver->getNumRows();
    p->numElements = solver->getNumElements();
    p->infty = solver->getInfinity();
    p->colLb = (double*) xmalloc(sizeof(double) * p->numCols);
    p->colUb = (double*) xmalloc(sizeof(double) * p->numCols);
    p->rhs = (double*) xmalloc(sizeof(double) * p->numRows);
    p->objCoef = (double*) xmalloc(sizeof(double) * p->numCols);
    p->rowSense = (char*) xmalloc(sizeof(char) * p->numRows);
    p->colType = (char*) xmalloc(sizeof(char) * p->numCols);
    p->colName = (char**) xmalloc(sizeof(char*) * p->numCols);
    p->rowName = (char**) xmalloc(sizeof(char*) * p->numRows);
    
    p->idxsByRow = (int**) xmalloc(sizeof(int*) * p->numRows);
    p->rowNElements = (int*) xmalloc(sizeof(int) * p->numRows);
    p->coefsByRow = (double**) xmalloc(sizeof(double*) * p->numRows);

    p->idxsByCol = (int**) xmalloc(sizeof(int*) * p->numCols);
    p->colNElements = (int*) xmalloc(sizeof(int) * p->numCols);
    p->coefsByCol = (double**) xmalloc(sizeof(double*) * p->numCols);

    const CoinPackedMatrix *mRow = solver->getMatrixByRow(), *mCol = solver->getMatrixByCol();
    const double *rhs = solver->getRightHandSide();
    const double *colLb = solver->getColLower(), *colUb = solver->getColUpper();
    const double *rowLb = solver->getRowLower(), *rowUb = solver->getRowUpper();
    const char *rowSense = solver->getRowSense();
    const double *objCoef = solver->getObjCoefficients();
    const char *colType = solver->getColType();
    int i, j;

    for(i = 0; i < p->numCols; i++)
    {
        p->colLb[i] = colLb[i];
        p->colUb[i] = colUb[i];
        p->objCoef[i] = objCoef[i];
        p->colType[i] = colType[i];

        p->colName[i] = (char*) xmalloc(sizeof(char) * MAX_NAME_SIZE);
        strncpy(p->colName[i], solver->getColName(i).c_str(), MAX_NAME_SIZE);
    }

    for(i = 0; i < p->numRows; i++)
    {
        const CoinShallowPackedVector &row = mRow->getVector(i);
        const int *idxs = row.getIndices();
        const double *coefs = row.getElements();

        p->rowNElements[i] = row.getNumElements();
        p->idxsByRow[i] = (int*) xmalloc(sizeof(int) * p->rowNElements[i]);
        p->coefsByRow[i] = (double*) xmalloc(sizeof(double) * p->rowNElements[i]);
        p->rowName[i] = (char*) xmalloc(sizeof(char) * MAX_NAME_SIZE);
        strncpy(p->rowName[i], solver->getRowName(i).c_str(), MAX_NAME_SIZE);
        p->rhs[i] = rhs[i];
        p->rowSense[i] = rowSense[i];

        for(j = 0; j < p->rowNElements[i]; j++)
        {
            p->idxsByRow[i][j] = idxs[j];
            p->coefsByRow[i][j] = coefs[j];
        }
    }

    for(i = 0; i < p->numCols; i++)
    {
        const CoinShallowPackedVector &col = mCol->getVector(i);
        const int *idxs = col.getIndices();
        const double *coefs = col.getElements();

        p->colNElements[i] = col.getNumElements();
        p->idxsByCol[i] = (int*) xmalloc(sizeof(int) * p->colNElements[i]);
        p->coefsByCol[i] = (double*) xmalloc(sizeof(double) * p->colNElements[i]);

        for(j = 0; j < p->colNElements[i]; j++)
        {
            p->idxsByCol[i][j] = idxs[j];
            p->coefsByCol[i][j] = coefs[j];
        }
    }

    return p;
}

void problem_free(Problem **p)
{
    int i;

    free((*p)->colLb);
    free((*p)->colUb);
    free((*p)->rhs);
    free((*p)->objCoef);
    free((*p)->rowSense);
    free((*p)->colType);
    
    for(i = 0; i < (*p)->numCols; i++)
    {
        free((*p)->idxsByCol[i]);
        free((*p)->coefsByCol[i]);
        free((*p)->colName[i]);
    }
    free((*p)->idxsByCol);
    free((*p)->coefsByCol);
    free((*p)->colName);
    free((*p)->colNElements);


    for(i = 0; i < (*p)->numRows; i++)
    {
        free((*p)->idxsByRow[i]);
        free((*p)->coefsByRow[i]);
        free((*p)->rowName[i]);
    }
    free((*p)->idxsByRow);
    free((*p)->coefsByRow);
    free((*p)->rowName);
    free((*p)->rowNElements);
    

    free(*p);
    *p = NULL;
}

int problem_num_cols(const Problem *p) { return p->numCols; }
int problem_num_rows(const Problem *p) { return p->numRows; }
int problem_num_elements(const Problem *p) { return p->numElements; }
double problem_get_infinity(const Problem *p) { return p->infty; }

int problem_row_size(const Problem *p, int idxRow)
{
    assert(idxRow >= 0 && idxRow < p->numRows);
    return p->rowNElements[idxRow];
}

char problem_row_sense(const Problem *p, int idxRow)
{
    assert(idxRow >= 0 && idxRow < p->numRows);
    return p->rowSense[idxRow];
}

const int* problem_row_idxs(const Problem *p, int idxRow)
{
    assert(idxRow >= 0 && idxRow < p->numRows);
    return p->idxsByRow[idxRow];
}

const double* problem_row_coefs(const Problem *p, int idxRow)
{
    assert(idxRow >= 0 && idxRow < p->numRows);
    return p->coefsByRow[idxRow];
}

double problem_row_rhs(const Problem *p, int idxRow)
{
    assert(idxRow >= 0 && idxRow < p->numRows);
    return p->rhs[idxRow];
}

const char* problem_row_name(const Problem *p, int idxRow)
{
    assert(idxRow >= 0 && idxRow < p->numRows);
    return p->rowName[idxRow];
}

char problem_var_is_binary(const Problem *p, int idxVar)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    return (p->colType[idxVar] == BINARY);
}

char problem_var_type(const Problem *p, int idxVar)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    return p->colType[idxVar];   
}

const char* problem_var_name(const Problem *p, int idxVar)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    return p->colName[idxVar];
}

double problem_var_lower_bound(const Problem *p, int idxVar)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    return p->colLb[idxVar];
}

double problem_var_upper_bound(const Problem *p, int idxVar)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    return p->colUb[idxVar];
}

double problem_var_obj_coef(const Problem *p, int idxVar)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    return p->objCoef[idxVar];
}

int problem_var_n_rows(const Problem *p, int idxVar)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    return p->colNElements[idxVar];
}

const int* problem_var_rows_idxs(const Problem *p, int idxVar)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    return p->idxsByCol[idxVar];
}

const double* problem_var_rows_coefs(const Problem *p, int idxVar)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    return p->coefsByCol[idxVar];
}

const double* problem_vars_lower_bound(const Problem *p) { return p->colLb; }

const double* problem_vars_upper_bound(const Problem *p) { return p->colUb; }

const double* problem_vars_obj_coefs(const Problem *p) { return p->objCoef; }

void problem_set_num_cols(Problem *p, double numCols) { p->numCols = numCols; }
void problem_set_num_rows(Problem *p, double numRows) { p->numRows = numRows; }
void problem_set_num_elements(Problem *p, int nElements) { p->numElements = nElements; }
void problem_set_infinity(Problem *p, double infty) { p->infty = infty; }

void problem_var_set_lower_bound(Problem *p, int idxVar, double value)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    p->colLb[idxVar] = value;
}

void problem_var_set_upper_bound(Problem *p, int idxVar, double value)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    p->colUb[idxVar] = value;
}

void problem_var_set_obj_coef(Problem *p, int idxVar, double value)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    p->objCoef[idxVar] = value;
}

void problem_var_set_type(Problem *p, int idxVar, char type)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    p->colType[idxVar] = type;  
}

void problem_var_set_name(Problem *p, int idxVar, const char *name)
{
    assert(idxVar >= 0 && idxVar < p->numCols);
    strncpy(p->colName[idxVar], name, MAX_NAME_SIZE);
}

void problem_set_row(Problem *p, int idxRow, const int *idxs, const double *coefs, int nElements, double rhs, char sense)
{
    assert(idxRow >= 0 && idxRow < p->numRows && nElements > 0);

    int i;

    if(p->rowNElements[idxRow] == 0)
    {
        p->idxsByRow[idxRow] = (int*) xmalloc(sizeof(int) * nElements);
        p->coefsByRow[idxRow] = (double*) xmalloc(sizeof(double) * nElements);
        p->numElements += nElements;
    }
    else
    {
        for(i = 0; i < p->rowNElements[idxRow]; i++)
        {
            const int idx = p->idxsByRow[idxRow][i];
            p->colNElements[idx]--;
        }
        p->idxsByRow[idxRow] = (int*) xrealloc(p->idxsByRow[idxRow], sizeof(int) * nElements);
        p->coefsByRow[idxRow] = (double*) xrealloc(p->coefsByRow[idxRow], sizeof(double) * nElements);  
        p->numElements += (nElements - p->rowNElements[idxRow]);
    }

    p->rowNElements[idxRow] = nElements;
    p->rhs[idxRow] = rhs;
    p->rowSense[idxRow] = sense;

    for(i = 0; i < nElements; i++)
    {
        p->idxsByRow[idxRow][i] = idxs[i];
        p->coefsByRow[idxRow][i] = coefs[i];
        p->colNElements[idxs[i]]++;
    }
}

void problem_row_set_name(Problem *p, int idxRow, const char *name)
{
    assert(idxRow >= 0 && idxRow < p->numRows);
    strncpy(p->rowName[idxRow], name, MAX_NAME_SIZE);
}

void problem_update_matrices_by_col(Problem *p)
{
    int i, j, position[p->numCols];

    memset(position, 0, sizeof(int) * p->numCols);

    for(i = 0; i < p->numCols; i++)
    {
        p->idxsByCol[i] = (int*) xmalloc(sizeof(int) * p->colNElements[i]);
        p->coefsByCol[i] = (double*) xmalloc(sizeof(double) * p->colNElements[i]);
    }

    for(i = 0; i < p->numRows; i++)
        for(j = 0; j < p->rowNElements[i]; j++)
        {
            const int idx = p->idxsByRow[i][j];
            const double coef = p->coefsByRow[i][j];
            p->idxsByCol[idx][position[idx]] = i;
            p->coefsByCol[idx][position[idx]] = coef;
            position[idx]++;
        }
}

OsiSolverInterface* problem_convert_to_osi(Problem *p)
{
    int i;
    double rowLb, rowUb;
    OsiSolverInterface *solver = new OsiClpSolverInterface();
    CoinBuild cb;

    solver->setIntParam(OsiNameDiscipline, 2);
    solver->messageHandler()->setLogLevel(0);
    solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);

    for(i = 0; i < p->numCols; i++)
    {
        solver->addCol(0, NULL, NULL, p->colLb[i], p->colUb[i], p->objCoef[i]);
        solver->setColName(i, p->colName[i]);
        
        if(p->colType[i] == CONTINUOUS)
            solver->setContinuous(i);
        else 
            solver->setInteger(i);
    }

    for(i = 0; i < p->numRows; i++)
    {
        switch(p->rowSense[i])
        {
            case 'E':
                rowLb = p->rhs[i];
                rowUb = p->rhs[i];
            break;

            case 'L':
                rowLb = -p->infty;
                rowUb = p->rhs[i];
            break;

            case 'G':
                rowLb = p->rhs[i];
                rowUb = p->infty;
            break;

            default:
                fprintf(stderr, "Error: invalid type of constraint!\n");
                exit(EXIT_FAILURE);
        }

        cb.addRow(p->rowNElements[i], p->idxsByRow[i], p->coefsByRow[i], rowLb, rowUb);
    }

    solver->addRows(cb);
    
    for(i = 0; i < p->numRows; i++)
        solver->setRowName(i, p->rowName[i]);

    return solver;
}