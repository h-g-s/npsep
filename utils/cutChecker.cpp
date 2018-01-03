#include <cstdio>
#include <cstdlib>
#include <OsiClpSolverInterface.hpp>
#include <vector>

using namespace std;

#define EPS 1e-6

using namespace std;

OsiSolverInterface *solver;
OsiClpSolverInterface *realSolver;

double readSolution(const char* filename, vector<int>& idxs, vector<double>& sol, vector<string>& names)
{
	FILE *file = fopen(filename, "rt");
	char line[500], varName[100];
	int varIndex, numCols = solver->getNumCols();;
	double obj = .0, varValue, varCost;
   const double *objCoefs = solver->getObjCoefficients();
   const vector<string> origNames = solver->getColNames();

	if(!fgets(line, 500, file)) //discarding the first line of cbc solution file
	{
		printf("Empty file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(line, 500, file) != NULL)
	{
      bool found = false;
		sscanf(line, "%d %s %lf %lf", &varIndex, varName, &varValue, &varCost);

      for(int i = 0; i < numCols; i++) //find the variable (the indexes can be different depending the solver)
         if(origNames[i] == varName)
         {
            found = true;
            obj += (varValue * objCoefs[i]);
            idxs.push_back(i);
            sol.push_back(varValue);
            names.push_back(origNames[i]);
            break;
         }

      if(!found)
      {
         printf("Variable not found!\n");
         exit(EXIT_FAILURE);
      }
	}

	fclose(file);

	return obj;
}

int main( int argc, char **argv )
{
	if(argc != 3)
	{
		printf("Invalid number of parameters! Required: 3.\n");
		exit(EXIT_FAILURE);
	}

	solver = NULL;
	realSolver = new OsiClpSolverInterface();
	solver = (OsiSolverInterface*) realSolver;
	solver->setIntParam(OsiNameDiscipline, 2);
	solver->readLp(argv[1]);
	solver->setIntParam(OsiNameDiscipline, 2);
  	solver->messageHandler()->setLogLevel(1);
	solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);

	const int numCols = solver->getNumCols();
	double obj;
	vector<int> idxs;
	vector<double> sol;
	vector<string> names;
	
	obj = readSolution(argv[2], idxs, sol, names);
	assert(idxs.size() == sol.size() && idxs.size() == names.size());

	for(int i = 0; i < numCols; i++)
	{
		solver->setColLower(i, .0);
		solver->setColUpper(i, .0);
	}

	for(int i = 0; i < (int)names.size(); i++) //fixing variables with values that appear in the solution file
	{
      	int idx = idxs[i];
		double value = sol[i];
		string name = names[i];

		solver->setColLower(idx, value);
		solver->setColUpper(idx, value);
	}

   solver->initialSolve();

   if(!solver->isProvenOptimal())
   {
      printf("Relaxation says problem is infeasible!\n");
      exit(EXIT_FAILURE);
   }

	solver->branchAndBound();

   if(!solver->isProvenOptimal())
   {
      printf("B&B says problem is infeasible!\n");
      exit(EXIT_FAILURE);
   }

	if(fabs(solver->getObjValue() - obj) > EPS)
	{
		printf("Different values of objective function: %lf \t %lf\n", solver->getObjValue(), obj);
		exit(EXIT_FAILURE);
	}

	delete realSolver;

	printf("\nSUCCESS!\n");

	return EXIT_SUCCESS;
}
