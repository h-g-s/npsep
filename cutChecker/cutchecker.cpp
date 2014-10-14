#include <cstdio>
#include <cstdlib>
#include <OsiCbcSolverInterface.hpp>
#include <vector>

using namespace std;

#define EPS 1e-6

using namespace std;

double readSolution(const char* filename, vector<int>& idxs, vector<double>& sol, vector<string>& names)
{
	FILE *file = fopen(filename, "rt");
	char line[500], varName[100];
	int varIndex;
	double obj = .0, varValue, varCost;

	if(!fgets(line, 500, file)) //discarding the first line of cbc solution file
	{
		printf("Empty file!\n");
		exit(EXIT_FAILURE);
	}

   	while(fgets(line, 500, file) != NULL)
	{
		sscanf(line, "%d %s %lf %lf", &varIndex, varName, &varValue, &varCost);
		//printf("%d %s %lf %lf	\n", varIndex, varName, varValue, varCost);
		obj += (varValue * varCost);
		idxs.push_back(varIndex);
		sol.push_back(varValue);
		names.push_back(varName);
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

	OsiSolverInterface *solver = NULL;
	OsiCbcSolverInterface *realSolver = new OsiCbcSolverInterface();
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

   	for(int i = 0; i < numCols; i++) //all variables are assigned with zero
   	{
   		solver->setColLower(i, .0);
   		solver->setColUpper(i, .0);
   	}

   	for(int i = 0; i < (int)idxs.size(); i++) //fixing variables with values that appear in the solution file
   	{
   		const int idx = idxs[i];
   		const double value = sol[i];
   		const string name = names[i];
   		string realName = solver->getColName(idx);

   		if(name != realName)
   		{
   			printf("Different variables: %s \t %s!\n", realName.c_str(), name.c_str());
   			exit(EXIT_FAILURE);
   		}

   		solver->setColLower(idx, value);
   		solver->setColUpper(idx, value);
   	}

   	solver->initialSolve();
   	solver->branchAndBound();

   	if(!solver->isProvenOptimal() || fabs(solver->getObjValue() - obj) > EPS)
   	{
   		printf("Invalid cuts were inserted!\n");
   		exit(EXIT_FAILURE);
   	}

   	delete realSolver;

	return EXIT_SUCCESS;
}
