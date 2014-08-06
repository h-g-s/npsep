#ifndef MIPSEP_HPP
#define MIPSEP_HPP

extern "C"
{
#include "cgraph.h"
}

class OsiSolverInterface;
class OsiClpSolverInterface;
class MPIndex;
class CbcModel;
class CglPreProcess;
#include <vector>
#include <set>
#include <string>

class MIPSep
{
   public:
      MIPSep( CGraph *_cgraph, const int _w[],
              const std::vector<std::string> &_colNames, const char *_probName );

      void writeLP();

      const std::set< std::set< int > > *cliquesFound() { return &cliques; }

      char problemName[256];

      void solve();

      virtual ~MIPSep();

   private:
      // best Solution
      const std::vector< int > &getSelectedNodes() { return selectedNodes; }

      void saveSolutions( CbcModel *cbc, CglPreProcess *preProc );

      double timeRelax;

      std::set< std::set< int > > cliques;

      std::string probName;

      std::vector< int > selectedNodes;

      CGraph *cgraph;

      int *w;

      void createXVars();
      void createYVars();
      void createConsActivation();

      void addBinVars( std::vector< std::string > &colNames, std::vector< double > &obj );
      void makeIntegers( int colStart, int nCols );


      OsiSolverInterface *lp;
      OsiClpSolverInterface *clp;

      std::vector< std::string > colNames;

      MPIndex *colIndex;
};

#endif // MIPSEP_HPP
