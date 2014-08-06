#ifndef CGLECLIQUE_HPP
#define CGLECLIQUE_HPP

#include <CglCutGenerator.hpp>
#include <CglCutGenerator.hpp>
extern "C"
{
#include "clique_separation.h"
#include "cgraph.h"
}


class CglEClique : public CglCutGenerator
{
public:
   CglEClique();

   /// Copy constructor
   CglEClique(const CglEClique& rhs);

   /// Clone
   virtual CglCutGenerator * clone() const;

   /* calling setCGraph before branch and bound
      speedup a lot separation passes since the same
      conflict graph is used and not rebuilt every time */
   void setCGraph( const CGraph *cg );

   void parseParameters( int _argc, const char **_argv );

   void setGenOddHoles( const bool gen );

   //virtual void generateCuts( const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info=CglTreeInfo() ) const ;

   virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
         const CglTreeInfo info = CglTreeInfo() );


   virtual ~CglEClique();

   /* for debugging purposes */

   const std::vector< std::string > *colNames;

private:
   const CGraph *_cgraph;

   int argc;
   char **argv;

   bool genOddHoles;
};

#endif // CGLECLIQUE_HPP

