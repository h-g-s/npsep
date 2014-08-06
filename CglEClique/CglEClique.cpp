#include "CglEClique.hpp"
#include <cstdlib>
#include <cstdio>
#include <algorithm>
extern "C"
{
#include "memory.h"
#include "oddhs.h"
}
#include "osi_cgraph.h"

#define MIN_VIOL 0.02

using namespace std;

CglEClique::CglEClique() :
   colNames(NULL),
   _cgraph(NULL),
   genOddHoles(true),
   argc(0),
   argv(NULL)
{
}

CglEClique::CglEClique(const CglEClique& rhs)
{
   this->_cgraph = rhs._cgraph;
   this->genOddHoles = rhs.genOddHoles;
   this->colNames = rhs.colNames;
}

CglCutGenerator* CglEClique::clone() const
{
   CglEClique *cglec = new CglEClique();

   cglec->_cgraph = this->_cgraph;
   cglec->genOddHoles = this->genOddHoles;
   cglec->colNames = this->colNames;

   return static_cast<CglCutGenerator*>(cglec);
}


void CglEClique::generateCuts( const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info )
{
   //printf("starting eclique generation ... "); fflush(stdout); fflush(stderr);
   //clock_t start = clock();

   const CGraph *cgraph;

   if (_cgraph)
      cgraph = _cgraph;
   else
      cgraph = osi_build_cgraph( (void*) &si );

   CliqueSeparation *sep = clq_sep_create( cgraph );

   if (argc>0)
      clq_sep_set_params_parse_cmd_line( sep, argc, (const char **)argv );

   //clq_sep_set_verbose( sep, 0 );
   clq_sep_set_rc( sep, si.getReducedCost() );
   clq_sep_separate( sep, si.getColSolution() );

   const CliqueSet *clqSet = clq_sep_get_cliques( sep );

   double *ones = (double*) xmalloc( sizeof(double)*si.getNumCols() );
   fill( ones, ones+si.getNumCols(), 1.0);

   /* adding cuts */
   int i;
   for ( i=0 ; (i<clq_set_number_of_cliques(clqSet)) ; ++i )
   {
      const int size = clq_set_clique_size( clqSet, i );
      const int *el = clq_set_clique_elements( clqSet, i );

      OsiRowCut osrc;

      if ( clq_sep_get_verbose( sep ) >= 2 )
      {
         double lhs = 0.0;
         printf("cut:\n");
         for (int i=0 ; (i<size) ; i++)
         {
           printf("(%d %s %g) ", el[i], si.getColName(i).c_str(), ones[i]);
           lhs += si.getColSolution()[el[i]];
         }
         printf("\nlhs: %.6f minViol: %.6f\n\n", lhs, clq_sep_get_min_viol(sep) );
      }

      /* checking violation */
      {
         double lhs = 0.0;
         for (int i=0 ; (i<size) ; i++)
           lhs += si.getColSolution()[el[i]];
         if ((lhs-1.0)<clq_sep_get_min_viol(sep))
            continue;
      }

      osrc.setRow( size, el, &(ones[0]) );
      //osrc.setLb( -si.getInfinity() );
      osrc.setUb( 1.0 );
      osrc.setGloballyValid();

      CoinAbsFltEq equal(1.0e-12);
      cs.insertIfNotDuplicate(osrc,equal);
   }

   /* searching for odd holes */
   if (genOddHoles)
   {
      OddHoleSep *oddhs = oddhs_create();
      oddhs_search_odd_holes( oddhs, si.getNumCols(), si.getColSolution(), si.getReducedCost(), cgraph );

      /* adding odd holes */
      for ( int i=0 ; (i<oddhs_get_odd_hole_count(oddhs)) ; ++i )
      {
         const int *oddEl = oddhs_get_odd_hole( oddhs, i );
         const int oddSize = oddhs_get_odd_hole( oddhs, i+1 ) - oddEl;
         double viol = oddhs_viol( oddSize, oddEl, si.getColSolution() );
         if ( viol < MIN_VIOL )
            continue;

         //printf("\nODD HOLE.\n");

         const int centerSize = oddhs_get_nwc_doh( oddhs, i );
         const int *centerIdx = oddhs_get_wc_doh( oddhs, i );

         const int cutSize = oddSize+centerSize;
         vector< int > idx; idx.reserve( cutSize );
         vector< double > coef; coef.reserve( cutSize );

         idx.insert( idx.end(), oddEl, oddEl+oddSize );
         idx.insert( idx.end(), centerIdx, centerIdx+centerSize );

         coef.insert( coef.end(), ones, ones+oddSize );

         const double rhs = oddhs_rhs( oddSize );

         if ( centerSize )
         {
            vector< double > wcCoefs( centerSize, rhs );
            coef.insert( coef.end(), wcCoefs.begin(), wcCoefs.end() );
         }

/*
         if (colNames)
         {
            for ( int i=0 ; (i<coef.size()) ; ++i )
               printf("%g %s   ", coef[i], (*colNames)[idx[i]].c_str() );
            printf("\n");
            for ( int i=0 ; (i<coef.size()) ; ++i )
               printf("%g  ", si.getColSolution()[idx[i]] );

            printf("\n");
         } */

         OsiRowCut orc;
         orc.setUb( rhs );
         orc.setGloballyValid();
         orc.setRow( cutSize, &(idx[0]), &(coef[0]) );
         cs.insertIfNotDuplicate(orc, CoinAbsFltEq (1.0e-12) );
      }

      oddhs_free( &oddhs );
   }

   free( ones );
   clq_sep_free( &sep );

   //clock_t end = clock();
   //double seconds = ((double)end-start) / ((double)CLOCKS_PER_SEC);
   //printf("done in %g seconds.\n", seconds ); fflush(stdout); fflush(stderr);
}

void CglEClique::setCGraph( const CGraph *cg )
{
   this->_cgraph = cg;
}

void CglEClique::setGenOddHoles( const bool gen )
{
   genOddHoles = gen;
}

void CglEClique::parseParameters( int _argc, const char **_argv )
{
   argc = _argc;
   argv = (char **) _argv;
}

CglEClique::~CglEClique()
{

}

