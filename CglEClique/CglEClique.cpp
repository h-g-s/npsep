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

void fillColSolution(const OsiSolverInterface &si, double colSol[])
{
   const int numCols = si.getNumCols();
   const double* origColSol = si.getColSolution();

   for(int i = 0; i < numCols; i++)
   {
      colSol[i] = origColSol[i];
      colSol[i+numCols] = 1 - origColSol[i];
   }
}

void fillReducedCost(const OsiSolverInterface &si, double rCost[])
{
   const int numCols = si.getNumCols();
   const double* origRCost = si.getReducedCost();

   for(int i = 0; i < numCols; i++)
   {
      rCost[i] = origRCost[i];
      rCost[i+numCols] = 0.0;
   }
}

void CglEClique::generateCuts( const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info )
{
   //printf("starting eclique generation ... "); fflush(stdout); fflush(stderr);
   //clock_t start = clock();

   const CGraph *cgraph;
   const int numCols = si.getNumCols(), numRows = si.getNumRows();

   if (_cgraph)
      cgraph = _cgraph;
   else
      cgraph = osi_build_cgraph( (void*) &si );

   CliqueSeparation *sep = clq_sep_create( cgraph );

   if (argc>0)
      clq_sep_set_params_parse_cmd_line( sep, argc, (const char **)argv );

   //clq_sep_set_verbose( sep, 0 );

   double colSol[numCols*2], rCost[numCols*2];
   fillColSolution(si, colSol);
   fillReducedCost(si, rCost);
   clq_sep_set_rc( sep, rCost );
   clq_sep_separate( sep, colSol );

   const CliqueSet *clqSet = clq_sep_get_cliques( sep );

   /* adding cuts */
   int i;
   int varCount[numCols];               //counting how many times a variable(and its complement)
   fill(varCount, varCount+numCols, 0); //appears in the current clique
   for ( i=0 ; (i<clq_set_number_of_cliques(clqSet)) ; ++i )
   {
      const int size = clq_set_clique_size( clqSet, i );
      const int *el = clq_set_clique_elements( clqSet, i );
      double coefs[size];
      int maxFrequency = 0;
      int idxs[size];
      OsiRowCut osrc;

      if ( clq_sep_get_verbose( sep ) >= 2 )
      {
         double lhs = 0.0, rhs = 1.0;
         printf("cut:\n");
         for (int j=0 ; (j<size) ; j++)
         {
            if(el[j] < numCols)
            {
               printf("(%d %s %g) ", el[j], (*colNames)[el[j]].c_str(), 1.0);
               lhs += colSol[el[j]];
            }
            else
            {
               printf("(%d %s %g) ", el[j], (*colNames)[el[j]].c_str(), -1.0);
               lhs -= colSol[el[j]-numCols];  
               rhs -= 1.0;
            }
         }
         printf("\nlhs: %.6f rhs: %.6f minViol: %.6f\n\n", lhs, rhs, clq_sep_get_min_viol(sep) );
      }

      double lhs = 0.0, rhs = 1.0;
      /* checking violation */
      {
         for (int j=0 ; (j<size) ; j++)
         {
            if(el[j] < numCols)
            {
               coefs[j] = 1.0;
               lhs += colSol[el[j]];
               idxs[j] = el[j];
               maxFrequency = max(maxFrequency, ++varCount[el[j]]);
            }
            else
            {
               coefs[j] = -1.0;
               lhs -= colSol[el[j]-numCols];
               rhs -= 1.0;
               idxs[j] = el[j]-numCols;
               maxFrequency = max(maxFrequency, ++varCount[el[j]-numCols]);
            }
         }

         assert(maxFrequency >= 0 && maxFrequency <= 2);

         if(maxFrequency == 2)
         {
            for (int j=0 ; (j<size) ; j++)
            {
               assert(varCount[idxs[j]] >= 0 && varCount[idxs[j]] <= 2);
               if(varCount[idxs[j]] == 1)
               {
                  int varIdx[] = {idxs[j]};
                  double varCoef[] = {coefs[j]};
                  osrc.setRow( 1, varIdx, &(varCoef[0]) );
                  //osrc.setLb( -si.getInfinity() );
                  osrc.setUb( 0.0 ); 
                  osrc.setGloballyValid();
                  CoinAbsFltEq equal(1.0e-12);
                  cs.insertIfNotDuplicate(osrc,equal);
               }
               varCount[idxs[j]] = 0;
            }
            continue;
         }

         for (int j=0 ; (j<size) ; j++)
            varCount[idxs[j]] = 0;

         if ((lhs-rhs)<clq_sep_get_min_viol(sep))
            continue;
      }

      osrc.setRow( size, idxs, &(coefs[0]) );
      //osrc.setLb( -si.getInfinity() );
      osrc.setUb( rhs ); 
      osrc.setGloballyValid();
      CoinAbsFltEq equal(1.0e-12);
      cs.insertIfNotDuplicate(osrc,equal);
   }

   /* searching for odd holes */
   if (genOddHoles)
   {
      double *ones = (double*) xmalloc( sizeof(double)*numCols*2 );
      fill( ones, ones+(numCols*2), 1.0);

      OddHoleSep *oddhs = oddhs_create();
      oddhs_search_odd_holes( oddhs, si.getNumCols(), si.getColSolution(), si.getReducedCost(), cgraph );

      /* adding odd holes */
      for ( int j=0 ; (j<oddhs_get_odd_hole_count(oddhs)) ; ++j )
      {
         const int *oddEl = oddhs_get_odd_hole( oddhs, j );
         const int oddSize = oddhs_get_odd_hole( oddhs, j+1 ) - oddEl;
         double viol = oddhs_viol( oddSize, oddEl, si.getColSolution() );
         if ( viol < MIN_VIOL )
            continue;

         //printf("\nODD HOLE.\n");

         const int centerSize = oddhs_get_nwc_doh( oddhs, j );
         const int *centerIdx = oddhs_get_wc_doh( oddhs, j );

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
            for ( int j=0 ; (j<coef.size()) ; ++j )
               printf("%g %s   ", coef[j], (*colNames)[idx[j]].c_str() );
            printf("\n");
            for ( int j=0 ; (j<coef.size()) ; ++j )
               printf("%g  ", si.getColSolution()[idx[j]] );

            printf("\n");
         } */

         OsiRowCut orc;
         orc.setUb( rhs );
         orc.setGloballyValid();
         orc.setRow( cutSize, &(idx[0]), &(coef[0]) );
         cs.insertIfNotDuplicate(orc, CoinAbsFltEq (1.0e-12) );
      }

      oddhs_free( &oddhs );
      free( ones );
   }

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

