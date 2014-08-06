#include <cstring>
#include <cmath>
#include <limits>
#include <cfloat>
#include <algorithm>
#include <climits>
#include <OsiSolverInterface.hpp>
#include <CoinPackedMatrix.hpp>
#include "osi_cgraph.h"
extern "C"
{
#include "memory.h"
#include "node_heap.h"
}

#define oo  (INT_MAX/2)
#define WORST_PRIORITY_ROW 10000   // to be considered

using namespace std;

#define EPS 1e-6

/* mininum size for a row to be considered a clique row */
#define MIN_CLIQUE_ROW 250

#define IS_NEGATIVE_DBL( v ) ( v < -1e-5 )

#define DBL_EQUAL( v1, v2 ) ( fabs(v1-v2) < EPS )

// conflict vector dimensions
#define CVEC_CAP   1800000
#define CVEC_FLUSH 1700000

#define MAX_NEGATIVE_COEFS 200

const double LARGE_CONST = std::min( DBL_MAX/10.0, 1e20 );

struct sort_sec_pair {
   bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
      if ( left.second != right.second )
         return left.second < right.second;

      return left.first < right.first;
   }
};

void processClique( const int n, const int *idx, CGraph *cgraph, vector< pair<int, int> > &cvec, const double *colLb, const double *colUb );

/* if size is large enough fech conflicts */
void fetchConflicts( vector< pair<int,int> > &cvec, const bool lastTime, CGraph *cgraph, vector<int> &neighs );

/* returns how much a variable can negatively contribute */
double mostNegativeContribution( const double coef, const char c, const double colLb, const double colUb );

/* returns how much a variable can positively contribute */
inline double unitaryContribution( const double coef, const char c, const double colLb, const double colUb );

/* skip rows with large values to avoid overflows */
bool canYeldLargeValue( const double coef, const double colLb, const double colUb );

int compute_row_priority( OsiSolverInterface *lp, const int r );

CGraph *osi_build_cgraph( void *_lp )
{
   clock_t start = clock();
   OsiSolverInterface *lp = (OsiSolverInterface *)_lp;

   if (lp->getNumIntegers()<2)
      return 0;
   
   int cgraphSize = lp->getNumCols();
#ifdef CONSIDER_BINARY_COMPLEMENT
   cgraphSize *= 2;
   printf("cgraph: considering binary complement.\n");
#endif
   CGraph *cgraph = cgraph_create( cgraphSize );

   const char *ctype = lp->columnType();

   const CoinPackedMatrix *M = lp->getMatrixByRow();

   const double *rhs = lp->getRightHandSide();

   const double *colLb = lp->getColLower();
   const double *colUb = lp->getColUpper();

   const char *sense = lp->getRowSense();

   /* conflict vector */
   vector< pair< int, int > > cvec;
   cvec.reserve( CVEC_CAP );

   /* used i fetch conflicts */
   vector< int > neighs; neighs.reserve( 8192 );

   NodeHeap *nhr = nh_create( lp->getNumRows(), oo );

   for ( int i=0 ; (i<lp->getNumRows()) ; ++i )
   {
      int p = compute_row_priority( lp, i );
      nh_update( nhr, i, p );
   }

   int idxRow, costR;
   while ( (costR=nh_remove_first(nhr,&idxRow)) <= WORST_PRIORITY_ROW )
   {
#define FIXED_IN_ZERO( idx ) ( (fabs(colLb[idx])<EPS) && (fabs(colUb[idx])<EPS) )
      const CoinShallowPackedVector &row = M->getVector(idxRow);
      const int *idx = row.getIndices();
      const double *coefs = row.getElements();

      double minCoef = numeric_limits<double>::max();
      double maxCoef = numeric_limits<double>::min();

      int nInts = 0; // number of integer variables
      int nCont = 0; // number of continuous variables in this row
      int nPos  = 0; // variables which may assume only zero or some positive value
      int nNeg  = 0; // variables which may assume negative values

      if ( (row.getNumElements()<2) || (fabs(rhs[idxRow])>=LARGE_CONST) )
         continue;

      if ( sense[idxRow] == 'R' )  // lets not consider ranged constraints by now
      {
         //printf("TODO: CHECK FOR RANGED CONSTRAINT (%s) rhs is %g\n", lp->getRowName(idxRow).c_str(), rhs[idxRow] );
         continue;
      }

      /* checking number of integers and limits for variables */
      int j;
      double mult = 1.0;
      if ( sense[idxRow] == 'G' )
         mult = -1.0;

      double thisRhs = rhs[idxRow] * mult;

      for ( j=0 ; (j<row.getNumElements()) ; ++j )
      {
         const int cidx = idx[j];

         if (FIXED_IN_ZERO(cidx))
            continue;

         if (ctype[cidx])
            nInts++;
         else
            nCont++;

         /* variable which only accepts positive (or zero) value */
         if ( colLb[cidx] >= -EPS )
            ++nPos;

         if ( colLb[cidx] <= -EPS )
            ++nNeg;

         if ( (nNeg>MAX_NEGATIVE_COEFS) || (canYeldLargeValue(coefs[j],colLb[cidx],colUb[cidx])) )
            break;

         const double neg = mostNegativeContribution( coefs[j], sense[idxRow], colLb[cidx], colUb[cidx] );

         minCoef = min( minCoef, coefs[j]*mult );
         maxCoef = max( maxCoef, coefs[j]*mult );

         thisRhs -= neg;
      } /* all columns */

      /* if there are not enough integers or if some condition
       * triggered premature skip for this row*/
      if ( (j<row.getNumElements()) || (nInts<=1) )
         continue;

#ifdef DEBUG_CONF
      int newConflicts = 0;
      vector< pair<int,int> > confS;
      confS.clear();
#endif // debug
      /* special case (to ignore), when all variables are fixed to zero */
      if ( (fabs(thisRhs)<=EPS) && DBL_EQUAL(maxCoef,minCoef)
         && DBL_EQUAL(maxCoef,1.0) )
         continue;

      /* special case: ready to use cliques */
      if ( DBL_EQUAL( minCoef, maxCoef ) &&  DBL_EQUAL( maxCoef, thisRhs ) &&
            (minCoef>EPS) && ((sense[idxRow]=='E') || (sense[idxRow]=='L'))
            && (nInts==row.getNumElements()) && (nPos==row.getNumElements()) )
      {
         processClique( row.getNumElements(), (const int *)idx, cgraph, cvec, colLb, colUb );
      }
      else
      {
         /* normal case, checking pairwise conflicts
          * (these are computed increasing RHS in thisRhs) */
         int j1, j2, cm1 = row.getNumElements()-1;
         for ( j1=0 ; (j1<cm1) ; ++j1 )
         {
            const int cidx1 = idx[j1];

            if ((!ctype[cidx1])||(FIXED_IN_ZERO(cidx1)))
               continue;

            const double pos1 = unitaryContribution( coefs[j1], sense[idxRow], colLb[cidx1], colUb[cidx1] );

            /* variable cannot be activated anyway, this is the case of integer variables 
             * with strictly positive lower bounds, e.g. miplib 2010 30n20b8, variable s4 */
            if ( pos1-EPS >= thisRhs )   
               continue;

            for ( j2=j1+1 ; (j2<row.getNumElements()) ; ++j2 )
            {
               const int cidx2 = idx[j2];

               if ((!ctype[cidx2]) || (FIXED_IN_ZERO(cidx2)))
                  continue;

               const double pos2 = unitaryContribution( coefs[j2], sense[idxRow], colLb[cidx2], colUb[cidx2] );

               if ( pos2-EPS >= thisRhs )   
                  continue;

               if (pos1+pos2>thisRhs+0.001)
               {
                  cvec.push_back( pair<int,int>(cidx1,cidx2) );
#ifdef DEBUG_CONF
                  newConflicts++;
                  confS.push_back( pair<int,int>(cidx1,cidx2) );
#endif
               }
#ifdef CONSIDER_BINARY_COMPLEMENT
               TODO
#endif // CONSIDER_BINARY_COMPLEMENT
            }
            fetchConflicts( cvec, false, cgraph, neighs );
         }
#ifdef DEBUG_CONF
         if (newConflicts)
         {
            printf("%s,%d,%c,%g : ", lp->getRowName(idxRow).c_str(), nNeg, sense[idxRow], rhs[idxRow]);
            for ( int cc=0 ; (cc<confS.size()) ; ++cc )
               printf("(%s,%s) ", lp->getColName(confS[cc].first).c_str(), lp->getColName(confS[cc].second).c_str() );
            printf("\n");
         }
#endif
      }
#undef FIXED_IN_ZERO
   }

   fetchConflicts( cvec, true, cgraph, neighs );

   nh_free( &nhr );

   clock_t end = clock();

   double ftime = ((double(end-start))/((double)CLOCKS_PER_SEC));
   printf("osi_cgraph took %.3f seconds.\n", ftime);

   return cgraph;
}

void processClique( const int n, const int *idx, CGraph *cgraph, vector< pair<int, int> > &cvec, const double *colLb, const double *colUb )
{
   if ( n >= MIN_CLIQUE_ROW )
   {
      vector< int > vidx; vidx.resize(n);

      memcpy( &(vidx[0]), idx, sizeof(int)*n );

      cgraph_add_clique( cgraph, &(vidx[0]), n );
   }
   else
   {
      int i1, i2, nm1 = n-1;
      for ( i1=0 ; (i1<nm1) ; ++i1 )
      {
         if ( fabs( colUb[idx[i1]] <= EPS ) && fabs( colLb[idx[i1]] <= EPS ) )
            continue;
         for ( i2=i1+1 ; (i2<n) ; ++i2 )
         {
            if ( fabs( colUb[idx[i2]] <= EPS ) && fabs( colLb[idx[i2]] <= EPS ) )
               continue;
            cvec.push_back( pair<int,int>(idx[i1],idx[i2]) );
         }
      }
   }
}

void fetchConflicts( vector< pair<int,int> > &cvec, const bool lastTime, CGraph *cgraph, vector<int> &neighs )
{
   if ( (!lastTime) && (cvec.size()<CVEC_FLUSH) )
      return;

   if ( cvec.size() == 0 )
      return;
   pair< int, int > last;
   int currNode;
   vector< pair<int,int> >::const_iterator vIt;

   cvec.push_back( pair<int,int>( INT_MAX, INT_MAX ) );

   neighs.clear();
   sort( cvec.begin(), cvec.end() );
   currNode = cvec.begin()->first;
   last = pair< int,int >( -1, -1 );

   for ( vIt=cvec.begin() ; (vIt!=cvec.end()) ; last=*vIt,++vIt )
   {
      if (last == *vIt)  // skipping repetitions
         continue;

      if ( vIt->first!=currNode )
      {
         if (neighs.size())
         {
            cgraph_add_node_conflicts_no_sim( cgraph, currNode, &(neighs[0]), neighs.size() );
            neighs.clear();
         }
         currNode = vIt->first;
      }
      neighs.push_back( vIt->second );
   }

   neighs.clear();
   sort( cvec.begin(), cvec.end(), sort_sec_pair() );
   currNode = cvec.begin()->second;
   last = pair< int,int >( -1, -1 );
   for ( vIt=cvec.begin() ; (vIt!=cvec.end()) ; last=*vIt,++vIt )
   {
      if (last == *vIt)  // skipping repetitions
         continue;

      if ( vIt->second!=currNode )
      {
         if (neighs.size())
         {
            cgraph_add_node_conflicts_no_sim( cgraph, currNode, &(neighs[0]), neighs.size() );
            neighs.clear();
         }
         currNode = vIt->second;
      }
      neighs.push_back( vIt->first );
   }

   cvec.clear();
}


double mostNegativeContribution( const double coef, const char s, const double colLb, const double colUb )
{
   double mult = 1.0;

   if (s=='G')
      mult = -1.0;

   double c = coef*mult;

   return min( c*colLb, c*colUb );
}

double unitaryContribution( const double coef, const char s, const double colLb, const double colUb )
{
   if ( colUb <= EPS )
      return 0.0;

   if (s=='G')
      return coef*-1.0;
   else
      return coef;

   // some dumb compilers print warnings if
   // function does not have a return outside ifs
   return coef;
}

bool canYeldLargeValue( const double coef, const double colLb, const double colUb )
{
   if ( fabs( coef ) >= LARGE_CONST )
      return true;

   if ( fabs( colLb ) >= LARGE_CONST )
      return true;

   if ( fabs( colUb ) >= LARGE_CONST )
      return true;

   if ( fabs(colUb)*fabs(coef) >= LARGE_CONST )
      return true;

   if ( fabs(colLb)*fabs(coef) >= LARGE_CONST )
      return true;

   return false;
}

int compute_row_priority( OsiSolverInterface *lp, const int r )
{
   int nInteger = 0;
   int nContinuous = 0;

   const CoinPackedMatrix *M = lp->getMatrixByRow();
   const CoinShallowPackedVector &row = M->getVector(r);
   const int *idx = row.getIndices();
   const double *coefs = row.getElements();
   int nPositive = 0;
   int nNegative = 0;
   const char sense = lp->getRowSense()[r];
   const double rhs = lp->getRightHandSide()[r];
   double minCoef = DBL_MAX;
   double maxCoef = -DBL_MAX;

   const int n = row.getNumElements();
   bool allSameSense;
   int result = 0, i;

   if ( n == 1 )
   {
      result = 10000;
      goto END;
   }

   for ( i = 0; (i<n) ; i++)
   {
      if (lp->isInteger(idx[i]))
         nInteger++;
      else
         nContinuous++;

      if (coefs[i]>EPS)
         nPositive++;
      else
         nNegative++;

      minCoef = std::min( minCoef, coefs[i] );
      maxCoef = std::max( maxCoef, coefs[i] );
   }

   allSameSense = ((!nPositive)||(!nNegative));

   /* checking if this is a "perfect" set packing constraint  */
   if ((
       ( (sense=='L')||(sense=='E') ) &&
       ( nInteger==row.getNumElements() ) &&
       ( allSameSense ) &&
       ( nPositive ) &&
       ( DBL_EQUAL( minCoef, maxCoef) ) &&
       ( DBL_EQUAL( rhs, maxCoef) )
      ))
      result += std::max( (1000 - n), 10 );
   else
      result += 5000;

   /* definitely a constraint were conflicts cannot be extracted */
   if ( (!nNegative) && (sense=='G') )
      result += 6000;

   if (nContinuous)
   {
      if (nContinuous < 10)
         result += 100*nContinuous;
      else
      {
         result = oo;
         goto END;
      }
   }

   if (!allSameSense)
   {
      result += 100;

      /* probably a flow conservation constraint if
         there are many coefficients with both signs */
      if ( (nPositive > 1) and (nNegative > 1) )
         result += min( nPositive, nNegative ) * 1000;   // if both are high values, penalizing

      if ( (nPositive>5) && (nNegative>5) )
      {
         result = oo/2;
         goto END;
      }
   }

END:
   return result;
}
