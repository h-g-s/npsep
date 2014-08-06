#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include "cgraph.h"
#include "memory.h"
#include "clique.h"
#include "vint_set.h"
#include "clique_extender.h"
#include "vectormgm.h"
#include "macros.h"

#define MAX_CLIQUE_SIZE 128

/* additional space for candidates */
#define CANDIDATES_SLACK 100

struct _CliqueExtender
{
   const CGraph *cgraph;
   CliqueSet *clqSet;
   int *candidates;
   int *newClique;
   int newCliqueSize;
   int candidatesCap;

   NeighIterator *nit;

   int costsCap;
   int *costs;
};

void clqe_check_nodes_cap( CliqueExtender *clqe )
{
   const CGraph *cgraph = clqe->cgraph;

   const int nodes = cgraph_size( cgraph );

   if ( clqe->candidatesCap < nodes )
   {
      clqe->candidatesCap = nodes;
      if ( !clqe->candidates )
      {
         clqe->candidates = xmalloc( sizeof(int)*nodes*CANDIDATES_SLACK );
         clqe->newClique = xmalloc( sizeof(int)*nodes );
      }
      else
      {
         clqe->candidates = xrealloc( clqe->candidates, sizeof(int)*clqe->candidatesCap*CANDIDATES_SLACK );
         clqe->newClique = xrealloc( clqe->newClique, sizeof(int)*clqe->candidatesCap );
      }
   }
}

CliqueExtender *clqe_create()
{
   CliqueExtender *clqe = xmalloc( sizeof(CliqueExtender) );

   clqe->cgraph  = NULL;
   clqe->clqSet  = clq_set_create();
   clqe->candidatesCap = 0;
   clqe->candidates = NULL;
   clqe->newClique = NULL;

   clqe->costs = NULL;
   clqe->costsCap = 0;

   clqe->nit = nit_create();

   return clqe;
}

int clqe_get_best_candidates_clique_insertion( CliqueExtender *clqe, const IntSet *clique )
{
   /* node with the smallest degree */
   int nodeSD = -1, degree = INT_MAX, i;

   const int cliqueSize = vint_set_size( clique ), *cliqueEl = vint_set_get_elements( clique );

   const CGraph *cgraph = clqe->cgraph;

   /* picking node with the smallest degree */
   for ( i=0 ; (i<cliqueSize) ; ++i )
   {
      if ( cgraph_degree( cgraph, cliqueEl[i] ) < degree )
      {
         degree = cgraph_degree( cgraph, cliqueEl[i] );
         nodeSD = cliqueEl[i];
      }
   }

   const int *costs = clqe->costs;

#ifdef DEBUG
   assert( costs );
   int previousCost = INT_MIN;
#endif

   NeighIterator *nit = clqe->nit;
   int selected = -1;

   nit_start( nit, cgraph, nodeSD, costs  );

   int nCandidates = 0;
   int *candidates = clqe->candidates;
   while ( ( (selected=nit_next(nit))!=INT_MAX ) && (nCandidates<MAX_CLIQUE_SIZE) )
   {
#ifdef DEBUG
      int curCost = costs[selected];
      assert( curCost>=previousCost );
      previousCost = curCost;
#endif
      /* need to have conflict with all nodes in clique an all others inserted */
      for ( i=0 ; (i<cliqueSize) ; ++i )
         if ( (!cgraph_conflicting_nodes( cgraph, cliqueEl[i], selected )) || (selected==cliqueEl[i]) )
            break;
      if (i<cliqueSize)
         continue;
      for ( i=0 ; (i<nCandidates) ; ++i )
         if (!cgraph_conflicting_nodes( cgraph, candidates[i], selected ))
            break;
      if (i<nCandidates)
         continue;

      candidates[nCandidates++] = selected;
   }

   return nCandidates;
}


int clqe_extend( CliqueExtender *clqe, const CGraph *cgraph, const IntSet *clique,
                 const int weight, const CliqueExtendingMethod clqem )
{
   int result = 0;

#ifdef DEBUG
   assert( (clqe) && (cgraph) && (clique) && ((vint_set_size(clique))) );
   int en1, en2;
   if (!clq_validate( cgraph, vint_set_size(clique), vint_set_get_elements(clique), &en1, &en2 ))
   {
      fprintf( stderr, "ERROR clqe_extend : Nodes %d and %d are not in conflict.\n", en1, en2 );
      exit( EXIT_FAILURE );
   }
#endif


   clqe->cgraph = cgraph;

   clqe_check_nodes_cap( clqe );

   int nCandidates = clqe_get_best_candidates_clique_insertion( clqe, clique );
   if (!nCandidates)
      goto TERMINATE;

   /* too many candidates, filtering */
   if (nCandidates > (MAX_CLIQUE_SIZE*2))
      nCandidates = (MAX_CLIQUE_SIZE*2);

   /* clique can be extended, starting to fill new clique */
   memcpy( clqe->newClique, vint_set_get_elements( clique ), sizeof(int)*vint_set_size(clique) );
   clqe->newCliqueSize = vint_set_size( clique );

   int idxSelected = -1, selectedNode, i, removals;
   int *candidates = clqe->candidates;

INSERT_CANDIDATE:
   if ( ( clqem == CLQEM_RANDOM ) || (clqe->costs == NULL) || (clqe->costsCap<cgraph_size(cgraph)) )
   {
      idxSelected = ((int)((((double)rand()) / (((double)RAND_MAX)+1.0)) * ((double)nCandidates)));
      if ( clqem == CLQEM_PRIORITY_GREEDY )
         fprintf( stderr, "Warning: using random selection for extension since no costs were informed.\n");
   }
   else
   {
      /* costs informed, picking the one with the lowest cost */
      int i, lowestCost = INT_MAX;
      for ( i=0 ; (i<nCandidates) ; ++i )
      {
         if ( clqe->costs[ candidates[i] ] < lowestCost )
         {
            lowestCost = clqe->costs[ candidates[i] ];
            idxSelected = i;
         }
      }
   }

#ifdef DEBUG
   assert( idxSelected >= 0 );
   assert( idxSelected < nCandidates );
#endif

   selectedNode = candidates[ idxSelected ];

   assert( selectedNode>=0 && selectedNode < cgraph_size(cgraph) );

   clqe->newClique[ clqe->newCliqueSize++ ] = selectedNode;

   /* removing from candidates those which do not conflict with new inserted node */
   removals = 0;
   for ( i=0 ; (i<nCandidates); ++i )
   {
      if ( ( !cgraph_conflicting_nodes(cgraph, selectedNode, candidates[i] ) ) || (candidates[i]==selectedNode) )
      {
         candidates[i] = INT_MAX;
         ++removals;
      }
   }

   qsort( candidates, nCandidates, sizeof(int), vint_set_cmp_int );

   nCandidates -= removals;

   if ( ( nCandidates ) && (clqe->newCliqueSize<MAX_CLIQUE_SIZE) )
      goto INSERT_CANDIDATE;

   /* at this point we have an extended clique */
   result += clq_set_add( clqe->clqSet, clqe->newCliqueSize, clqe->newClique, weight );

TERMINATE:

   /*
   {
      int n1, n2;
      int res = clq_validate( cgraph, clqe->newCliqueSize, clqe->newClique, &n1, &n2 );
      assert(res);
      int i;
      printf("CLIQUE BEF ADD SIZE:\n");
      for ( i=0 ; (i<clqe->newCliqueSize) ; ++i )
         printf( "%d ", clqe->newClique[i] );

      int sizeS = clq_set_clique_size( clqe->clqSet, clq_set_number_of_cliques(clqe->clqSet)-1 );
      const int *el = clq_set_clique_elements( clqe->clqSet, clq_set_number_of_cliques(clqe->clqSet)-1 );
      printf("\n-> -> CLIQUE BEF ADD SIZE: <- <-\n");
      for ( i=0 ; (i<sizeS) ; ++i )
         printf( "%d ", el[i] );

      printf("\n");
   } */

   return result;
}

void clqe_set_clear( CliqueExtender *clqe )
{
   clq_set_clear( clqe->clqSet );
}

const CliqueSet *clqe_get_cliques( CliqueExtender *clqe )
{
   return clqe->clqSet;
}

void clqe_set_costs( CliqueExtender *clqe, const int costs[], const int n )
{
   vmg_adjust_vector_capacity( (void**)&(clqe->costs), &(clqe->costsCap), n, sizeof(int) );

   memcpy( clqe->costs, costs, sizeof(int)*n );

   if ( clqe->costsCap > n )
      memset( clqe->costs+n, 0, (clqe->costsCap-n)*sizeof(int) );
}

const int *clqe_get_costs( CliqueExtender *clqe )
{
   return clqe->costs;
}

void clqe_free( CliqueExtender  **clqe )
{
   if ( (*clqe)->newClique )
      free( (*clqe)->newClique );
   if ( (*clqe)->candidates )
      free( (*clqe)->candidates );

   clq_set_free( &(((*clqe)->clqSet) ) );

   if ( (*clqe)->costs )
      free( (*clqe)->costs );

   nit_free( &((*clqe)->nit) );

   free( *clqe );
   (*clqe) = NULL;
}
