#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "clique.h"
#include "vint_set.h"
#include "vectormgm.h"
#include "memory.h"
#include "macros.h"
#include "strUtils.h"

#define INI_CAP 1024

#define HASH_SIZE 8192

static const size_t hash_numbers[] = { 37, 31, 29, 17, 13, 11, 7, 1 };
static const size_t n_hash_numbers = 8;

size_t int_vector_hash_code( const int n, const int idx[] );

int int_vector_equals( const int n1, const int idx1[], const int n2, const int idx2[] );

/***
 * returns 1 if clique was already inserted,
 * 0 otherwise
 ***/
int clq_set_clique_already_inserted ( const CliqueSet *clqSet, const int n, const int idx[]  );

struct _CliqueSet
{
   IntSet tmpClique;
   IntSet tmpCliqueOI;
   IntSet *cliques;
   int *W;

   /* indicates the position of each clique */
   IntSet *hash;

   int cliquesCap;
   int numberOfCliques;
   int weightSum;
};

CliqueSet *clq_set_load( const char *fileName )
{
#define LINE_SIZE 65536
#define INT_CHARS 10

   char *line = xmalloc( LINE_SIZE*INT_CHARS );
   int *elements = xmalloc( LINE_SIZE*sizeof(int) );
   int size = 0;
   int weight = 0;
   char **strVector;
   CREATE_STRING_VECTOR( strVector, LINE_SIZE, INT_CHARS );

   CliqueSet *clqSet = clq_set_create();

   FILE *f = fopen( fileName, "r" );
   if (!f)
   {
      fprintf( stderr, "Could not open file %s", &(fileName[0]) );
      exit( EXIT_FAILURE );
   }

   while (fgets(line, LINE_SIZE, f))
   {
      if ( line[strlen(line)-1] == '\n' )
         line[strlen(line)-1] = '\0';

      if ( line[0] == '[' )
      {
         if ( size > 0 )
         {
            clq_set_add( clqSet, size, elements, weight );
            size = 0;
            weight = 0;
         }

         sscanf( line+1, "%d", &weight );
         char *start;
         if ( (start=strstr( line+1, "]" )) == NULL )
         {
            fprintf( stderr, "Invalid file format.\n" );
            exit( EXIT_FAILURE );
         }
         size = splitString( strVector, start+1, ' ', LINE_SIZE, INT_CHARS, 1 );
         {
            int i, nextEl = 0;
            for ( i=0 ; (i<size) ; ++i )
            {
               if ( (strlen(strVector[i])>0) && (isdigit(strVector[i][0])) )
                  elements[nextEl++] = atoi( strVector[i] ) - 1;
            }
            size = nextEl;
         }
      }
      else
      {
         if ( (strlen(line)==0) || (digitsInLine(line, LINE_SIZE)==0) )
            continue;
         if ( size == 0 )
         {
            fprintf( stderr, "Invalid file format.\n" );
            exit( EXIT_FAILURE );
         }
         int newElements = splitString( strVector, line, ' ', LINE_SIZE, INT_CHARS, 1 );
         {
            int i;
            for ( i=0 ; (i<newElements) ; ++i )
               elements[size+i] = atoi( strVector[i] );
         }
         size += newElements;
      }
   }
   if ( size > 0 )
      clq_set_add( clqSet, size, elements, weight );

   FREE_STRING_VECTOR( strVector );
   fclose( f );
   free( line );
   free( elements );


   return clqSet;
#undef LINE_SIZE
#undef INT_CHARS
}


CliqueSet *clq_set_create( )
{
   int i;

   CliqueSet *clqSet = xmalloc( sizeof(CliqueSet) );

   clqSet->hash = xmalloc( sizeof(IntSet)*HASH_SIZE );
   for ( i=0 ; (i<HASH_SIZE) ; ++i )
      vint_set_init( clqSet->hash+i );

   clqSet->numberOfCliques = 0;
   clqSet->cliquesCap = INI_CAP;
   clqSet->cliques = xmalloc( sizeof(IntSet)*clqSet->cliquesCap );
   for ( i=0 ; (i<clqSet->cliquesCap) ; ++i )
      vint_set_init( clqSet->cliques+i );

   clqSet->W = xmalloc( sizeof(int)*clqSet->cliquesCap ) ;
   clqSet->weightSum = 0;
   vint_set_init( &(clqSet->tmpClique) );
   vint_set_init( &(clqSet->tmpCliqueOI) );

   return clqSet;
}

int clq_set_add( CliqueSet *clqSet, const int size, const int nodes[], const int w )
{
   vint_set_clear( &(clqSet->tmpClique) );
   vint_set_add( &(clqSet->tmpClique), nodes, size );

   const int *snodes = vint_set_get_elements( &(clqSet->tmpClique) );

   if ( clq_set_clique_already_inserted( clqSet, size, snodes ) )
      return 0;

   if ( clqSet->cliquesCap < clqSet->numberOfCliques+1 )
   {
      int i;
      const int oldCap = clqSet->cliquesCap;
      clqSet->cliquesCap *= 2;
      clqSet->W = xrealloc( clqSet->W, sizeof(int)*clqSet->cliquesCap );
      clqSet->cliques = xrealloc( clqSet->cliques, sizeof(IntSet)*clqSet->cliquesCap );
      for ( i=oldCap ; (i<clqSet->cliquesCap) ; ++i )
         vint_set_init( clqSet->cliques+i );
   }

   vint_set_add( clqSet->cliques+clqSet->numberOfCliques, snodes, size );
   clqSet->W[clqSet->numberOfCliques] = w;
   clqSet->weightSum += w;

   /* inserting into hash table */
   size_t hash_code = int_vector_hash_code( size, snodes );
#ifdef DEBUG
   assert( hash_code<HASH_SIZE );
#endif

   vint_set_add( clqSet->hash + hash_code, &(clqSet->numberOfCliques)  , 1 );

   clqSet->numberOfCliques++;

   return 1;
}

int clq_set_weight( const CliqueSet *clqSet, const int clique )
{
   return clqSet->W[clique];
}


int clq_set_clique_size( const CliqueSet *clqSet, const int clique )
{
   return ( vint_set_size( clqSet->cliques+clique ) );
}

const int *clq_set_clique_elements( const CliqueSet *clqSet, const int clique )
{
   return vint_set_get_elements( clqSet->cliques+clique );
}

void clq_set_free( CliqueSet **clqSet )
{
   int i;
   for ( i=0 ; (i<(*clqSet)->cliquesCap) ; ++i )
      vint_set_clean( (*clqSet)->cliques+i );

   for ( i=0 ; (i<HASH_SIZE) ; ++i )
      vint_set_clean( (*clqSet)->hash+i );

   vint_set_clean( &((*clqSet)->tmpClique) );
   vint_set_clean( &((*clqSet)->tmpCliqueOI) );

   free( (*clqSet)->cliques );
   free( (*clqSet)->hash );
   free( (*clqSet)->W );
   free( (*clqSet) );
   (*clqSet) = NULL;
}


int clq_validate( const CGraph *cgraph, const int size, const int nodes[], int *n1, int *n2 )
{
   int i,j;

   *n1 = -1;
   *n2 = -1;
   const int sizeM1 = size-1;
   for ( i=0 ; (i<sizeM1) ; ++i )
      for ( j=i+1 ; (j<size) ; ++j )
      {
         if ( (!cgraph_conflicting_nodes( cgraph, nodes[i], nodes[j] )) || (nodes[i] == nodes[j]) )
         {
            *n1 = nodes[i];
            *n2 = nodes[j];
            return 0;
         }
      }

   return 1;
}

int clq_comp_int( const void *v1, const void *v2 )
{
   return (*((const int *)v1)) - (*((const int *)v2));
}

int clq_conflicts_with_all( const CGraph *cgraph, const int node, const int size, const int nodes[] )
{
   int i;
   for ( i=0 ; (i<size) ; ++i )
      if (!cgraph_conflicting_nodes( cgraph, node, nodes[i] ))
         return 0;

   return 1;
}

int clq_set_number_of_cliques( const CliqueSet *clqSet )
{
   return clqSet->numberOfCliques;
}

void clq_set_print( const CliqueSet *clqSet )
{
   int i;
   for ( i=0 ; (i<clq_set_number_of_cliques( clqSet)) ; ++i )
   {
      printf("[%d] " , clq_set_weight( clqSet, i) );
      const int *el = clq_set_clique_elements( clqSet, i );
      int j;
      for ( j=0 ; (j<clq_set_clique_size(clqSet, i) ) ; ++j )
         printf("%d ", el[j]+1 );
      printf("\n");
   }
}

int clq_set_weight_sum( CliqueSet *clqSet )
{
   return clqSet->weightSum;
}

size_t int_vector_hash_code( const int n, const int idx[] )
{
   size_t code =  0;

   code += n * hash_numbers[0];

   code += idx[0] * hash_numbers[1];

   size_t i;
   for ( i=1 ; (i<n) ; ++i )
      code += hash_numbers[ i%n_hash_numbers ] * idx[i];

   code = code % HASH_SIZE;

#ifdef DEBUG
   assert( code >=0 ); assert( code < HASH_SIZE );
#endif

   return code;
}

int int_vector_equals( const int n1, const int idx1[], const int n2, const int idx2[] )
{
   if ( n1 != n2 )
      return 0;

   register int i;
   for ( i=0 ; (i<n1) ; ++i )
      if ( idx1[i] != idx2[i] )
      return 0;

   return 1;
}

int clq_set_clique_already_inserted ( const CliqueSet *clqSet, const int n, const int idx[]  )
{
   const size_t hash_code = int_vector_hash_code( n, idx );
   const IntSet *hashEntry = clqSet->hash +  hash_code ;
   const int hashEntrySize = vint_set_size( hashEntry );
   const int *hashEntryElements = vint_set_get_elements( hashEntry );

   int i;
   for ( i=0 ; (i<hashEntrySize) ; ++i )
   {
      const int cliqueIndex = hashEntryElements[i];
      const IntSet *otherClique = clqSet->cliques + cliqueIndex;

      if ( int_vector_equals( n, idx, vint_set_size(otherClique), vint_set_get_elements(otherClique) ) )
         return 1;
   }

   return 0;
}

int clq_set_clique_has_element( const CliqueSet *clqSet, const int clique, const int element )
{
   if (vint_set_find( clqSet->cliques + clique, element ))
      return 1;

   return 0;
}

void clq_set_save( const CGraph *cgraph, const CliqueSet *clqSet, const char *fileName )
{
   FILE *f = fopen( fileName, "w" );

   int i;
   for ( i=0 ; (i<clqSet->numberOfCliques) ; ++i )
   {
      const int size = clq_set_clique_size( clqSet, i );
      const int *elements = clq_set_clique_elements( clqSet, i );
      const int w = clq_set_weight( clqSet, i );
      fprintf( f, "[%d]", w );
      int j;
      for ( j=0 ; (j<size) ; ++j )
      {
         fprintf( f, " %d", elements[j]+1 );
         if (strlen(cgraph_get_node_name( cgraph, elements[j]))>0)
            fprintf( f, "(%s)", cgraph_get_node_name( cgraph, elements[j]) );
         fprintf( f, " " );
      }
      fprintf( f, "\n" );
   }

   fclose( f );
}

void clq_set_clear( CliqueSet *clqSet )
{
   /* clearing hash contents */
   {
      int i;
      for (  i=0 ; (i<clq_set_number_of_cliques( clqSet )) ; ++i )
      {
         const int n = clq_set_clique_size( clqSet, i );
         const int *el = clq_set_clique_elements(  clqSet, i );
         size_t entry = int_vector_hash_code( n, el );
         vint_set_clear( clqSet->hash+entry );
      }
   }
   /* clearing cliques */
   {
      int i;
      for ( i=0 ; (i<clqSet->cliquesCap) ; ++i )
         vint_set_clear( clqSet->cliques+i );
   }

   clqSet->numberOfCliques = 0;
   clqSet->weightSum = 0;
   memset( clqSet->W, 0, sizeof(int)*clqSet->cliquesCap );
}

void clq_set_cpy( CliqueSet *clqs_target, const CliqueSet *clqs_source )
{
   clq_set_clear( clqs_target );
   clq_set_add_cliques( clqs_target, clqs_source );
}

int clq_set_add_cliques( CliqueSet *clqs_target, const CliqueSet *clqs_source )
{
   int result = 0;
   {
      int i;
      for ( i=0 ; ( i<clq_set_number_of_cliques( clqs_source ) ) ; ++i )
      {
         const int *el = clq_set_clique_elements( clqs_source, i );
         const int weight = clq_set_weight( clqs_source, i );

         result += clq_set_add( clqs_target, clq_set_clique_size( clqs_source,i ), el, weight );
      }
   }
   return result;
}

void clq_set_add_using_original_indexes( CliqueSet *target, const CliqueSet *source, const int orig[] )
{
   IntSet *tmps = &(target->tmpCliqueOI);

   {
      int i;
      for ( i=0 ; (i<clq_set_number_of_cliques(source))  ; ++i )
      {
         vint_set_clear( tmps );
         const IntSet *clique = clq_set_get_clique( source, i);
         const int weight = clq_set_weight( source, i );
         const int size = vint_set_size( clique );
         const int *el = vint_set_get_elements( clique );
         vint_set_add_using_original_indexes( tmps, el, size, orig );
         clq_set_add( target, vint_set_size(tmps), vint_set_get_elements(tmps), weight );
      }
   }
}

const IntSet *clq_set_get_clique( const CliqueSet *clqSet, const int idx )
{
   return clqSet->cliques + idx;
}

int clq_dominates( const IntSet *a, const IntSet *b )
{
   const int sizeA = vint_set_size( a );
   const int sizeB = vint_set_size( b );

   if ( sizeA <= sizeB )
      return 0;

   /* a is bigger, checking if it has all elements of b */
   const int *elB = vint_set_get_elements( b );
   int i;
   for ( i=0 ; (i<sizeB) ; ++i )
      if ( vint_set_find( a, elB[i] )==NULL )
         return 0;

   return 1;
}


