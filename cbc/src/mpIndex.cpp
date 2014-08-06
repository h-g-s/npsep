#include "mpIndex.hpp"
#include <climits>
#include <algorithm>
#include <cassert>
#include <cstdio>

using namespace std;

void MPIndex::associate( const MPName &mpName, const int index )
{
   mpNameIdx[ mpName ] = index;
   map< MPName, int >::const_iterator it = mpNameIdx.find( mpName );

   mpIdxName[ index ] = (MPName*) &(it->first) ;

   std::vector< std::map< int, std::set< int > > > &varCatalog = catalog[ mpName.getPrefix() ];

   const size_t nIndexes = mpName.getIndexes().size();
   if ( varCatalog.size() < nIndexes  )
      varCatalog.resize( nIndexes );

   for ( size_t i=0 ; (i<nIndexes) ; ++i )
   {
      const int value = mpName.getIndexes()[i];

      if ( value == MPAnyIdx )
         continue;

      varCatalog[i][value].insert( index );
   }
}

void MPIndex::associate( const char *name, const int index )
{
   associate( MPName(name, delimiter), index );
}

int MPIndex::indexOf( const char *name ) const
{
   return indexOf( MPName(name, delimiter) );
}

vector< int > MPIndex::matching( const MPName &mpName ) const
{
   vector< int > result;
   result.reserve( 500 );

   vector< const set<int> * > invSets;
   involvedSets( mpName, invSets );
   processIntersection( invSets, result );

   return result;
}

void MPIndex::matching( const MPName &mpName, MPIndexedResult &vmatches) const
{
   vmatches.reserve( 500 );

   register std::map<MPName, int >::const_iterator it;
   for ( it=mpNameIdx.begin() ; (it!=mpNameIdx.end()) ; it++ )
   {
      if (mpName.matches( it->first ))
      {
        vmatches.push_back( pair <const MPName*, int > (&(it->first), it->second));
      }
   }
}

const MPName *MPIndex::nameOf( const int index ) const
{
   std::map< int, MPName* >::const_iterator it;
   it = mpIdxName.find( index );
   if (it == mpIdxName.end())
      return 0;

   return it->second;
}

size_t MPIndex::size() const
{
   return mpNameIdx.size();
}

const set< int > * MPIndex::getSetFor( const std::string &prefix, int dim, int value ) const
{
   map< string, vector< map< int, set< int > > > >::const_iterator cIt = catalog.find( prefix );
   if ( cIt == catalog.end() )
      return NULL;

   assert( cIt->second.size() > (size_t)dim );

   const map< int, set< int > > &mapSet = (cIt->second)[dim];
   map< int, set< int > >::const_iterator sIt = mapSet.find( value );
   if ( sIt == mapSet.end() )
      return NULL;

   return ( const set< int > * ) &(sIt->second);
}

bool MPIndex::involvedSets( const MPName &mpNF, vector< const set< int >* > &sets ) const
{
   const vector< int > &idxs = mpNF.getIndexes();
   sets.reserve( idxs.size() );
   vector< int >::const_iterator vIt = idxs.begin();
   int dim = 0;
   for ( ; (vIt != idxs.end()) ; vIt++,++dim )
   {
      const int idx = *vIt;
      if (idx == MPAnyIdx)
         continue;

      const set< int > *cs = getSetFor( mpNF.getPrefix(), dim, idx );
      assert( cs != NULL );

      sets.push_back( cs );
   } // all vector indexes

   return (sets.size()>0);
}

bool SortBySetSize ( const std::set< int >* s1, const std::set< int >* s2 )
{
   return (s1->size()) < (s2->size());
}

void MPIndex::processIntersection( vector< const set< int >* > &sets, vector< int > &result ) const
{
   sort( sets.begin(), sets.end(), SortBySetSize );

#ifdef DEBUG
   if (sets.size()>1)
   {
      vector< const set< int >* >::iterator vIt = sets.begin();

      const set< int > *s1 = (*vIt);
      const set< int > *s2 = (*(++vIt));
      assert( s1->size() <= s2->size() );
   }
#endif

   const set< int > *s = (*(sets.begin()));
   result.reserve( s->size() ); result.clear();

   const vector< const set< int >* >::const_iterator secondSetIt = (sets.begin()+1);

   // interating into the first set
   for ( set< int >::const_iterator sIt=s->begin() ; (sIt!=s->end()) ; ++sIt )
   {
      const int idx = *sIt;

      // it must exist also in all other sets
      vector< const set< int >* >::const_iterator osIt = secondSetIt;
      for ( ; (osIt != sets.end()) ; ++osIt )
      {
         if ( (*osIt)->find( idx ) == (*osIt)->end() )
            break;
      } // checking existence in other sets
      if ( osIt == sets.end() )
         result.push_back( idx );
   } // all elements in the smallest set
}

void MPIndex::print()
{
   map< MPName, int >::const_iterator it;
   for ( it=mpNameIdx.begin() ; (it!=mpNameIdx.end()) ; ++it )
      printf("%s %d\n", it->first.getName().c_str(), it->second );
}

int MPIndex::indexOf( const MPName &mpName )  const
{
   std::map< MPName, int >::const_iterator it;
   it = mpNameIdx.find( mpName );
   if (it == mpNameIdx.end())
      return INT_MIN;

   return it->second;
}
