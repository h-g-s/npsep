#ifndef MPINDEX_H
#define MPINDEX_H

/**
 * Class to store and index mathematical
 * programming object names
 **/

#include <map>
#include <set>
#include <vector>
#include <utility>
#include "mpName.hpp"

typedef std::vector < std::pair<const MPName*, int> > MPIndexedResult;

class MPIndex
{
public:
   MPIndex() :
   delimiter('_') {}

   void associate( const MPName &mpName, const int index );
   void associate( const char *name, const int index );

   /**
    * returns the index of a given name
    * or a negative value if
    * not found
    **/
   int indexOf( const MPName &mpName ) const;
   int indexOf( const char *name ) const;

   /**
    * returns the name relatex to an
    * index of NULL if not found
    **/
   const MPName *nameOf( const int index ) const;

   /**
    * returns all indexes matching
    * a given pattern
    **/
   std::vector< int > matching( const MPName &mpName ) const;
   void matching( const MPName &mpName, MPIndexedResult &vmatches) const;

   size_t size() const;

   void print();

private:
   // returns true if query can be filtered using set intersection
   bool involvedSets( const MPName &mpNF, std::vector< const std::set< int >* > &sets ) const;
   void processIntersection( std::vector< const std::set< int >* > &sets, std::vector< int > &result )  const;

   std::map< MPName, int > mpNameIdx;
   std::map< int, MPName* > mpIdxName;

   // for each var,
   //   in each index position
   //      each value
   //          has a set of associated variable indexes
   std::map< std::string, std::vector< std::map< int, std::set< int > > > > catalog;
   const std::set< int > * getSetFor( const std::string &prefix, int dim, int value ) const;

   char delimiter;
};

#endif // MPINDEX_H
