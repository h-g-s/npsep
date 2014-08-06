#ifndef MPNAME_H
#define MPNAME_H

/**
 * Class used to store the structured
 * identification of a mathematical
 * programming object in the form:
 *    prefix_index1_index2_...indexn
 * masks can be used to query for a given
 * a MPName object can have "wildcards":
 * in the string this value is specified by "*"
 * which indicates that any index
 * matches in this position
 **/

#include <string>
#include <vector>
#include <climits>

const int MPAnyIdx = INT_MIN;
#ifdef MPNAME
#error MPNAME definition previously defined. Used in mpName.h
#endif

#define MPNAME( name, indexes ) MPName( name, indexes, sizeof(indexes)/sizeof(int) )

class MPName
{
public:
   MPName( const char *name, const char delimiter );
   MPName( const std::string name, const char delimiter );
   MPName( const char *_prefix, const int _indexes[], const int count );

   const std::string &getPrefix() const { return prefix; }

   const std::vector<int> &getIndexes() const { return indexes; }

   bool matches(const MPName &other) const;

   std::string getName() const;

   static MPName create( const char *name, ... );

   void operator = ( const MPName &other );
   bool operator < ( const MPName &other ) const;
private:
   std::string prefix;

   std::vector<int> indexes;

   inline bool indexesMatch( const std::vector<int> &leftIdxs, const std::vector<int> &rightIdxs ) const;
};

#endif
