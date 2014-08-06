#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <exception>
#include <limits>
#include <iostream>
#include <cstdarg>
#include "mpName.hpp"

using namespace std;

MPName::MPName( const char *name, const char delimiter )
{
   // counting how many indexes we have
   int seps = 0;
   register const char *s = name-1;
   while ( *(++s) != '\0' )
      if (*s == delimiter)
	 seps++;

   const int nIdxs = seps;
   indexes.reserve(nIdxs);

   char strDelimiter[2];
   sprintf(strDelimiter, "%c", delimiter);

   char tempName[32];
   strncpy( tempName, name, 32 );

   // separating pattern in tokens
   char *str = NULL, *last;
   str=strtok_r(tempName, strDelimiter, &last );
   prefix = string(str);

   while ( (s=strtok_r(NULL, strDelimiter, &last )) )
   {
      if (strlen(s)==0)
      {
         cerr << "Error at: " << __FILE__ << ":" << __LINE__ << " delimiter separating empty index in MP pattern: " << name << endl;
         throw exception();
      }
      if ( s[0] == '*' )
         indexes.push_back( MPAnyIdx );
      else
         indexes.push_back( atoi(s) );
   } // token separation
}

MPName::MPName( const std::string name, const char delimiter )
{
   MPName( name.c_str(), delimiter);
}

MPName::MPName( const char *_prefix, const int _indexes[], const int count )
{
   prefix = string(_prefix);
   indexes.reserve(count);
   for ( int i=0 ; (i<count) ; i++ )
      indexes.push_back( _indexes[i] );
}

bool MPName::matches(const MPName &other) const
{
   return ( (prefix == other.getPrefix()) &&
            indexesMatch(indexes, other.getIndexes()) );
}

bool MPName::indexesMatch( const std::vector<int> &leftIdxs, const std::vector<int> &rightIdxs ) const
{
   if (leftIdxs.size() != rightIdxs.size())
      return false;

   register std::vector<int>::const_iterator
      vItL = leftIdxs.begin(),
      vItR = rightIdxs.begin();
   for ( ; (vItL!=leftIdxs.end()) ; vItL++, vItR++ )
   {
      if ( (*vItL)!=(*vItR) &&
           ((*vItL)!=MPAnyIdx) && ((*vItR)!=MPAnyIdx) )
           break;
   }

   return ( vItL == leftIdxs.end() );
}

void MPName::operator = ( const MPName &other )
{
   prefix = other.getPrefix();
   indexes = other.getIndexes();
}

bool MPName::operator < ( const MPName &other ) const
{
   if ( prefix != other.getPrefix() )
      return (prefix < other.getPrefix());

   return ( indexes < other.getIndexes() );
}

MPName MPName::create( const char *name, ... )
{
   vector< int > indexes; indexes.reserve( 4 );

   va_list arguments;
   va_start ( arguments, name );
   int v;

   while ( (v=va_arg(arguments,int))!=-1 )
      indexes.push_back( v );

   va_end ( arguments );

   return MPName( name, &(indexes[0]), indexes.size() );
}

string MPName::getName() const
{
   char result[256];
   sprintf(result, "%s", prefix.c_str());
   vector<int>::const_iterator it;
   for ( it=indexes.begin() ; (it!=indexes.end()) ; it++ )
   {
      char strIdx[20];
      sprintf(strIdx, "_%d", *it);
      strcat(result, strIdx);
   }

   return string( result );
}
