#include "GreedyCliqueExtender.hpp"

using namespace std;

int GreedyCliqueExtender::extend( CGraph *cgraph, std::set< int > *nodes )
{
   int result = 0;

   for ( int i=0 ; ( i<cgraph_size(cgraph) ) ; ++i )
   {
      if ( nodes->find( i ) == nodes->end() ) /// not inserted yet
      {
         /// if has conflict with everyone
         set< int >::const_iterator nIt=nodes->begin();
         for (  ; (nIt!=nodes->end()) ; ++nIt )
         {
            if ( !( cgraph_conflicting_nodes( cgraph, i, *nIt ) ) )
               break;
         }
         if ( nIt == nodes->end() )
         {
            nodes->insert( i );
            ++result;
         }
      }
   }

   return result;
}

