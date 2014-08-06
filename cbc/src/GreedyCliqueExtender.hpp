#ifndef GREEDYCLIQUEEXTENDER_H
#define GREEDYCLIQUEEXTENDER_H

extern "C"
{
#include "cgraph.h"
}
#include <set>

class GreedyCliqueExtender
{
   public:
      static int extend( CGraph *cgraph, std::set< int > *nodes );
   private:
};

#endif // GREEDYCLIQUEEXTENDER_H
