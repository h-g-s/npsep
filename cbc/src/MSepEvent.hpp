#ifndef MSEPEVENT_HPP
#define MSEPEVENT_HPP

#include "cgraph.h"
#include "coin/CbcEventHandler.hpp"
#include <set>

using namespace std;

class MSepEvent : public CbcEventHandler
{
   public:
      /// _cliques points to the clique set which will be updated
      /// origCols must be used if this is a preprocessed model
      MSepEvent( CbcModel *_model, CGraph *_cgraph, set< set<int> > *_cliques, const int *_origCols );

      virtual CbcEventHandler::CbcAction event(CbcEvent whichEvent);

      virtual CbcEventHandler * clone() const;

      virtual ~MSepEvent();

   protected:
      CbcModel *model;
      CGraph *cgraph;
      set< set<int> > *cliques;
      const int *origCols;
};

#endif // MSEPEVENT_HPP
