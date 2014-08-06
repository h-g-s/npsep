#include "MSepEvent.hpp"
#include "coin/CbcModel.hpp"
#include "GreedyCliqueExtender.hpp"

MSepEvent::MSepEvent( CbcModel *_model, CGraph *_cgraph, set< set<int> > *_cliques, const int *_origCols )
  : CbcEventHandler( _model ),
    model( _model ), cgraph( _cgraph ), cliques( _cliques ), origCols(_origCols)
{
}


CbcEventHandler::CbcAction MSepEvent::event(CbcEvent whichEvent)
{
#ifdef DEBUG
   assert( model!=NULL );
   assert( cgraph!=NULL );
#endif
   switch (whichEvent)
   {
      case solution :
      case heuristicSolution :
         {
            set< int > clique;
            printf("\n -> -> solution found: cost %g number of saved solutions: %d \n", model->getObjValue(), model->numberSavedSolutions() );
            for ( int i=0 ; (i<model->getNumCols()) ; ++i )
               if ( model->bestSolution()[i] >= 0.99 )
                  clique.insert( origCols[i] );
            printf("trying to extend clique with %zu elements ... ", clique.size() );
            GreedyCliqueExtender::extend( cgraph, &clique );
            printf("clique now has %zu elements.\n", clique.size() );
            cliques->insert( clique );
         }

         break;
      default:
         break;
   };

   return noAction;
}

CbcEventHandler *MSepEvent::clone() const
{
   return new MSepEvent( this->model, this->cgraph, this->cliques, this->origCols );
}

MSepEvent::~MSepEvent()
{

}
