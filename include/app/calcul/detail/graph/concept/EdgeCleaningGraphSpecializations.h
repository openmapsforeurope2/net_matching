#ifndef _APP_CALCUL_DETAIL_GRAPH_CONCEPT_EDGECLEANINGGRAPHSPECIALIZATIONS_H_
#define _APP_CALCUL_DETAIL_GRAPH_CONCEPT_EDGECLEANINGGRAPHSPECIALIZATIONS_H_

// SOCLE
#include <ign/graph/concept/direction.h>

// APP
#include <app/calcul/detail/graph/EdgeCleaningEdge.h>



namespace ign{
namespace graph{
namespace concept{


    //specialisation pour le concept direction
    template <>
    AllowedDirection direction( app::calcul::detail::graph::EdgeCleaningEdge const& properties );
}
}
}

#endif
