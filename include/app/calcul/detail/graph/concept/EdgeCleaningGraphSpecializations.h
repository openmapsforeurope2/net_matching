#ifndef _APP_CALCUL_DETAIL_GRAPH_CONCEPT_EDGECLEANINGGRAPHSPECIALIZATIONS_H_
#define _APP_CALCUL_DETAIL_GRAPH_CONCEPT_EDGECLEANINGGRAPHSPECIALIZATIONS_H_

// SOCLE
#include <ign/graph/concept/direction.h>

// EPG
#include <epg/graph/concept/geom.h>

// APP
#include <app/calcul/detail/graph/EdgeCleaningEdge.h>
#include <app/calcul/detail/graph/EdgeCleaningGraph.h>


namespace ign{
namespace graph{
namespace concept{

    //--
    template <>
    AllowedDirection direction( app::calcul::detail::graph::EdgeCleaningEdge const& properties );

}
}
}

namespace epg{
namespace graph{
namespace concept{

    //--
    template <>
    ign::geometry::Point vertexGeom( app::calcul::detail::graph::EdgeCleaningGraph const& graph, app::calcul::detail::graph::EdgeCleaningGraph::vertex_descriptor d  );

    //--
    template <>
    ign::geometry::LineString edgeGeom( app::calcul::detail::graph::EdgeCleaningGraph const& graph, app::calcul::detail::graph::EdgeCleaningGraph::edge_descriptor d  );

}
}
}

#endif
