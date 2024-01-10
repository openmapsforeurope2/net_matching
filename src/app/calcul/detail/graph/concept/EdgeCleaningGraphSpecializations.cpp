// APP
#include <app/calcul/detail/graph/concept/EdgeCleaningGraphSpecializations.h>


namespace ign{
namespace graph{
namespace concept{

    //specialisation pour le concept direction
    template <>
    AllowedDirection direction( app::calcul::detail::graph::EdgeCleaningEdge const& properties )
    {
        return ALLOWED_BOTH;
    }

}
}
}

namespace epg{
namespace graph{
namespace concept{

    //specialisation du concept geometry
	template<>
	ign::geometry::Point vertexGeom( app::calcul::detail::graph::EdgeCleaningGraph const& graph, app::calcul::detail::graph::EdgeCleaningGraph::vertex_descriptor d )
	{	
		return graph.getGeometry(d) ;	
	}

    template<>
	ign::geometry::LineString edgeGeom( app::calcul::detail::graph::EdgeCleaningGraph const& graph, app::calcul::detail::graph::EdgeCleaningGraph::edge_descriptor d )
	{	
		return graph.getGeometry(d) ;	
	}
    
}
}
}

