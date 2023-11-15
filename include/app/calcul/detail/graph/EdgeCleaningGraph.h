#ifndef _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGGRAPH_H_
#define _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGGRAPH_H_


//SOCLE
#include <ign/geometry/graph/GeometryGraph.h>

//APP
#include <app/calcul/detail/graph/EdgeCleaningEdge.h>

namespace app{
namespace calcul{
namespace detail{
namespace graph{

	class EdgeCleaningGraph : public ign::geometry::graph::GeometryGraph< ign::geometry::graph::PunctualVertexProperties, app::calcul::detail::graph::EdgeCleaningEdge > {
    public:
        /// \brief constructor 
		EdgeCleaningGraph() {};

		/// \brief destructor
		~EdgeCleaningGraph() {};

		///\brief return graph Name
		virtual std::string getName() const {return "EdgeCleaningGraph" ; }
    };
}
}
}
}

#endif