#ifndef _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGGRAPH_H_
#define _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGGRAPH_H_

// SOCLE
#include <ign/geometry/graph/GeometryGraph.h>

// APP
#include <app/calcul/detail/graph/EdgeCleaningEdge.h>
#include <app/calcul/detail/graph/EdgeCleaningVertex.h>


namespace app{
namespace calcul{
namespace detail{
namespace graph{

	/// @brief Graph de travail
	class EdgeCleaningGraph : public ign::geometry::graph::GeometryGraph< app::calcul::detail::graph::EdgeCleaningVertex, app::calcul::detail::graph::EdgeCleaningEdge > {
    public:
        /// \brief Constructeur
		EdgeCleaningGraph() {};

		/// \brief Destructeur
		~EdgeCleaningGraph() {};

		/// @brief 
		/// @return Nom du graph
		virtual std::string getName() const {return "EdgeCleaningGraph" ; }
    };
}
}
}
}

#endif