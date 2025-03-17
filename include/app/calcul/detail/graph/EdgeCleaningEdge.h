#ifndef _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGEDGE_H_
#define _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGEDGE_H_


// SOCLE
#include <ign/geometry/Point.h>
#include <ign/geometry/graph/LinearEdgeProperties.h>



namespace app{
namespace calcul{
namespace detail{
namespace graph{

	/// @brief Propriétés des arcs pour le graph EdgeCleaningGraph
	struct EdgeCleaningEdge : public ign::geometry::graph::LinearEdgeProperties {

        /// \brief constructeur
        EdgeCleaningEdge(){};

        /// \brief destructeur
        ~EdgeCleaningEdge(){};

        //--
        double        weight;
    };
}
}
}
}

#endif