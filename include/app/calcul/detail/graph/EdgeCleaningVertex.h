#ifndef _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGVERTEX_H_
#define _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGVERTEX_H_


// SOCLE
#include <ign/geometry/Point.h>
#include <ign/geometry/graph/PunctualVertexProperties.h>



namespace app{
namespace calcul{
namespace detail{
namespace graph{

	/// @brief Propriétés des sommets pour le graph EdgeCleaningGraph
	struct EdgeCleaningVertex : public ign::geometry::graph::PunctualVertexProperties {

        /// \brief constructeur
        EdgeCleaningVertex():isCp(false){};

        /// \brief destructeur
        ~EdgeCleaningVertex(){};

        //--
        bool        isCp;
    };
}
}
}
}

#endif