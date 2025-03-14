#ifndef _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGEDGE_H_
#define _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGEDGE_H_


// SOCLE
#include <ign/geometry/Point.h>
#include <ign/geometry/graph/LinearEdgeProperties.h>



namespace app{
namespace calcul{
namespace detail{
namespace graph{

	/// @brief 
	struct EdgeCleaningEdge : public ign::geometry::graph::LinearEdgeProperties {

        /// \brief constructor
        EdgeCleaningEdge(){};

        /// \brief destructor
        ~EdgeCleaningEdge(){};

        //--
        double        weight;
    };
}
}
}
}

#endif