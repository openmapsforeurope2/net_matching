#ifndef _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGEDGE_H_
#define _APP_CALCUL_DETAIL_GRAPH_EDGECLEANINGEDGE_H_


//SOCLE
#include <ign/geometry/Point.h>
#include <ign/geometry/graph/LinearEdgeProperties.h>



namespace app{
namespace calcul{
namespace detail{
namespace graph{

	struct EdgeCleaningEdge : public ign::geometry::graph::LinearEdgeProperties {
        /// \brief
        EdgeCleaningEdge(){};

        /// \brief
        ~EdgeCleaningEdge(){};

        //--
        double        weight;
    };
}
}
}
}

#endif