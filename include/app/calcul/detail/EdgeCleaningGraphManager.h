#ifndef _APP_CALCUL_DETAIL_EDGECLEANINGGRAPHMANAGER_H_
#define _APP_CALCUL_DETAIL_EDGECLEANINGGRAPHMANAGER_H_

// SOCLE
#include <ign/geometry/graph/GeometryGraph.h>
#include <ign/geometry/graph/tools/SnapRoundPlanarizer.h>

namespace app{
namespace calcul{
namespace detail{

    struct EdgeCleaningEdge {
        /// \brief
        EdgeCleaningEdge(std::string country_, bool isCl_ = false):isCl(isCl_), country(country_){};

        /// \brief
        ~EdgeCleaningEdge(){};

        //--
        bool		  isCl;
        std::string   country;
    };

    class EdgeCleaningGraphManager {

        public:
            typedef ign::geometry::graph::PunctualVertexProperties  vertexProperties;
            typedef ign::geometry::graph::LinearEdgeProperties      edgeProperties;

            typedef ign::geometry::graph::GeometryGraph<vertexProperties, edgeProperties>   GraphType;
            typedef typename GraphType::edge_descriptor             edge_descriptor;
            typedef typename GraphType::oriented_edge_descriptor    oriented_edge_descriptor;
            typedef typename GraphType::face_descriptor             face_descriptor;
            typedef typename GraphType::face_iterator               face_iterator;
            typedef typename GraphType::vertex_descriptor           vertex_descriptor;
		    typedef typename GraphType::vertex_iterator             vertex_iterator;
        
        private:

            std::map<std::string, EdgeCleaningEdge>                        _mEdges;
            GraphType                                                      _graph;
            ign::geometry::graph::tools::SnapRoundPlanarizer<GraphType>*   _builder;

        public:


            //--
            EdgeCleaningGraphManager() {
                _builder = new ign::geometry::graph::tools::SnapRoundPlanarizer<GraphType>(_graph);
            };

            //--
            ~EdgeCleaningGraphManager() {
                delete _builder;
            };

            inline GraphType const& getGraph() const {
                return _graph;
            }

            //--
            bool addEdge( 
                ign::geometry::LineString const& ls,
                std::string const& idOrigin,
                EdgeCleaningEdge const& edgeProperties
			) {
                _mEdges.insert(std::make_pair(idOrigin, edgeProperties));
                return _builder->addEdge(ls, idOrigin);
            };

            //--
            void planarize() {
                _builder->planarize();
            };

            //--
            void createFaces() {
                _graph.createFaces();
            };

            //--
            bool isCl(edge_descriptor e) {
                std::vector< std::string > const& vOrigins = _graph.origins(e);
                for (std::vector< std::string >::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit){
                    std::map<std::string, EdgeCleaningEdge>::const_iterator mit = _mEdges.find(*vit);
                    if ( mit != _mEdges.end() && mit->second.isCl ) return true;
                }
                return false;
            };

            //--
            std::string getCountry(edge_descriptor e) {
                std::vector< std::string > const& vOrigins = _graph.origins(e);
                std::map<std::string, EdgeCleaningEdge>::const_iterator mit = _mEdges.find(vOrigins.front());
                if ( mit != _mEdges.end() ) return mit->second.country;
                return "";
            };
    };

}
}
}

#endif