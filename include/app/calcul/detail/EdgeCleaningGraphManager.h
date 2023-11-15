#ifndef _APP_CALCUL_DETAIL_EDGECLEANINGGRAPHMANAGER_H_
#define _APP_CALCUL_DETAIL_EDGECLEANINGGRAPHMANAGER_H_

// SOCLE
#include <ign/geometry/graph/tools/SnapRoundPlanarizer.h>
#include <ign/geometry/graph/builder/SimpleGraphBuilder.h>

// APP
#include <app/calcul/detail/graph/EdgeCleaningGraph.h>

namespace app{
namespace calcul{
namespace detail{

    	struct OriginEdgeProperties {
        /// \brief
        OriginEdgeProperties(std::string country_, bool isCl_ = false):isCl(isCl_), country(country_){};

        /// \brief
        ~OriginEdgeProperties(){};

        //--
        bool		  isCl;
        std::string   country;
    };

    class EdgeCleaningGraphManager {

        public:
            typedef typename graph::EdgeCleaningGraph               GraphType;
            typedef typename GraphType::edge_descriptor             edge_descriptor;
            typedef typename GraphType::edge_iterator               edge_iterator;
            typedef typename GraphType::oriented_edge_descriptor    oriented_edge_descriptor;
            typedef typename GraphType::face_descriptor             face_descriptor;
            typedef typename GraphType::face_iterator               face_iterator;
            typedef typename GraphType::vertex_descriptor           vertex_descriptor;
		    typedef typename GraphType::vertex_iterator             vertex_iterator;
            typedef typename GraphType::edges_path                  edges_path;
		    typedef typename GraphType::edges_path_const_iterator   edges_path_const_iterator;
        
        private:

            std::map<std::string, OriginEdgeProperties>                    _mEdges;
            GraphType                                                      _graph;
            ign::geometry::graph::tools::SnapRoundPlanarizer<GraphType>*   _builder;
            ign::geometry::graph::builder::SimpleGraphBuilder<GraphType>*  _simpleBuilder;

        public:


            //--
            EdgeCleaningGraphManager() {
                _builder = new ign::geometry::graph::tools::SnapRoundPlanarizer<GraphType>(_graph);
                _simpleBuilder = new ign::geometry::graph::builder::SimpleGraphBuilder<GraphType>(_graph, 1e-5);
            };

            //--
            ~EdgeCleaningGraphManager() {
                delete _builder;
                delete _simpleBuilder;
            };

            inline GraphType const& getGraph() const {
                return _graph;
            }

            inline GraphType & getGraph() {
                return _graph;
            }

            //--
            bool addEdge( 
                ign::geometry::LineString const& ls,
                std::string const& idOrigin,
                OriginEdgeProperties const& edgeProperties
			) {
                _mEdges.insert(std::make_pair(idOrigin, edgeProperties));
                return _builder->addEdge(ls, idOrigin);
            };

            //--
            oriented_edge_descriptor addEdge(
                vertex_descriptor vSource,
                vertex_descriptor vTarget,
                std::vector< ign::geometry::Point > const& intermediatePoints,
                std::string const& idOrigin,
                OriginEdgeProperties const& edgeProperties
			) {
                std::map<std::string, OriginEdgeProperties>::iterator mit = _mEdges.find(idOrigin);
                if (mit != _mEdges.end()) mit->second = edgeProperties;
                else _mEdges.insert(std::make_pair(idOrigin, edgeProperties));
                return _graph.addEdge(vSource, vTarget, intermediatePoints);
            };

            //--
            oriented_edge_descriptor addEdgeSimple( 
                ign::geometry::LineString const& ls,
                std::string const& idOrigin,
                OriginEdgeProperties const& edgeProperties
			) {
                _mEdges.insert(std::make_pair(idOrigin, edgeProperties));
                return _simpleBuilder->addEdge(ls, idOrigin);
            };

            //--
            void planarize() {
                _builder->planarize();
            };

            //--
            void initWeight() {
                edge_iterator eit, eend;
                for (_graph.edges(eit, eend) ; eit != eend ; ++eit) {
                    _graph[*eit].weight = _graph.getGeometry(*eit).length();
                }
            };

            //--
            void createFaces() {
                _graph.createFaces();
            };

            //--
            void clear() {
                _graph.clear();
            };

            //--
            bool isCl(edge_descriptor e) const {
                std::vector< std::string > const& vOrigins = _graph.origins(e);
                for (std::vector< std::string >::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit){
                    std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(*vit);
                    if ( mit != _mEdges.end() && mit->second.isCl ) return true;
                }
                return false;
            };

            //--
            bool isTouchingCl(vertex_descriptor v) const {
                std::vector< edge_descriptor > vIncidents = _graph.incidentEdges(v);
                for (std::vector<edge_descriptor>::const_iterator vit = vIncidents.begin() ; vit != vIncidents.end() ; ++vit) {
                    if(this->isCl(*vit)) {
                        return true;
                    }
                }
                return false;
            };

            //--
            std::string getCountry(edge_descriptor e) const {
                std::vector< std::string > const& vOrigins = _graph.origins(e);
                std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(vOrigins.front());
                if ( mit != _mEdges.end() ) return mit->second.country;
                return "";
            };
    };

}
}
}

#endif