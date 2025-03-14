#ifndef _APP_CALCUL_DETAIL_EDGECLEANINGGRAPHMANAGER_H_
#define _APP_CALCUL_DETAIL_EDGECLEANINGGRAPHMANAGER_H_

// SOCLE
#include <ign/geometry/graph/tools/SnapRoundPlanarizer.h>
#include <ign/geometry/graph/builder/SimpleGraphBuilder.h>

// APP
#include <app/calcul/detail/graph/EdgeCleaningGraph.h>

// EPG
#include <epg/tools/StringTools.h>


namespace app{
namespace calcul{
namespace detail{

    /// @brief 
    struct OriginEdgeProperties {

        /// \brief
        OriginEdgeProperties(std::string country_, bool isCl_ = false):isCl(isCl_), country(country_), wTag(""){};

        /// \brief
        OriginEdgeProperties(std::string country_, std::string wTag_, bool isCl_ = false):isCl(isCl_), country(country_), wTag(wTag_){};

        /// \brief
        ~OriginEdgeProperties(){};

        //--
        bool		  isCl;
        std::string   country;
        std::string   wTag;
    };

    /// @brief 
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
        typedef typename GraphType::linear_origin_iterator      linear_origin_iterator;

    public:

        /// @brief 
        EdgeCleaningGraphManager() :
            _builder(0),
            _simpleBuilder(0),
            _simplifiedPlanarization(false)
        {
        };

        /// @brief 
        ~EdgeCleaningGraphManager() {
            if (_builder) delete _builder;
            if (_simpleBuilder) delete _simpleBuilder;
        };

        /// @brief 
        /// @param simplified 
        void setSimplifiedPlanarization(bool simplified) {
            _simplifiedPlanarization = simplified;
        }

        /// @brief 
        /// @return 
        inline GraphType const& getGraph() const {
            return _graph;
        }

        /// @brief 
        /// @return 
        inline GraphType & getGraph() {
            return _graph;
        }

        /// @brief 
        /// @param ls 
        /// @param idOrigin 
        /// @param edgeProperties 
        /// @return 
        bool addEdge( 
            ign::geometry::LineString const& ls,
            std::string const& idOrigin,
            OriginEdgeProperties const& edgeProperties
        ) {
            if ( !_builder ) _builder = new ign::geometry::graph::tools::SnapRoundPlanarizer<GraphType>(_graph);
            _mEdges.insert(std::make_pair(idOrigin, edgeProperties));
            return _builder->addEdge(ls, idOrigin);
        };

        /// @brief 
        /// @param vSource 
        /// @param vTarget 
        /// @param intermediatePoints 
        /// @param idOrigin 
        /// @param edgeProperties 
        /// @return 
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

        /// @brief 
        /// @param ls 
        /// @param idOrigin 
        /// @param edgeProperties 
        /// @return 
        oriented_edge_descriptor addEdgeSimple( 
            ign::geometry::LineString const& ls,
            std::string const& idOrigin,
            OriginEdgeProperties const& edgeProperties
        ) {
            if ( !_simpleBuilder ) _simpleBuilder = new ign::geometry::graph::builder::SimpleGraphBuilder<GraphType>(_graph, 1e-5);
            _mEdges.insert(std::make_pair(idOrigin, edgeProperties));
            return _simpleBuilder->addEdge(ls, idOrigin);
        };

        /// @brief 
        void planarize() {
            if (_builder) _builder->planarize();
        };

        /// @brief 
        void initWeight() {
            edge_iterator eit, eend;
            for (_graph.edges(eit, eend) ; eit != eend ; ++eit) {
                _graph[*eit].weight = _graph.getGeometry(*eit).length();
            }
        };

        /// @brief 
        void createFaces() {
            _graph.createFaces();
        };

        /// @brief 
        void clear() {
            _graph.clear();
            delete _builder;
            delete _simpleBuilder;
            _builder = 0;
            _simpleBuilder = 0;
        };

        /// @brief 
        /// @param e 
        /// @return 
        bool isCl(edge_descriptor e) const {
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            for (std::vector< std::string >::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit){
                std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(*vit);
                if ( mit != _mEdges.end() && mit->second.isCl ) return true;
            }
            return false;
        };

        /// @brief 
        /// @param v 
        /// @return 
        bool isTouchingCl(vertex_descriptor v) const {
            std::vector< edge_descriptor > vIncidents = _graph.incidentEdges(v);
            for (std::vector<edge_descriptor>::const_iterator vit = vIncidents.begin() ; vit != vIncidents.end() ; ++vit) {
                if(this->isCl(*vit)) {
                    return true;
                }
            }
            return false;
        };

        /// @brief 
        /// @param e 
        /// @return 
        std::string getWTag(edge_descriptor e) const {
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(vOrigins.front());
            if ( mit != _mEdges.end() ) return mit->second.wTag;
            return "";
        };

        /// @brief 
        /// @param e 
        /// @return 
        std::string getCountry(edge_descriptor e) const {
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(vOrigins.front());
            if ( mit != _mEdges.end() ) return mit->second.country;
            return "";
        };

        /// @brief 
        /// @param e 
        /// @return 
        std::set<std::string> getCountries(edge_descriptor e) const {
            std::set<std::string> sCountry;
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            for(std::vector< std::string >::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit) {
                std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(*vit);
                if ( mit != _mEdges.end() ) sCountry.insert(mit->second.country) ;
            }
            return sCountry;
        };

        /// @brief 
        /// @param e 
        /// @return 
        std::set<std::string> getSingleCountries(edge_descriptor e) const {
            std::set<std::string> sCountry = getCountries(e);
            std::set<std::string> sSingleCountry;

            for(std::set<std::string>::const_iterator sit = sCountry.begin() ; sit != sCountry.end() ; ++sit) {
                std::vector<std::string> vCountry;
                epg::tools::StringTools::Split(*sit, "#", vCountry);

                sSingleCountry.insert(vCountry.begin(), vCountry.end());
            }
            return sSingleCountry;
        };
        
    private:
        //--
        std::map<std::string, OriginEdgeProperties>                    _mEdges;
        //--
        GraphType                                                      _graph;
        //--
        ign::geometry::graph::tools::SnapRoundPlanarizer<GraphType>*   _builder;
        //--
        ign::geometry::graph::builder::SimpleGraphBuilder<GraphType>*  _simpleBuilder;
        //--
        bool                                                           _simplifiedPlanarization;
    };

    

}
}
}

#endif