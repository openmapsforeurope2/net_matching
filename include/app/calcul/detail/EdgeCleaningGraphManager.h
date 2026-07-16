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

    /// @brief Structure permettant de conserver en mémoire les propriétés des données d'origines
    struct OriginEdgeProperties {

        /// \brief
        OriginEdgeProperties(std::string country_, bool isCl_ = false):isCl(isCl_), country(country_), wTag(""), natId(""){};

        /// \brief
        OriginEdgeProperties(std::string country_, std::string wTag_, bool isCl_ = false):isCl(isCl_), country(country_), wTag(wTag_), natId(""){};

        /// \brief
        OriginEdgeProperties(std::string country_, std::string wTag_, std::string natId_, bool isCl_ = false):isCl(isCl_), country(country_), wTag(wTag_), natId(natId_){};

        /// \brief
        ~OriginEdgeProperties(){};

        //--
        bool		  isCl;
        std::string   country;
        std::string   wTag;
        std::string   natId;
    };

    /// @brief Classe utilitaire facilitant la manipulation du graph de travail
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

        /// @brief Constructeur
        EdgeCleaningGraphManager() :
            _builder(0),
            _simpleBuilder(0),
            _simplifiedPlanarization(false)
        {
        };

        /// @brief Destructeur
        ~EdgeCleaningGraphManager() {
            if (_builder) delete _builder;
            if (_simpleBuilder) delete _simpleBuilder;
        };

        /// @brief Setter pour indiquer si la planarisation doit être simplifiée
        /// @param simplified Booléen
        void setSimplifiedPlanarization(bool simplified) {
            _simplifiedPlanarization = simplified;
        }

        /// @brief Retourne le graph encapsulé
        /// @return Graph encapsulé
        inline GraphType const& getGraph() const {
            return _graph;
        }

        /// @brief Retourne le graph encapsulé
        /// @return Graph encapsulé
        inline GraphType & getGraph() {
            return _graph;
        }

        /// @brief Ajoute un arc au graph
        /// @param ls Géométrie de l'arc
        /// @param idOrigin Identifiant d'origine
        /// @param edgeProperties Propriété de l'arc d'origine
        /// @return Booléen indiquant si l'arc a bien été ajouté
        bool addEdge( 
            ign::geometry::LineString const& ls,
            std::string const& idOrigin,
            OriginEdgeProperties const& edgeProperties
        ) {
            if ( !_builder ) _builder = new ign::geometry::graph::tools::SnapRoundPlanarizer<GraphType>(_graph);
            _mEdges.insert(std::make_pair(idOrigin, edgeProperties));
            return _builder->addEdge(ls, idOrigin);
        };

        /// @brief supprime un arc du graph
        /// @param e arc à supprimer
        void removeEdge( 
            edge_descriptor e
        ) {
            std::vector< std::string > const& vOrigins = _graph.origins(e);

            _graph.removeEdge(e, true);

            for (std::vector< std::string >::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit){
                if( _graph.isLinearOrigin(*vit) ) continue; //not removed
                
                _mEdges.erase(*vit);
            }
        };

        /// @brief Ajoute un arc au graph
        /// @param vSource Identifiant du sommet source
        /// @param vTarget Identifiant du sommet cible
        /// @param intermediatePoints Points intermédiaires constituants la géométrie de l'arc
        /// @param idOrigin Identifiant d'origine
        /// @param edgeProperties Propriété de l'arc d'origine
        /// @return Identifiant de l'arc ajouté
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

        /// @brief Ajoute un arc au graph en utiisant le SimpleGraphBuilder
        /// @param ls Géométrie de l'arc
        /// @param idOrigin Identifiant d'origine
        /// @param edgeProperties Propriété de l'arc d'origine
        /// @return Identifiant de l'arc ajouté
        oriented_edge_descriptor addEdgeSimple( 
            ign::geometry::LineString const& ls,
            std::string const& idOrigin,
            OriginEdgeProperties const& edgeProperties
        ) {
            if ( !_simpleBuilder ) _simpleBuilder = new ign::geometry::graph::builder::SimpleGraphBuilder<GraphType>(_graph, 1e-5);
            _mEdges.insert(std::make_pair(idOrigin, edgeProperties));
            return _simpleBuilder->addEdge(ls, idOrigin);
        };

        /// @brief Lance la planarisation du graph
        void planarize() {
            if (_builder) _builder->planarize();
        };

        /// @brief Lance l'affectation du poids des arcs (longueur)
        void initWeight() {
            edge_iterator eit, eend;
            for (_graph.edges(eit, eend) ; eit != eend ; ++eit) {
                _graph[*eit].weight = _graph.getGeometry(*eit).length();
            }
        };

        /// @brief Lance la création des faces sur le graph planaire
        void createFaces() {
            _graph.createFaces();
        };

        /// @brief Nettoyage de la classe (graph et builders)
        void clear() {
            _graph.clear();
            delete _builder;
            delete _simpleBuilder;
            _builder = 0;
            _simpleBuilder = 0;
        };

        /// @brief Indique si le sommet est un connecting point
        /// @param v Identifiant du sommet
        /// @return Booléen
        bool isCp(vertex_descriptor v) const {
            return _graph[v].isCp;
        };

        /// @brief Indique si l'arc est une connecting line
        /// @param e Identifiant de l'arc
        /// @return Booléen
        bool isCl(edge_descriptor e) const {
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            for (std::vector< std::string >::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit){
                std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(*vit);
                if ( mit != _mEdges.end() && mit->second.isCl ) return true;
            }
            return false;
        };

        /// @brief Retourne les linéaires d'origine qui sont des CL
        /// @param e Identifiant de l'arc
        /// @return Booléen
        std::vector<std::string> getClOrigins(edge_descriptor e) const {
            std::vector<std::string> vCl;
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            for (std::vector< std::string >::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit){
                std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(*vit);
                if ( mit != _mEdges.end() && mit->second.isCl ) vCl.push_back(*vit);
            }
            return vCl;
        };

        /// @brief Indique si le sommet touche une connecting line
        /// @param v Indentifiant du sommet
        /// @return Booléen
        bool isTouchingCl(vertex_descriptor v) const {
            std::vector< edge_descriptor > vIncidents = _graph.incidentEdges(v);
            for (std::vector<edge_descriptor>::const_iterator vit = vIncidents.begin() ; vit != vIncidents.end() ; ++vit) {
                if(this->isCl(*vit)) {
                    return true;
                }
            }
            return false;
        };

        /// @brief Retourne les linéaires incidents
        /// @param v Indentifiant du sommet
        /// @return Booléen
        // std::vector<std::string> getIncidentOrigins(vertex_descriptor v) const {
        //     std::vector<std::string> vIncidentOrigins;
        //     std::vector< edge_descriptor > vIncidents = _graph.incidentEdges(v);
        //     for (std::vector<edge_descriptor>::const_iterator vit = vIncidents.begin() ; vit != vIncidents.end() ; ++vit) {
        //         std::vector< std::string > const& vOrigins = _graph.origins(*vit);
        //         for(std::vector< std::string >::const_iterator vit2 = vOrigins.begin() ; vit2 != vOrigins.end() ; ++vit2) {
        //             std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdges = _graph.getInducedEdges(*vit2);
        //             if(!foundInducedEdges.first)
        //                 continue;
                    
        //             if(_graph.source(foundInducedEdges.second.begin) == v || _graph.target(foundInducedEdges.second.rbegin) == v)
        //                 vIncidentOrigins.push_back(*vit2);
        //             }
        //         }
        //     }
        //     return vIncidentOrigins;
        // };

        /// @brief Retourne la valeur du champ de travail W_TAG_NAME de l'arc
        /// @param e Identifiant de l'arc
        /// @return Valeur du champ W_TAG_NAME
        std::string getWTag(edge_descriptor e) const {
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(vOrigins.front());
            if ( mit != _mEdges.end() ) return mit->second.wTag;
            return "";
        };

        /// @brief Retourne la valeur du champ de travail W_NATIONAL_IDENTIFIER_NAME de l'arc
        /// @param e Identifiant de l'arc
        /// @return Valeur du champ W_NATIONAL_IDENTIFIER_NAME
        std::string getNatId(edge_descriptor e) const {
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(vOrigins.front());
            if ( mit != _mEdges.end() ) return mit->second.natId;
            return "";
        };

        /// @brief Retourne la valeur du champ de travail W_NATIONAL_IDENTIFIER_NAME du linéaire
        /// @param id Identifiant du linéaire
        /// @return Valeur du champ W_NATIONAL_IDENTIFIER_NAME
        std::string getNatId(std::string const& id) const {
            std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(id);
            if ( mit != _mEdges.end() ) return mit->second.natId;
            return "";
        };

        /// @brief Retourne un code pays de l'arc
        /// @param e Identifiant de l'arc
        /// @return Code pays
        std::string getCountry(edge_descriptor e) const {
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(vOrigins.front());
            if ( mit != _mEdges.end() ) return mit->second.country;
            return "";
        };

        /// @brief Retourne un code pays du linéaire
        /// @param id Identifiant du linéaire
        /// @return Code pays
        std::string getCountry(std::string const& id) const {
            std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(id);
            if ( mit != _mEdges.end() ) return mit->second.country;
            return "";
        };

        /// @brief Retourne les codes pays de l'arc (plusieurs codes sont possibles si l'arc du graph possède plusieurs linéaires d'origine)
        /// @param e Identifiant de l'ars
        /// @return Codes pays
        std::set<std::string> getCountries(edge_descriptor e) const {
            std::set<std::string> sCountry;
            std::vector< std::string > const& vOrigins = _graph.origins(e);
            for(std::vector< std::string >::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit)
            {
                std::map<std::string, OriginEdgeProperties>::const_iterator mit = _mEdges.find(*vit);
                if ( mit != _mEdges.end() ) sCountry.insert(mit->second.country) ;
            }
            return sCountry;
        };

        /// @brief Retourne les codes pays simples de l'arc (en supprimant les '#')
        /// @param e Identifiant de l'arc
        /// @return Codes pays simples
        std::set<std::string> getSingleCountries(edge_descriptor e) const {
            std::set<std::string> sCountry = getCountries(e);
            std::set<std::string> sSingleCountry;

            for(std::set<std::string>::const_iterator sit = sCountry.begin() ; sit != sCountry.end() ; ++sit)
            {
                std::vector<std::string> vCountry;
                epg::tools::StringTools::Split(*sit, "#", vCountry);

                sSingleCountry.insert(vCountry.begin(), vCountry.end());
            }
            return sSingleCountry;
        };

        /// @brief Retourne la geometrie de l'arc d'origine reconstruite à partir des edges induits.
        /// Si des edges induits ont été supprimés et que la reconstruction de la géométrie résulte
        /// en plusieurs polylignes, c'est la polyligne la plus longue qui est retournée.
        /// @param id identifiant de la géométrie d'origine
        /// @return une paire avec un booléen indiquant si la géomtrie à pu être reconstruite et la géométrie reconstruite
        std::pair<bool, ign::geometry::LineString> getLongestInducedLineString(std::string const& id) const {
            std::vector<ign::geometry::LineString> vInducedLs(1, ign::geometry::LineString());
            
            std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdges = _graph.getInducedEdges(id);
            if(!foundInducedEdges.first)
                return std::make_pair(false, ign::geometry::LineString());
            
            std::vector< oriented_edge_descriptor > const& vEdges = foundInducedEdges.second;
            
            vertex_descriptor vdTarget = 0;
            for ( size_t i = 0 ; i < vEdges.size() ; ++i )
            {
                ign::graph::EdgeDirection direction = vEdges[ i ].direction;
                edge_descriptor ed = vEdges[ i ].descriptor;
                
                vertex_descriptor vdSource = ( direction == ign::graph::DIRECT )? _graph.source( ed ) : _graph.target( ed );
                
                if ( i != 0 && vdSource != vdTarget ) {
                    vInducedLs.back().addPoint(_graph[vdTarget].point);
                    vInducedLs.push_back(ign::geometry::LineString());
                }
                
                vdTarget = ( direction == ign::graph::DIRECT )? _graph.target( ed ) : _graph.source( ed );
                
                vInducedLs.back().addPoint( _graph[vdSource].point );
                
                std::vector< ign::geometry::Point > const& vPoints = _graph[ ed ].intermediatePoints;
                
                if( direction == ign::graph::DIRECT )
                {
                    std::vector< ign::geometry::Point >::const_iterator it = vPoints.begin();
                    for( ; it != vPoints.end() ; ++it )
                         vInducedLs.back().addPoint( *it );
                }
                else{
                    std::vector< ign::geometry::Point >::const_reverse_iterator it = vPoints.rbegin();
                    for( ; it != vPoints.rend() ; ++it )
                         vInducedLs.back().addPoint( *it );
                }
            }
            vInducedLs.back().addPoint( _graph[vdTarget].point );

            size_t idLongest = 0;
            if (vInducedLs.size() > 1) {
                double longestLength = 0;
                for ( size_t i = 0 ; i < vInducedLs.size() ; ++i ) {
                    double length = vInducedLs[i].length();
                    if( length > longestLength) {
                        longestLength = length;
                        idLongest = i;
                    }
                }
            }
            
            return std::make_pair(true, vInducedLs[idLongest]);
        }
        
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