#ifndef _APP_CALCUL_EDGECLEANINGOP_H_
#define _APP_CALCUL_EDGECLEANINGOP_H_

//SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>
#include <ign/tools/stringtools.h>

//EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>

//APP
#include <app/calcul/detail/EdgeCleaningGraphManager.h>
#include <app/tools/StringTools.h>

namespace app{
namespace calcul{

	class EdgeCleaningOp {

	public:

		typedef app::calcul::detail::EdgeCleaningGraphManager::GraphType   GraphType;
		typedef typename GraphType::edge_descriptor                        edge_descriptor;
		typedef typename GraphType::oriented_edge_descriptor               oriented_edge_descriptor;
		typedef typename GraphType::edge_iterator                          edge_iterator;
		typedef typename GraphType::face_descriptor                        face_descriptor;
		typedef typename GraphType::face_iterator                          face_iterator;
		typedef typename GraphType::vertex_descriptor                      vertex_descriptor;
		typedef typename GraphType::vertex_iterator                        vertex_iterator;
		typedef typename GraphType::edges_path                             edges_path;
		typedef typename GraphType::edges_path_const_iterator              edges_path_const_iterator;
		typedef app::calcul::detail::OriginEdgeProperties                  OriginEdgeProperties;

		/// \brief
		EdgeCleaningOp( 
			std::string borderCode,
            bool verbose 
        );

		/// \brief
		~EdgeCleaningOp();

		/// \brief
		void cleanFaces() const;

		/// \brief
		void cleanFaces2ByCountry() const;

		/// \brief
		void cleanFacesAndAntennaByCountry(std::set<std::string> & sTreatedFeatures, ign::feature::FeatureFilter filter = ign::feature::FeatureFilter()) const;

		/// \brief
		bool cleanFaces2(ign::feature::FeatureFilter filter = ign::feature::FeatureFilter()) const;

		/// \brief
		void cleanPathsOutOfCountry() const;

        /// \brief
		bool cleanAntennas(std::set<std::string> & sTreatedFeatures) const;

		/// \brief
		void cleanParalelleEdges() const;

	private:
		//--
		ign::feature::sql::FeatureStorePostgis*              _fsCp;
		//--
		ign::feature::sql::FeatureStorePostgis*              _fsEdge;
		//--
		std::map<std::string, ign::geometry::GeometryPtr>    _mCountryGeomPtr;
		//--
		std::map<std::string, ign::geometry::GeometryPtr>    _mCountryGeomWithBuffPtr;
		//--
		epg::log::EpgLogger*                                 _logger;
		//--
		epg::log::ShapeLogger*                               _shapeLogger;
		//--
		std::string                                          _countryCode;
		//--
		bool                                                 _verbose;

	private:

		//--
		void _init();

		//--
		void _loadGraph(
			app::calcul::detail::EdgeCleaningGraphManager & graphManager,
			bool planarize,
			bool simplifiedPlanarization = false,
			ign::feature::FeatureFilter filter = ign::feature::FeatureFilter()
		) const;

		//--
		std::pair<double, double> _getLengths( ign::geometry::Geometry const& geom, ign::geometry::Point const* startPoint = 0 ) const;

		//--
		double _getAntennaLength(GraphType const& graph, std::list<oriented_edge_descriptor> const& lEdges) const;

		//--
		std::pair<double, double> _getRatioAndLengthFirstPart(GraphType const& graph, std::string country, std::list<oriented_edge_descriptor> const& path) const;

		//--
		template < typename ContainerType >
		double _getRatio(
			GraphType const& graph, 
			std::string country,
			ContainerType const& container
		) const {
			bool isCl = country.find("#") != std::string::npos;
			if (isCl) return 1;
			
            double lengthInCountry = 0;
            double length = 0;
            typename ContainerType::const_iterator lit = container.begin();
            for ( ; lit != container.end() ; ++lit) {
                ign::geometry::LineString edgeGeom = graph.getGeometry(*lit);
                _addLengths(country, edgeGeom, lengthInCountry, length);
            }
			if (length == 0)
				return 0;

            return lengthInCountry / length;
        }

		//--
		double _getLengthWithBuff( ign::geometry::Geometry const& geom ) const;

		//--
		double _getRatioWithBuff(GraphType const& graph, std::string country, std::list<oriented_edge_descriptor> const& lEdges) const;

		//--
        void _removePath(
			GraphType & graph, 
			std::list<oriented_edge_descriptor> const& path, 
			std::set<edge_descriptor>& sEdge2Remove
		) const;

		//--
		template < typename ContainerType >
        void _removeEdges(GraphType & graph, ContainerType const& container, std::set<edge_descriptor>& sEdge2Remove) const
        {
			std::set<std::string> sFeature2Delete;

			std::string currentOrigin = "";

            std::list<edge_descriptor>::const_iterator lit = container.begin();
            for ( ; lit != container.end() ; ++lit) {
                if (graph.origins(*lit).size() > 1) {
                    _logger->log(epg::log::WARN, "Edge with multiple origins [edge id] "+tools::StringTools::toString(graph.origins(*lit)));
                }

				std::set<std::string> sOrigins;
				for (size_t i = 0 ; i < graph.origins(*lit).size() ; ++i)
					sOrigins.insert(graph.origins(*lit)[i]);

				if (sOrigins.find(currentOrigin) != sOrigins.end()) continue;
				
				currentOrigin = *sOrigins.begin();

				if (ign::tools::StringManip::FindSubString(currentOrigin,"CONNECTINGLINE")) {
					_logger->log(epg::log::WARN, "Edge has a cl as origin [cl id] "+currentOrigin);
					continue;
				}
				
				sFeature2Delete.insert(currentOrigin);
                
				if (graph.origins(*lit).size() == 1 )
					sEdge2Remove.insert(*lit);
            }

			for( std::set<std::string>::const_iterator sit = sFeature2Delete.begin(); sit != sFeature2Delete.end() ; ++sit) {
				ign::feature::Feature dFeat;
				_fsEdge->getFeatureById(*sit, dFeat);
				if (!dFeat.getId().empty())
					_shapeLogger->writeFeature("ecl_deleted_edges", dFeat);

				_fsEdge->deleteFeature(*sit);
			}
			
        }

		//--
		void _removePathAndGraphEdges(
            GraphType & graph,
            std::list<oriented_edge_descriptor> const& path
        ) const;

		//--
		template < typename ContainerType >
        void _removeEdgesAndGraphEdges(GraphType & graph, ContainerType const& container) const
        {
			std::set<edge_descriptor> sEdge2Remove;
			_removeEdges(graph, container, sEdge2Remove);

			
			for ( std::set<edge_descriptor>::const_iterator sit = sEdge2Remove.begin() ; sit != sEdge2Remove.end() ; ++sit )
				graph.removeEdge(*sit);
		}

		//--
		std::pair<double, double> _addLengths(
			std::string country,
			ign::geometry::LineString const& ls,
			double & lengthInCountry,
			double & length
		) const;

		//--
		void _addLengthsWithBuff(
			std::string country,
			ign::geometry::LineString const& ls,
			double & lengthInCountry,
			double & length
		) const;

		//--
		bool _vertexIsConnected2Cl(detail::EdgeCleaningGraphManager const& graphManager, vertex_descriptor v) const;

		//--
		bool _vertexIsCp(GraphType const& graph, vertex_descriptor v) const;

		//--
		bool _cleanAntenna(
            GraphType & graph,
            std::string const& country,
            std::list<oriented_edge_descriptor> const& lAntennas,
            bool bAntennaIsConnected2CF
        ) const;
		
		//--
		bool _isSlimSurface( 
            ign::geometry::Polygon const& poly, 
            double maxWidth,
			ign::geometry::Point const ** p1 = 0,
            ign::geometry::Point const ** p2 = 0
		) const;

		//--
		bool _isSlimSurface( 
            ign::geometry::LineString const& closedLs, 
            double maxWidth
		) const;

		//--
        bool _getFacePaths(
            detail::EdgeCleaningGraphManager const& graphManager, 
            face_descriptor fd, 
            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> & vpCountryEdges
        ) const;

		//--
		double _getPathLength(
            GraphType const& graph, 
            std::list<oriented_edge_descriptor> const& path
        ) const;

		//--
		std::set<std::string> _mergeFacePaths(std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> & vpCountryEdges) const;

		//--
		bool _cleanFaces2(detail::EdgeCleaningGraphManager & graphManager) const;

		//--
		bool _cleanGraphFaces(detail::EdgeCleaningGraphManager & graphManager) const;

		//--
		void _cleanFacesAntennas(
            detail::EdgeCleaningGraphManager & graphManager, 
            std::set<vertex_descriptor> const& sVertices
        ) const ;

        //--
        std::pair<bool, std::list<app::calcul::detail::EdgeCleaningGraphManager::oriented_edge_descriptor>> _getAntenna(
            detail::EdgeCleaningGraphManager const& graphManager,
            vertex_descriptor v,
            std::set<std::string> & sTreatedFeatures,
            bool isPlanarGraph
        ) const;

		//--
		bool _cleanAntennas(
			detail::EdgeCleaningGraphManager & graphManager,
			// std::set<vertex_descriptor> & sTreatedDangles,
			std::set<std::string> & sTreatedFeatures,
			bool isPlanarGraph = false
		) const;

		//--
		bool _cleanGraphAntennas(
            detail::EdgeCleaningGraphManager & graphManager,
			// std::set<vertex_descriptor> & sTreatedDangles,
			std::set<std::string> & sTreatedFeatures,
			bool isPlanarGraph = false
        ) const;

		//--
		typename app::calcul::detail::EdgeCleaningGraphManager::GraphType::oriented_edge_descriptor _getNextEdge(GraphType const& graph, edge_descriptor e, vertex_descriptor vTarget, bool isPlanarGraph) const;

		//--
		template < typename ContainerType >
        ContainerType _getReversePath(ContainerType const& path) const {
            ContainerType reversedPath;

            for (typename ContainerType::const_reverse_iterator vit = path.rbegin() ; vit != path.rend() ; ++vit ) {
                ign::graph::OrientedEdge<edge_descriptor> newOe(vit->descriptor, vit->direction == ign::graph::DIRECT ? ign::graph::REVERSE : ign::graph::DIRECT);
                reversedPath.push_back(newOe);
            }

            return reversedPath;
        }

        //--
        void _addAntennaEdges(GraphType const& graph, oriented_edge_descriptor oe, std::list<oriented_edge_descriptor> & lEdges, bool isPlanarGraph) const;

        //--
        typename app::calcul::detail::EdgeCleaningGraphManager::GraphType::vertex_descriptor _getTarget(GraphType const& graph, oriented_edge_descriptor oe, bool isPlanarGraph) const;
    };

}
}

#endif