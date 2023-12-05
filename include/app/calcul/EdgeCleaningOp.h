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
		static void clean(
			std::string borderCode,
			bool verbose
		);

		/// \brief
		EdgeCleaningOp( 
			std::string borderCode,
            bool verbose 
        );

		/// \brief
		~EdgeCleaningOp();

		/// \brief
		void cleanAll() const;

		/// \brief
		void cleanFaces() const;

		/// \brief
		void cleanFaces2() const;

		/// \brief
		void cleanPathsOutOfCountry() const;

        /// \brief
		void cleanAntennas() const;

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
		void _init( 
			std::string borderCode
        );

		//--
		void _loadGraph(
			app::calcul::detail::EdgeCleaningGraphManager & graphManager,
			bool planarize,
			ign::feature::FeatureFilter filter = ign::feature::FeatureFilter()
		) const;

		//--
		double _getLength( ign::geometry::Geometry const& geom ) const;

		//--
		double _getAntennaLength(GraphType const& graph, std::list<edge_descriptor> const& lEdges) const;

		//--
		double _getRatio(GraphType const& graph, std::string country, std::list<oriented_edge_descriptor> const& path) const;

		//--
		double _getRatio(GraphType const& graph, std::string country, std::list<edge_descriptor> const& lEdges) const;

		//--
		double _getLengthWithBuff( ign::geometry::Geometry const& geom ) const;

		//--
		double _getRatioWithBuff(GraphType const& graph, std::string country, std::list<edge_descriptor> const& lEdges) const;

		//--
        void _removePath(
			GraphType & graph, std::list<oriented_edge_descriptor> const& path, 
			std::list<edge_descriptor>& lEdge2Remove
		) const;

		//--
		template < typename ContainerType >
        void _removeEdges(GraphType & graph, ContainerType const& container, std::list<edge_descriptor>& lEdge2Remove) const
        {
            std::list<edge_descriptor>::const_iterator lit = container.begin();
            for ( ; lit != container.end() ; ++lit) {
                if (graph.origins(*lit).size() > 1) {
                    _logger->log(epg::log::WARN, "Edge with multiple origins [edge id] "+tools::StringTools::toString(graph.origins(*lit)));
                }
                for (size_t i = 0 ; i < graph.origins(*lit).size() ; ++i) {
                    std::string edgeId = graph.origins(*lit)[i];

                    if (ign::tools::StringManip::FindSubString(edgeId,"CONNECTINGLINE")) {
                        _logger->log(epg::log::WARN, "Edge has a cl as origin [cl id] "+edgeId);
                        continue;
                    }

                    ign::feature::Feature dFeat;
                    _fsEdge->getFeatureById(edgeId, dFeat);
                    _shapeLogger->writeFeature("ecl_deleted_edges", dFeat);

                    _fsEdge->deleteFeature(edgeId);
                }

				lEdge2Remove.push_back(*lit);
            }
        }

		//--
		template < typename ContainerType >
        void _removeEdgesAndGraphEdges(GraphType & graph, ContainerType const& container) const
        {
			std::list<edge_descriptor> lEdge2Remove;
			_removeEdges(graph, container, lEdge2Remove);

			for ( std::list<edge_descriptor>::const_iterator lit = lEdge2Remove.begin() ; lit != lEdge2Remove.end() ; ++lit )
				graph.removeEdge(*lit);
		}

		//--
		void _addLengths(
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
		void _cleanAntenna(
            GraphType & graph,
            std::string const& country,
            std::list<edge_descriptor> const& lAntennas,
            bool bAntennaIsConnected2CF
        ) const;

		//--
		std::pair<ign::geometry::LineString, ign::geometry::LineString> _getSubLineStrings(
			size_t id1, 
            size_t id2, 
			ign::geometry::LineString const& ls
		) const;

		//--
		bool _isSlimSurface( 
            ign::geometry::Polygon const& poly, 
            double maxWidth
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
    };

}
}

#endif