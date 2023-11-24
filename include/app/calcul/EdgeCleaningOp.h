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
			std::string edgeTable,
			std::string borderCode,
			bool verbose
		);

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
		EdgeCleaningOp( 
            std::string edgeTable,
			std::string borderCode,
            bool verbose 
        );

		//--
		~EdgeCleaningOp();

		//--
		void _init( 
            std::string edgeTable,
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
		double _getRatio(GraphType const& graph, std::string country, std::list<edge_descriptor> const& lEdges) const;

		//--
		double _getLengthWithBuff( ign::geometry::Geometry const& geom ) const;

		//--
		double _getRatioWithBuff(GraphType const& graph, std::string country, std::list<edge_descriptor> const& lEdges) const;

		//--
		template < typename ContainerType >
        void _removeEdges(GraphType const& graph, ContainerType const& container) const
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
                    _shapeLogger->writeFeature("edge_cleaning_deleted_edges", dFeat);

                    _fsEdge->deleteFeature(edgeId);
                }
            }
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
		void _cleanFaces() const;

		//--
		void _cleanPathsOutOfCountry() const;

        //--
		void _cleanAntennas() const;

        //--
		void _cleanParalelleEdges() const;

		//--
		void _clean() const;
    };

}
}

#endif