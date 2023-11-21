#ifndef _APP_CALCUL_EDGECONNECTOROP_H_
#define _APP_CALCUL_EDGECONNECTOROP_H_

// SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>

// EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>

// APP
#include <app/calcul/detail/EdgeCleaningGraphManager.h>
#include <app/geometry/tools/LengthIndexedLineString.h>


namespace app{
namespace calcul{

	class EdgeConnectorOp {

	public:

		typedef app::calcul::detail::EdgeCleaningGraphManager::GraphType   GraphType;
		typedef typename GraphType::edge_descriptor                        edge_descriptor;
		typedef typename GraphType::oriented_edge_descriptor               oriented_edge_descriptor;
		// typedef typename GraphType::edge_iterator                          edge_iterator;
		// typedef typename GraphType::face_descriptor                        face_descriptor;
		// typedef typename GraphType::face_iterator                          face_iterator;
		typedef typename GraphType::vertex_descriptor                      vertex_descriptor;
		typedef typename GraphType::vertex_iterator                        vertex_iterator;
		// typedef typename GraphType::edges_path                             edges_path;
		// typedef typename GraphType::edges_path_const_iterator              edges_path_const_iterator;
		typedef typename GraphType::linear_origin_iterator                 linear_origin_iterator;
		typedef app::calcul::detail::OriginEdgeProperties                  OriginEdgeProperties;

		/// \brief
		static void compute(
			std::string borderCode, 
			bool verbose
		);

		
		

	private:
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsLandmask;
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsEdge;
		//--
		// ign::feature::sql::FeatureStorePostgis*            _fsCp;
		// //--
		// ign::feature::sql::FeatureStorePostgis*            _fsCl;
		//--
		std::map<std::string, ign::geometry::GeometryPtr>    _mCountryGeomPtr;
		//--
		std::map<std::string, ign::geometry::GeometryPtr>    _mCountryGeomWithBuffPtr;
		//--
		epg::log::EpgLogger*                               _logger;
		//--
		epg::log::ShapeLogger*                             _shapeLogger;
		//--
		std::string                                        _borderCode;
		//--
		bool                                               _verbose;

	private:

		//--
		EdgeConnectorOp( 
            std::string borderCode, 
            bool verbose 
        );

		//--
		~EdgeConnectorOp();

		//--
		void _init();

		//--
		void _compute();

		//--
		void _loadGraph(
			app::calcul::detail::EdgeCleaningGraphManager & graphManager,
			bool planarize,
			ign::feature::FeatureFilter filter = ign::feature::FeatureFilter()
		) const;

		//--
		void _loadGraphAndPlanarize(
            app::calcul::detail::EdgeCleaningGraphManager & graphManager,
            std::string const& countryCode
        ) const;

		//--
		double _getZ( geometry::tools::LengthIndexedLineString const& lengthIndexedLs, ign::geometry::Point const& p) const;

    };

}
}

#endif