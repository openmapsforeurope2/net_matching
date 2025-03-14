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

	/// @brief 
	class EdgeConnectorOp {

	public:

		typedef app::calcul::detail::EdgeCleaningGraphManager::GraphType   GraphType;
		typedef typename GraphType::edge_descriptor                        edge_descriptor;
		typedef typename GraphType::oriented_edge_descriptor               oriented_edge_descriptor;
		typedef typename GraphType::vertex_descriptor                      vertex_descriptor;
		typedef typename GraphType::vertex_iterator                        vertex_iterator;
		typedef typename GraphType::linear_origin_iterator                 linear_origin_iterator;
		typedef app::calcul::detail::OriginEdgeProperties                  OriginEdgeProperties;

		/// @brief 
		/// @param borderCode 
		/// @param verbose 
		static void Compute(
			std::string borderCode, 
			bool verbose
		);

	private:
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsLandmask;
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsEdge;
		//--
		std::map<std::string, ign::geometry::GeometryPtr>  _mCountryGeomPtr;
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
		bool _inCountry(std::string const& country, ign::geometry::Point const& pt ) const;

		//--
		void _removeMcoordinate(ign::geometry::LineString & ls) const;

		//--
		bool _isCuttingPoint(
            app::calcul::detail::EdgeCleaningGraphManager & graphManager,
            vertex_descriptor v
        ) const;

		//--
		void _loadGraph(
			app::calcul::detail::EdgeCleaningGraphManager & graphManager,
			bool planarize,
			ign::feature::FeatureFilter filter = ign::feature::FeatureFilter()
		) const;

		//--
		void _loadGraphAndPlanarize(
            app::calcul::detail::EdgeCleaningGraphManager & graphManager/*,
            std::string const& countryCode*/
        ) const;

		//--
		double _getZ( geometry::tools::LengthIndexedLineString const& lengthIndexedLs, ign::geometry::Point const& p) const;

    };

}
}

#endif