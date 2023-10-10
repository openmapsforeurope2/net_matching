#ifndef _APP_CALCUL_EDGECLEANINGOP_H_
#define _APP_CALCUL_EDGECLEANINGOP_H_

//SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>

//EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>

//APP
#include <app/calcul/detail/EdgeCleaningGraphManager.h>

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
		typedef app::calcul::detail::EdgeCleaningEdge                      EdgeCleaningEdge;

		/// \brief
		static void clean(
			std::string edgeTable,
			std::string borderCode,
			bool verbose
		);

	private:
		//--
		ign::feature::sql::FeatureStorePostgis*              _fsEdge;
		//--
		ign::feature::sql::FeatureStorePostgis*              _fsCl;
		//--
		std::map<std::string, ign::geometry::GeometryPtr>    _mCountryGeomPtr;
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
		void _loadGraph(app::calcul::detail::EdgeCleaningGraphManager & graphManager, bool planarize = false) const;

		//--
		double _getLength( ign::geometry::Geometry const& geom ) const;

		//--
		double _getRatio(GraphType const& graph, std::string country, std::list<edge_descriptor> const& lEdges) const;

		//--
		void _removeEdges(GraphType const& graph, std::list<edge_descriptor> const& lEdges) const;

		//--
		void _addLengths(
			std::string country,
			ign::geometry::LineString const& ls,
			double & lengthInCountry,
			double & length
		) const;

		//--
		void _clean() const;
    };

}
}

#endif