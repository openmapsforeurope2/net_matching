#ifndef _APP_CALCUL_JUNCTIONMATCHINGOP_H_
#define _APP_CALCUL_JUNCTIONMATCHINGOP_H_

#include <ign/geometry/graph/GeometryGraph.h>

#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/tools/MultiLineStringTool.h>
#include <ign/geometry/graph/GeometryGraph.h>


namespace app{
namespace calcul{

	class JunctionMatchingOp {

	public:

		JunctionMatchingOp(std::string countryCodeDouble, bool verbose = false);
		~JunctionMatchingOp();

		typedef ign::geometry::graph::GeometryGraph< ign::geometry::graph::PunctualVertexProperties, ign::geometry::graph::LinearEdgeProperties >  GraphType;
		typedef typename GraphType::edge_descriptor edge_descriptor;
		typedef typename GraphType::vertex_descriptor vertex_descriptor;
		typedef typename GraphType::oriented_edge_descriptor oriented_edge_descriptor;

		static void MatchJunctions(std::string countryCodeDouble, bool verbose = false);
		static void DisplaceJunctions(std::string countryCodeDouble, bool verbose = false);

		
	private:

		void _init(std::string countryCodeDouble, bool verbose);

		void _matchJunctions();
		void _displaceJunctions();

		void _loadGraphEdges(std::string countryCodeSimple, GraphType& graphEdges);

		bool _IsSimilarIncidentsEdgesOnJunctions(std::set<double> sAnglEdgesJ1, std::set<double> sAnglEdgesJ2);

		double _getAngleEdgeIncident(GraphType& graphEdgCountry, oriented_edge_descriptor& oeit);

		void _getMatchedJunctBest(std::map< vertex_descriptor, vertex_descriptor>& mMatchedJuncRefWithBestJuncMatched, GraphType& graphRef, GraphType& graph2match);

		void _setNewGeomJunction(GraphType& graph, vertex_descriptor& vJunction, ign::geometry::Point ptNewGeomJunction, std::map<std::string, ign::feature::Feature>& mEdgesModifiedGeom);

	private:
		ign::feature::sql::FeatureStorePostgis* _fsEdge;
		ign::feature::sql::FeatureStorePostgis* _fsBoundary;
		//ign::feature::sql::FeatureStorePostgis* _fsBoundarySmoothed;
		ign::feature::sql::FeatureStorePostgis* _fsLandmask;

		epg::log::EpgLogger*                               _logger;
		//--
		epg::log::ShapeLogger*                             _shapeLogger;
		//--
		bool                                               _verbose;

		std::string                                        _countryCodeDouble;

		std::vector<std::string>						   _vCountriesCodeName;

	};

}
}

#endif