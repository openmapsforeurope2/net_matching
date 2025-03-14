#ifndef _APP_CALCUL_JUNCTIONMATCHINGOP_H_
#define _APP_CALCUL_JUNCTIONMATCHINGOP_H_

//SOCLE
#include <ign/geometry/graph/GeometryGraph.h>

//EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/tools/MultiLineStringTool.h>

//APP
#include <app/calcul/detail/EdgeCleaningGraphManager.h>


namespace app{
namespace calcul{

	/// @brief 
	class JunctionMatchingOp {

	public:

		typedef app::calcul::detail::EdgeCleaningGraphManager::GraphType   GraphType;
		typedef typename GraphType::edge_descriptor                        edge_descriptor;
		typedef typename GraphType::vertex_descriptor                      vertex_descriptor;
		typedef typename GraphType::oriented_edge_descriptor               oriented_edge_descriptor;
		typedef app::calcul::detail::OriginEdgeProperties                  OriginEdgeProperties;

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		JunctionMatchingOp(std::string const& countryCodeDouble, bool verbose = false);

		/// @brief 
		~JunctionMatchingOp();

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		static void MatchJunctions(std::string const& countryCodeDouble, bool verbose = false);

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		static void DisplaceJunctions(std::string const& countryCodeDouble, bool verbose = false);


	private:
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsEdge;
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsBoundary;
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsLandmask;
		//--
		epg::log::EpgLogger*                               _logger;
		//--
		epg::log::ShapeLogger*                             _shapeLogger;
		//--
		bool                                               _verbose;
		//--
		std::string                                        _countryCodeDouble;
		//--
		std::vector<std::string>						   _vCountriesCodeName;

		
	private:

		//--
		void _init(std::string const& countryCodeDouble, bool verbose);

		//--
		void _matchJunctions() const;

		//--
		void _displaceJunctions() const;

		//--
		void _loadGraph(
			std::string const& countryCodeSimple,
			app::calcul::detail::EdgeCleaningGraphManager & graphManager
		) const;

		bool _isFictitious(
			vertex_descriptor vJunction,
			app::calcul::detail::EdgeCleaningGraphManager const& graphManager
		) const;

		//--
		bool _IsSimilarIncidentsEdgesOnJunctions(
			std::set<double> const& sAnglEdgesJ1,
			std::set<double> const& sAnglEdgesJ2
		) const;

		//--
		void _getMatchedJunctBest(
			std::map< vertex_descriptor, vertex_descriptor> & mMatchedJuncRefWithBestJuncMatched,
			GraphType const& graphRef,
			GraphType const& graph2match
		) const;

		//--
		void _setNewGeomJunction(
			GraphType const& graph,
			vertex_descriptor vJunction,
			ign::geometry::Point const& ptNewGeomJunction,
			std::map<std::string, ign::feature::Feature> & mEdgesModifiedGeom
		) const;

	};

}
}

#endif