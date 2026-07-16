#ifndef _APP_CALCUL_JUNCTIONMATCHINGOP_H_
#define _APP_CALCUL_JUNCTIONMATCHINGOP_H_

// SOCLE
#include <ign/geometry/graph/GeometryGraph.h>

// EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/tools/MultiLineStringTool.h>

// APP
#include <app/calcul/detail/EdgeCleaningGraphManager.h>


namespace app{
namespace calcul{

	/// @brief Classe utilisée pour l'appairage des carrefours
	class JunctionMatchingOp {

	public:

		typedef app::calcul::detail::EdgeCleaningGraphManager::GraphType   GraphType;
		typedef typename GraphType::edge_descriptor                        edge_descriptor;
		typedef typename GraphType::vertex_descriptor                      vertex_descriptor;
		typedef typename GraphType::oriented_edge_descriptor               oriented_edge_descriptor;
		typedef app::calcul::detail::OriginEdgeProperties                  OriginEdgeProperties;

		/// @brief Constructeur
		/// @param borderCode Code pays double
		/// @param verbose Mode verbeux
		JunctionMatchingOp(std::string const& borderCode, bool verbose = false);

		/// @brief Destructeur
		~JunctionMatchingOp();

		/// @brief Appairage des carrefours. Si l'un des deux carrefours appairés est fictif, sa position est conservée 
		/// et c'est l'autre  carrefour qui est déplacé. Sinon les deux carrefours A et B sont déplacés au milieu 
		/// du segment [AB].
		/// @param borderCode Code pays double
		/// @param verbose Mode verbeux
		static void Compute(std::string const& borderCode, bool verbose = false);


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
		std::string                                        _borderCode;
		//--
		std::vector<std::string>						   _vCountry;

		
	private:

		//--
		void _init();

		//--
		void _matchJunctions() const;

		//--
		void _displaceJunctions() const;

		//--
		void _loadGraph(
			std::string const& countryCode,
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
		void _computeOrientedMatching(
			std::map< vertex_descriptor,vertex_descriptor> & mJunct1Junct2,
			GraphType const& graph1,
			GraphType const& graph2
		) const;

		//--
		void _moveVertex(
			GraphType const& graph,
			vertex_descriptor v,
			ign::geometry::Point const& newVertexGeom,
			std::map<std::string, ign::feature::Feature> & mEdges2Modify
		) const;

	};

}
}

#endif