#ifndef _APP_CALCUL_CLINAREAGENERATIONOP_H_
#define _APP_CALCUL_CLINAREAGENERATIONOP_H_

#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>

// SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>

//APP
#include <app/calcul/detail/EdgeCleaningGraphManager.h>

//EPG
#include <ome2/calcul/utils/AttributeMerger.h>


namespace app{
namespace calcul{
namespace detail{

	enum ENDING{
		START,
		END
	};

    struct IncidentFeature {
        /// \brief
        IncidentFeature(
			std::string originId_,
			ENDING ending_
		):
			originId(originId_),
			ending(ending_)
		{};

        /// \brief
        ~IncidentFeature(){};

        //--
        ENDING		          ending;
        std::string           originId;
		ign::geometry::Point  ptTarget;
    };
}
}
}

namespace app{
namespace calcul{

	class CLInAreaGenerationOp {

	public:

		typedef app::calcul::detail::EdgeCleaningGraphManager::GraphType   GraphType;
		typedef typename GraphType::edge_descriptor                        edge_descriptor;
		typedef typename GraphType::edge_iterator                          edge_iterator;
		typedef typename GraphType::oriented_edge_descriptor               oriented_edge_descriptor;
		typedef typename GraphType::face_descriptor                        face_descriptor;
		typedef typename GraphType::face_iterator                          face_iterator;
		typedef typename GraphType::vertex_descriptor                      vertex_descriptor;
		typedef typename GraphType::vertex_iterator                        vertex_iterator;
		typedef app::calcul::detail::OriginEdgeProperties                  OriginEdgeProperties;
		typedef typename GraphType::edges_path                             edges_path;
		typedef typename GraphType::edges_path_const_iterator              edges_path_const_iterator;

		typedef std::multimap<std::string, detail::IncidentFeature>::const_iterator  m_iterator;

		/// \brief
		static void compute(
			bool verbose
		);

	private:
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsEdge;
		//--
		epg::log::EpgLogger*                               _logger;
		//--
		epg::log::ShapeLogger*                             _shapeLogger;
		//--
		bool                                               _verbose;
		//--
		ome2::calcul::utils::AttributeMerger               _attrMerger;

	private:

		//--
		CLInAreaGenerationOp(bool verbose = false);

		//--
		~CLInAreaGenerationOp();

		//--
		void _init();

		//--
		void _compute() const;

		//--
		void _createCLOnFaces(
            detail::EdgeCleaningGraphManager const& graphManager,
            std::map<std::string, std::set<edge_descriptor>> & mFeatMergedEdges,
            std::multimap<std::string, detail::IncidentFeature> & mmIncidentFeatures
        ) const;

		//--
		void _createCLOnOverlappingEdges(
            GraphType const& graph,
            std::map<std::string, std::set<edge_descriptor>> & mFeatMergedEdges,
            std::multimap<std::string, detail::IncidentFeature> & mmIncidentFeatures
        ) const;

		//--
		void _cleanOverLappingCl() const;

		//--
		void _addIncidentFeatures(
            GraphType const& graph,
            vertex_descriptor v,
            std::set<edge_descriptor> const& sEdges,
            std::multimap<std::string, detail::IncidentFeature> & mmIncidentFeatures
        ) const;

		//--
		int _getIndex(
			ign::geometry::LineString const& ls,
            ign::geometry::Point const& pt
		) const;

		//--
		void _loadGraph(
			app::calcul::detail::EdgeCleaningGraphManager & graphManager,
			ign::feature::FeatureFilter filter = ign::feature::FeatureFilter(),
			bool planarize = true
		) const;

		//--
		void _mergeByWTag() const;

		//--
		bool _pathsGeomAreEqual(
			ign::geometry::Polygon const& poly,
            ign::geometry::LineString & path1geom,
            ign::geometry::LineString & path2geom,
            double maxWidth,
            bool useHausdorff = false
        ) const;

		//--
		ign::geometry::LineString _computeMeanPath(
            ign::geometry::LineString const& path1geom,
            ign::geometry::LineString const& path2geom
		) const;

		//--
		std::map<double, std::string> _getOriginEdges(
			GraphType const& graph,
			std::string const& country,
			std::list<oriented_edge_descriptor> const& path,
			std::map<std::string, std::set<edge_descriptor>> & mFeatMergedEdges,
			std::map<double, std::vector<detail::IncidentFeature>> & mIncidentFeatures
		) const;

		//--
		std::set<std::string> _mergeFacePaths(
            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> & vpCountryEdges
        ) const;

		//--
		std::list<app::calcul::detail::EdgeCleaningGraphManager::GraphType::oriented_edge_descriptor> _getExteriorRingEdges(
            GraphType const& graph, 
            face_descriptor fd
        ) const;

		//--
		bool _getFacePaths(
            detail::EdgeCleaningGraphManager const& graphManager, 
            face_descriptor fd, 
            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> & vpCountryEdges
        ) const;

		//--
		ign::geometry::LineString _convertPathToLineString(
            GraphType const& graph,
			std::string const& country,
            std::list<oriented_edge_descriptor> const& path
        ) const;

		//--
		bool _isFictitious(
            GraphType const& graph,
            std::string const& country,
            std::list<oriented_edge_descriptor> const& path
        ) const;

		//--
		std::string _getOrigin(
            GraphType const& graph,
            std::string const& country,
            edge_descriptor e
        ) const;

		//--
		void _setZ(
            ign::geometry::Point & pt,
            ign::geometry::LineString const& refGeom
        ) const;

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

	};

}
}

#endif