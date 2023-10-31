#ifndef _APP_CALCUL_CFEATCONNECTIONOP_H_
#define _APP_CALCUL_CFEATCONNECTIONOP_H_

//SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>

//EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/calcul/matching/detail/LineStringDeformer.h>

// SOCLE
#include <ign/geometry/graph/GeometryGraph.h>

// BOOST
#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>


namespace app{
namespace calcul{

	class CFeatConnectionOp {
		typedef boost::bimap<boost::bimaps::multiset_of<std::string>, boost::bimaps::set_of<std::string>> bimap_t;
        typedef bimap_t::value_type value_type;

	public:

		typedef ign::geometry::graph::GeometryGraph< ign::geometry::graph::PunctualVertexProperties, ign::geometry::graph::LinearEdgeProperties >  GraphType;
		typedef typename GraphType::edge_descriptor edge_descriptor;
		typedef typename GraphType::oriented_edge_descriptor oriented_edge_descriptor;
		typedef typename GraphType::vertex_descriptor vertex_descriptor;
        typedef typename GraphType::edge_iterator edge_iterator;

		/// \brief
		static void computeCp(
			std::string edgeTable, 
            std::string cpTable,
			std::string countryCode, 
			bool verbose
		);

		/// \brief
		static void computeCl(
			std::string edgeTable, 
            std::string clTable,
			std::string countryCode, 
			bool verbose
		);

		/// \brief
		static void computeCpCl(
			std::string edgeTable,
			std::string cpTable,
            std::string clTable,
			std::string countryCode, 
			bool verbose
		);

		/// \brief
		static void computeClImport(
			std::string edgeTable,
			std::string clTable,
			std::string countryCode,
			bool verbose
		);
		

	private:
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsLandmask;
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsEdge;
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsCp;
		//--
		ign::feature::sql::FeatureStorePostgis*            _fsCl;
		//--
		epg::log::EpgLogger*                               _logger;
		//--
		epg::log::ShapeLogger*                             _shapeLogger;
		//--
		std::string                                        _countryCode;
		//--
		bool                                               _verbose;

	private:

		//--
		CFeatConnectionOp( 
            std::string edgeTable, 
            std::string cpTable,
			std::string clTable,
            std::string countryCode, 
            bool verbose 
        );

		//--
		~CFeatConnectionOp();

		//--
		void _init( 
            std::string edgeTable, 
            std::string cpTable,
			std::string clTable
        );

		//--
		void _computeCpDisplacements(std::map<ign::geometry::Point, ign::math::Vec2d> & mDisplacements) const;

		//--
		void _computeClDisplacements(std::map<ign::geometry::Point, ign::math::Vec2d> & mDisplacements) const;

		//--
		void _computeCp();

		//--
		void _computeCl();

		//--
		void _computeCpCl();

		// //--
		// ign::geometry::LineString _mergecl(
		// 	ign::feature::Feature const& refClFeat,
		// 	std::string const& edgeLink,
		// 	std::set<std::string> & sTreatedCl
		// ) const;

		//--
		std::pair<bool, std::string> _getSingleValue(
            std::string edgeLinks,
            std::string countryCodes,
            std::string country
		) const;

		//--
		std::pair<bool, ign::feature::Feature> _getNearestChild(
            ign::geometry::Geometry const& refGeom,
            std::string const& parentFeatureId,
            bimap_t const& mParentChilds
		) const;

		//--
		std::pair<bool, std::string> _getNearestEdge(
			ign::geometry::Geometry const& refGeom,
            std::string countryCodeName,
            std::string edgeLinkName,
            std::string edgeLink,
			std::set<std::string> const& sEdge2Delete
		) const;

		//--
		void _loadEdgeGraph(GraphType & graph) const;

		//--
		std::pair<bool, vertex_descriptor> _getNearestVertex(
            GraphType const& graph, 
            ign::geometry::Point const& pt, 
            double searchDistance
        ) const;

		//--
		void _applyEdgeDisplacement(
            GraphType & graph,
            std::map<ign::geometry::Point, ign::math::Vec2d> const & mReferences,
			std::vector<edge_descriptor> & vDeformedEdges,
            std::set<std::string> & sCollapsedEdges
        ) const;

		//--
		void _applyDisplacement(
			GraphType & graph, 
			std::map< ign::geometry::Point, ign::math::Vec2d > const & mDisplacements,
			epg::calcul::matching::detail::LineStringDeformer const & lineStringDeformer,
			std::vector<edge_descriptor> & vDeformedEdges,
			std::set<std::string> & sCollapsedEdges,
			double influenceDist
		) const;

        //--
        void _persistEdgeDisplacement(
            GraphType & graph,
            std::vector<edge_descriptor> & vDeformedEdges
        ) const;

		//--
		ign::math::Vec2d _computeDisplacement( std::vector< ign::math::Vec2d > const& vVectors ) const;

		void _importCLintoEdgeTable();
    };

}
}

#endif