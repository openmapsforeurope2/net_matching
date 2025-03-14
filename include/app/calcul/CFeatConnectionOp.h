#ifndef _APP_CALCUL_CFEATCONNECTIONOP_H_
#define _APP_CALCUL_CFEATCONNECTIONOP_H_

// SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>

// EPG
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

	/// @brief 
	class CFeatConnectionOp {
		typedef boost::bimap<boost::bimaps::multiset_of<std::string>, boost::bimaps::set_of<std::string>> bimap_t;
        typedef bimap_t::value_type value_type;

	public:

		typedef ign::geometry::graph::GeometryGraph< ign::geometry::graph::PunctualVertexProperties, ign::geometry::graph::LinearEdgeProperties >  GraphType;
		typedef typename GraphType::edge_descriptor edge_descriptor;
		typedef typename GraphType::oriented_edge_descriptor oriented_edge_descriptor;
		typedef typename GraphType::vertex_descriptor vertex_descriptor;
        typedef typename GraphType::edge_iterator edge_iterator;

		/// @brief 
		/// @param countryCode 
		/// @param verbose 
		CFeatConnectionOp( 
            std::string countryCode, 
            bool verbose 
        );

		/// @brief 
		~CFeatConnectionOp();


		/// @brief 
		/// @param countryCode 
		/// @param verbose 
		static void ComputeCp( 
            std::string countryCode, 
            bool verbose 
        );

		/// @brief 
		/// @param countryCode 
		/// @param verbose 
		static void ComputeCl( 
            std::string countryCode, 
            bool verbose 
        );

		/// @brief 
		/// @param countryCode 
		/// @param verbose 
		static void ComputeClImport( 
            std::string countryCode, 
            bool verbose 
        );

		/// \brief
		void computeCp();

		/// \brief
		void computeCl();

		/// \brief
		void computeCpCl();

		/// \brief
		void computeClImport();


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
		void _init();

		//--
		void _computeCpDisplacements(std::map<ign::geometry::Point, ign::math::Vec2d> & mDisplacements, std::string const& country) const;

		//--
		void _computeClDisplacements(std::map<ign::geometry::Point, ign::math::Vec2d> & mDisplacements, std::string const& country) const;

		//--
		void _computeCp(std::string const& country);

		//--
		void _computeCl(std::string const& country);

		//--
		void _computeCpCl(std::string const& country);

		//--
		void _addDisplacement(
            ign::geometry::Point const& sourcePoint,
            ign::geometry::Point const& targetPoint,
            std::map<ign::geometry::Point, ign::math::Vec2d>& mDisplacements,
			std::map<ign::geometry::Point, ign::geometry::LineString>& mDisplacementCls,
            ign::feature::Feature const& fCl,
            std::string const& linkedFeatureId,
            std::string const& country
        ) const;

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
		void _loadEdgeGraph(GraphType & graph, std::string const& country) const;

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

		//--
		void _importCLintoEdgeTable();
    };

}
}

#endif