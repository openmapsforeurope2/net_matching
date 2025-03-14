#ifndef _APP_CALCUL_CFEATGENERATIONOP_H_
#define _APP_CALCUL_CFEATGENERATIONOP_H_

// SOCLE
#include <ign/geometry/graph/GeometryGraph.h>

// EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/sql/tools/IdGeneratorFactory.h>
#include <epg/tools/MultiLineStringTool.h>
#include <ome2/calcul/utils/AttributeMerger.h>

namespace app{
namespace calcul{

	/// @brief 
	class CFeatGenerationOp {
	private:
		typedef ign::geometry::graph::GeometryGraph< ign::geometry::graph::PunctualVertexProperties, ign::geometry::graph::LinearEdgeProperties >  GraphType;
		typedef typename GraphType::edge_descriptor edge_descriptor;
		typedef typename GraphType::vertex_descriptor vertex_descriptor;

	public:

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		CFeatGenerationOp(std::string countryCodeDouble, bool verbose = false);

		/// @brief 
		~CFeatGenerationOp();

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		static void ComputeCL(std::string countryCodeDouble, bool verbose = false);

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		static void ComputeCP(std::string countryCodeDouble, bool verbose = false);

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		static void GenerateConnectingLinesByCountry(std::string countryCodeDouble, bool verbose = false);

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		static void MergeConnectingLinesOnBorder(std::string countryCodeDouble, bool verbose = false);

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		static void SnapConnectingLines(std::string countryCodeDouble, bool verbose = false);

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		static void DeleteConnectingLines(std::string countryCodeDouble, bool verbose = false);

		/// @brief 
		/// @param countryCodeDouble 
		/// @param verbose 
		static void UpdateGeomConnectingLines(std::string countryCodeDouble, bool verbose = false);

	private:

		//--
		ign::feature::sql::FeatureStorePostgis*  _fsEdge;
		//--
		ign::feature::sql::FeatureStorePostgis*  _fsBoundary;
		//--
		ign::feature::sql::FeatureStorePostgis*  _fsLandmask;
		//--
		ign::feature::sql::FeatureStorePostgis*  _fsCP;
		//--
		ign::feature::sql::FeatureStorePostgis*  _fsCL;
		//--
		std::string                              _reqFilterEdges2generateCF;
		//--
		epg::log::EpgLogger*                     _logger;
		//--
		epg::log::ShapeLogger*                   _shapeLogger;
		//--
		std::string                              _countryCodeDouble;
		//--
		std::vector<std::string>                 _vCountriesCodeName;
		//--
		epg::sql::tools::IdGeneratorInterfacePtr _idGeneratorCP;
		//--
		epg::sql::tools::IdGeneratorInterfacePtr _idGeneratorCL;
		//--
		ome2::calcul::utils::AttributeMerger     _attrMergerOnBorder;
		//--
		epg::tools::MultiLineStringTool*         _mlsBorderSmoothed;
		//--
		std::set<std::string>                    _sFormwayValues4BigDist2Merge;
		//--
		bool                                     _verbose;
		
	private:

		//--
		void _computeCL() const;

		//--
		void _computeCP() const;

		//--
		void _generateConnectingLinesByCountry() const;

		//--
		void _mergeConnectingLinesOnBorder() const;

		//--
		void _snapConnectingLines() const;

		//--
		void _deleteConnectingLines() const;

		//--
		void _updateGeomConnectingLines() const;

		//--
		void _init(std::string countryCodeDouble, bool verbose);

		//--
		void _getBorderCutByAngle(
			ign::geometry::LineString & lsBorder,
			std::vector<ign::geometry::LineString> & vLsBorderCutByAngle,
			double angleMaxToCutBorder
		) const;

		//--
		void _getCLfromBorder(
			ign::geometry::LineString & lsBorder,
			ign::geometry::GeometryPtr & buffBorder,
			double distBuffer,
			double thresholdNoCL,
			double angleMax,
			double ratioInBuff,
			double snapOnVertexBorder
		) const;

		//--
		double _getAngleEdgeWithBorder(
			ign::geometry::LineString & lsEdge,
			ign::geometry::LineString & lsBorder
		) const;

		//--
		void _getGeomCL(
			ign::geometry::LineString & lsCL,
			epg::tools::MultiLineStringTool & mslBorder,
			ign::geometry::LineString& lsStart2EndToPrject,
			double distMaxBorder,
			double snapOnVertexBorder
		) const;

		//--
		void _addToUndershootNearBorder(
			ign::geometry::LineString & lsBorder,
			ign::geometry::GeometryPtr& buffBorder,
			double distUnderShoot
		) const;

		//--
		void _getCPfromIntersectBorder(
			ign::geometry::LineString & lsBorder,
			double distCLIntersected
		) const;

		//--
		void _snapCl2Cl(double distMaxClClosest) const;

		//--
		// @SK : est ce normal que ce ne soit pas la reference des objets passée en paramètre ?
		bool _hasClExtremityClose(
			double distMaxClClosest,
			ign::feature::Feature fClCurr,
			ign::geometry::Point ptClCurr,
			ign::feature::Feature & fCl2snap,
			bool isClosestStartCl2snap
		) const;

		//--
		void _snapCPNearBy(double snapOnVertexBorder) const;

		//--
		bool _areMergeable(
			ign::feature::Feature const& feat1,
			ign::feature::Feature const& feat2,
			double distance
		) const;

		//--
		bool _areDistanceTypeCompatible(
			ign::feature::Feature const& feat1,
			ign::feature::Feature const& feat2,
			double distance
		) const ;

		//--
		bool _areCollinear(
			ign::geometry::LineString const& ls1,
			ign::geometry::LineString const& ls2
		) const;

		//--
		// @SK est ce que les params pourraient etre const ?
		bool _isEdgeConnected2cl(
			ign::geometry::Geometry & geomObjNearCl,
			ign::geometry::Envelope & envArroundGeom,
			ign::feature::Feature & fCl2SnapOn,
			double distMinCl
		) const;

		//--
		void _snapCpOnClNearBy(
			double distCp2snapCl,
			double snapDistOnVertexFromCl,
			std::map< std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint> > & mClSplitedByCp
		) const;

		//--
		void _cutClByCp(
			std::map< std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint> > & mClSplittedByCp
		) const;

		//--
		bool _getNearestCP(
			ign::feature::Feature const& fCP,
			double distMergeCP,
			std::map < std::string, ign::feature::Feature> & mCPNear
		) const;

		//--
		void _deleteClByAngleAndDistEdges(
			double angleMax,
			double distMax,
			double snapOnVertexBorder
		) const;

		//--
		void _mergeIntersectingCL2(
			double distMergeCL,
			double snapOnVertexBorder
		) const;

		//--
		void _mergeIntersectingClWithGraph(
			double distMaxEdges,
			double snapProjCl2edge
		) const;
		
		//--
		// @SK : est ce normal que ce ne soit pas la reference des objets passée en paramètre ?
		bool _getCLToMerge(
			ign::feature::Feature fCL,
			double distMergeCL,
			std::map < std::string, ign::feature::Feature> & mCL2merge,
			std::set<std::string> & sCountryCode
		) const;

		//--
		void _getBorderFromEdge(
			ign::geometry::LineString & lsEdgeOnBorder,
			ign::geometry::LineString & lsBorder
		) const;

		//--
		bool _isNextEdgeInAntennas(
			ign::feature::Feature & fEdge,
			ign::geometry::Point & ptCurr,
			ign::feature::Feature & edgeNext,
			ign::geometry::Point & ptNext
		) const;

		//--
		void _updateGeomCL(double snapOnVertex) const;

		//--
		void _getGeomProjClOnEdge(
			ign::geometry::LineString & lsCl,
			ign::geometry::LineString & lsEdge,
			ign::geometry::LineString & lsprojClOnEdg,
			double snapOnVertex
		) const;

		//--
		void _getClDoublonGeom() const;

		//--
		void _loadGraphCL(GraphType & graphCL) const;

		//--
		void _loadGraphEdges(
			std::string countryCodeSimple,
			GraphType & graphEdges
		) const;

		//--
		bool _isConnectedEdges(
			GraphType & graph,
			std::string const& idEdge1,
			std::string const& idEdge2
		) const;

		//--
		std::pair<bool,std::pair<std::string, std::string>> _getClLinkedEdges(
			std::string const& linkedFeatIdName,
			GraphType & graphCL,
			GraphType::edge_descriptor eCl
		) const;

		//--
		bool _areParallelEdges(
			GraphType & graphCL,
			GraphType::edge_descriptor e1,
			GraphType::edge_descriptor e2
		) const;

		//--
		ign::geometry::Point _getLinkedEdgesConnectingPoint(
			GraphType const& graph,
			std::string const& idEdge1,
			std::string const& idEdge2
		) const;

		//--
		void _setContinuityCl(GraphType & graphCL) const;

		//--
		void _deleteCLUnderThreshold() const;

		//--
		void _getGeomCountry(
			std::string countryCodeSimple,
			ign::geometry::MultiPolygon & geomCountry
		) const;
 
		//--
		void _mergingEdgesByOrigin(GraphType & graph) const;
	};

}
}

#endif