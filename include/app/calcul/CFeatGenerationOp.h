#ifndef _APP_CALCUL_CFEATGENERATIONOP_H_
#define _APP_CALCUL_CFEATGENERATIONOP_H_

// SOCLE
#include <ign/geometry/graph/GeometryGraph.h>

// EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/tools/MultiLineStringTool.h>
#include <epg/tools/geometry/SegmentIndexedGeometry.h>
#include <epg/sql/tools/IdGeneratorFactory.h>

// OME2
#include <ome2/calcul/utils/AttributeMerger.h>


namespace app{
namespace calcul{

	/// @brief Classe utilitaire dédiée à la génération des connecting lines et des connecting points
	class CFeatGenerationOp {
	private:
		//--
		typedef ign::geometry::graph::GeometryGraph< ign::geometry::graph::PunctualVertexProperties, ign::geometry::graph::LinearEdgeProperties >  GraphType;
		typedef typename GraphType::edge_descriptor                        edge_descriptor;
		typedef typename GraphType::oriented_edge_descriptor               oriented_edge_descriptor;
		typedef typename GraphType::edge_iterator                          edge_iterator;
		typedef typename GraphType::edges_path                             edges_path;
		typedef typename GraphType::edges_path_const_iterator              edges_path_const_iterator;
		typedef typename GraphType::vertex_descriptor                      vertex_descriptor;
		//--
		enum ENDING{
			NONE,
			START,
			END,
			BOTH
		};

	public:

		/// @brief Constructeur
		/// @param borderCode Code pays double (avec '#')
		/// @param verbose Mode verbeux
		CFeatGenerationOp(
			std::string const& borderCode, 
			bool verbose = false
		);

		/// @brief 
		~CFeatGenerationOp();

		/// @brief Generation des connecting lines
		/// @param borderCode Code pays double (avec '#')
		/// @param verbose Mode verbeux
		static void ComputeCL(
			std::string const& borderCode,
			bool verbose = false
		);

		/// @brief Generation des connecting points
		/// @param borderCode Code pays double (avec '#')
		/// @param verbose Mode verbeux
		static void ComputeCP(
			std::string const& borderCode,
			bool verbose = false
		);

		/// @brief Generation des connecting lines pays par pays
		/// @param borderCode 
		/// @param verbose Mode verbeux
		static void GenerateConnectingLinesByCountry(
			std::string const& borderCode,
			bool verbose = false
		);

		/// @brief Fusion des connecting lines projetées sur les frontières
		/// @param borderCode Code pays double (avec '#')
		/// @param verbose Mode verbeux
		static void MergeConnectingLinesOnBorder(
			std::string const& borderCode,
			bool verbose = false
		);

		/// @brief Snap des connecting lines pour éviter des petites discontinuité 
		/// @param borderCode Code pays double (avec '#')
		/// @param verbose Mode verbeux
		static void SnapConnectingLines(
			std::string const& borderCode,
			bool verbose = false
		);

		/// @brief Suppression des connecting lines dont les objets d'origines ne sont pas cohérents (dist et angle)
		/// @param borderCode Code pays double (avec '#')
		/// @param verbose Mode verbeux
		static void DeleteConnectingLines(
			std::string const& borderCode,
			bool verbose = false
		);

		/// @brief Calcul de la géométrie des connecting lines par interpolation des géométries d'origine
		/// @param borderCode Code pays double (avec '#')
		/// @param verbose Mode verbeux
		static void UpdateGeomConnectingLines(
			std::string const& borderCode,
			bool verbose = false
		);

	private:

		//--
		ign::feature::sql::FeatureStorePostgis*  _fsEdge;
		//--
		ign::feature::sql::FeatureStorePostgis*  _fsBoundary;
		//--
		ign::feature::sql::FeatureStorePostgis*  _fsCP;
		//--
		ign::feature::sql::FeatureStorePostgis*  _fsCL;
		//--
		std::string                              _sqlFilterForCpGeneration;
		//--
		epg::log::EpgLogger*                     _logger;
		//--
		epg::log::ShapeLogger*                   _shapeLogger;
		//--
		std::string                              _borderCode;
		//--
		std::vector<std::string>                 _vCountryCode;
		//--
		epg::sql::tools::IdGeneratorInterfacePtr _idGeneratorCP;
		//--
		epg::sql::tools::IdGeneratorInterfacePtr _idGeneratorCL;
		//--
		ome2::calcul::utils::AttributeMerger     _attrMergerOnBorder;
		//--
		epg::tools::MultiLineStringTool*         _mlsBorderSmoothed;
		//--
		std::set<std::string>                    _sFormOfWayException;
		//--
		std::string                              _tagFromCl;
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
		void _init();

		//--
		void _getBorderCutByAngle(
			ign::geometry::LineString const& lsBorder,
			std::vector<ign::geometry::LineString> & vLsBorderCutByAngle,
			double angleMaxToCutBorder
		) const;

		//--
		void _getCLfromBorder(
			ign::geometry::LineString const& lsBorder,
			ign::geometry::Geometry const& buffBorder
		) const;

		//--
		double _getAngleEdgeWithBorder(
			ign::geometry::LineString const& lsEdge,
			ign::geometry::LineString const& lsBorder
		) const;

		//--
		double _getAngle(
			ign::geometry::Point const& ptSource1,
			ign::geometry::Point const& ptTarget1,
			ign::geometry::Point const& ptSource2,
			ign::geometry::Point const& ptTarget2
		) const;

		//--
		ign::geometry::LineString _getGeomCL(
			epg::tools::MultiLineStringTool & mslBorder,
			ign::geometry::LineString const& lsStart2EndToPrject,
			double distMaxBorder,
			double snapOnVertexBorder
		) const;

		//--
		ign::geometry::LineString _getLinestring(
			std::vector<ign::geometry::LineString> const& subEdges,
			size_t idFirst,
			size_t idLast
		) const;

		// //--
		// void _addToUndershootNearBorder(
		// 	ign::geometry::LineString const& lsBorder
		// ) const;

		//--
		void _getCPfromCl() const;

		//--
		ign::geometry::MultiLineString _getBorderGeom() const;

		//--
		bool _isDangle(
			ign::feature::Feature const& fEdge,
			CFeatGenerationOp::ENDING ending
		) const;

		//--
		void _getCPfromBorderUnderShoot(
			ign::geometry::Geometry const& borderGeom,
			epg::tools::geometry::SegmentIndexedGeometry const& segIndexBorder
		) const;

		//--
		void _getCPfromBorder(
			ign::geometry::Geometry const& borderGeom
		) const;

		//--
		void _removeDuplicateCP() const;

		//--
		std::pair<bool, ign::feature::Feature> _hasDuplicateCandidate(
			ign::feature::Feature const& fCp,
			bool fromCl
		) const;

		//--
		void _recordCp(
			ign::geometry::Geometry const& cpGeom,
			ign::feature::Feature const& linkedEdgeFeat,
			std::string tag = ""
		) const;

		// //--
		// void _getCPfromIntersectBorder(
		// 	ign::geometry::LineString const& lsBorder
		// ) const;

		//--
		void _snapCl2Cl(double distMaxClClosest) const;

		//--
		void _snapTo( 
			double distMaxClClosest,
			CFeatGenerationOp::ENDING ending,
			ign::feature::Feature & fCL,
			ign::geometry::LineString & newClGeom,
			std::map<std::string, std::pair<CFeatGenerationOp::ENDING, ign::feature::Feature>> & mFClModified
		) const;

		//--
		std::vector<std::pair<CFeatGenerationOp::ENDING, ign::feature::Feature>> _getClExtremityClose(
			double distMaxClClosest,
			CFeatGenerationOp::ENDING ending,
			ign::feature::Feature const& fClCurr
		) const;

		//--
		void _snapCPNearBy(
			epg::tools::geometry::SegmentIndexedGeometry const& segIndexBorder
		) const;

		//--
		std::pair<bool, ign::feature::Feature> _hasCPfromCL(
			std::list<std::string> const& lCp,
			std::map<std::string, ign::feature::Feature> const& mCPNear
		) const;

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
		std::vector<ign::geometry::LineString> _getIntersectingCls(
			ign::geometry::Geometry const& geom
		) const;

		// //--
		// bool _isEdgeConnected2cl(
		// 	ign::geometry::Geometry const& geomObjNearCl,
		// 	ign::feature::Feature & fCl2SnapOn,
		// 	double distMinCl
		// ) const;

		// //--
		// void _snapCpOnClNearBy(
		// 	std::map<std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint>> & mClSplitedByCp
		// ) const;

		// //--
		// void _cutClByCp(
		// 	std::map<std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint>> const& mClSplittedByCp
		// ) const;

		//--
		void _cutClByCp() const;

		//--
		ign::geometry::Point _getClosestGeometry(
			ign::geometry::MultiPoint const& mlp,
			ign::geometry::LineString const& ls,
			ign::geometry::Point const& pt
		) const;

		//--
		bool _getNearestCP(
			ign::feature::Feature const& fCP,
			double distMergeCP,
			std::map < std::string, ign::feature::Feature> & mCPNear
		) const;

		//--
		void _deleteClByAngleAndDistEdges() const;

		//--
		void _mergeIntersectingClWithGraph(
			double distMaxEdges,
			double snapProjCl2edge
		) const;

		// //--
		// bool _isNextEdgeInAntennas(
		// 	ign::feature::Feature const& fEdge,
		// 	ign::geometry::Point const& ptCurr,
		// 	ign::feature::Feature & edgeNext,
		// 	ign::geometry::Point & ptNext
		// ) const;

		//--
		void _updateGeomCL(double snapOnVertex) const;

		//--
		ign::geometry::LineString _computeMeanGeom(
			ign::geometry::Geometry const& geom1,
			ign::geometry::Geometry const& geom2
		) const;

		//--
		ign::geometry::Geometry* _getGeomProjClOnEdge(
			ign::geometry::LineString const& lsCl,
			ign::geometry::LineString const& lsEdge,
			double snapOnVertex
		) const;

		// //--
		// void _getClDoublonGeom() const;

		//--
		void _loadGraphCL(GraphType & graphCL) const;

		//--
		void _loadGraphEdges(
			std::string const& countryCodeSimple,
			GraphType & graphEdges
		) const;

		//--
		bool _isConnectedEdges(
			GraphType const& graph,
			std::string const& idEdge1,
			std::string const& idEdge2
		) const;

		//--
		std::pair<bool,std::pair<std::string, std::string>> _getClLinkedEdges(
			std::string const& linkedFeatIdName,
			GraphType const& graphCL,
			GraphType::edge_descriptor eCl
		) const;

		//--
		bool _areParallelEdges(
			GraphType const& graphCL,
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
		void _setContinuityCl(GraphType const& graphCL) const;

		//--
		void _deleteCLUnderThreshold() const;

		//--
		double _getLength(
			GraphType const& graph,
			edges_path const& path
		) const;
	
		//--
		edges_path _getPath(
			GraphType const& graph,
			edge_descriptor e
		) const;
 
		//--
		void _mergingEdgesByOrigin(GraphType & graph) const;
	};

}
}

#endif