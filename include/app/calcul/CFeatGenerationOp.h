#ifndef _APP_CALCUL_CFEATGENERATIONOP_H_
#define _APP_CALCUL_CFEATGENERATIONOP_H_

#include <ign/geometry/graph/GeometryGraph.h>

#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/sql/tools/IdGeneratorFactory.h>

namespace app{
namespace calcul{

	class CFeatGenerationOp {

	public:

		CFeatGenerationOp(std::string countryCodeDouble, bool verbose = false);
		~CFeatGenerationOp();

		typedef ign::geometry::graph::GeometryGraph< ign::geometry::graph::PunctualVertexProperties, ign::geometry::graph::LinearEdgeProperties >  GraphType;
		typedef typename GraphType::edge_descriptor edge_descriptor;
		typedef typename GraphType::vertex_descriptor vertex_descriptor;

		//static void compute(std::string countryCode, bool verbose);
		void computeCL();
		void computeCP();


	private:


		void _init(std::string countryCodeDouble, bool verbose);

		void _getCLfromBorder(ign::geometry::LineString & lsBorder, ign::geometry::GeometryPtr& buffBorder,  double distBuffer, double thresholdNoCL, double angleMax, double ratioInBuff, double snapOnVertexBorder);

		double _getAngleEdgeWithBorder(ign::geometry::LineString& lsEdge, ign::geometry::LineString& lsBorder);

		void _getGeomCL(ign::geometry::LineString& lsCL, ign::geometry::LineString& lsBorder, ign::geometry::Point ptStartToProject, ign::geometry::Point ptEndToProject, double snapOnVertexBorder);
	

		void _addToUndershootNearBorder(ign::geometry::LineString & lsBorder, ign::geometry::GeometryPtr& buffBorder, double distUnderShoot);

		void _getCPfromIntersectBorder(ign::geometry::LineString & lsBorder, double distCLIntersected);

		void _snapCl2Cl(double distMaxClClosest);

		bool _hasClExtremityClose(double distMaxClClosest, ign::feature::Feature fClCurr, ign::geometry::Point ptClCurr, ign::feature::Feature& fCl2snap, bool& isClosestStartCl2snap);

		//bool _isEdgeIntersectedPtWithCL(ign::feature::Feature& fEdge, ign::geometry::Point ptIntersectBorder, double distCLIntersected);


		//void mergeCPNearBy(double distMergeCP, double snapOnVertexBorder);
		void _snapCPNearBy(double distMergeCP, double snapOnVertexBorder);
		bool _areCollinear(
			ign::geometry::LineString const& ls1,
			ign::geometry::LineString const& ls2
		) const;

		void _snapCPNearBy2(double distMergeCP, double snapOnVertexBorder);

		void _snapCpOnClNearBy(double distCp2snapCl, double snapDistOnVertexFromCl, std::map< std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint> > & mClSplitedByCp);

		void _cutClByCp(std::map< std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint> > & mClSplittedByCp);

		bool _getNearestCP(ign::feature::Feature fCP,double distMergeCP, std::map < std::string, ign::feature::Feature>& mCPNear);

		void _addFeatAttributeMergingOnBorder(ign::feature::Feature& featMergedAttr, ign::feature::Feature& featAttrToAdd, std::string separator);

		void _deleteClByAngleAndDistEdges(double angleMax, double distMax, double snapOnVertexBorder);

		//void mergeCL(double distMergeCL, double snapOnVertexBorder);
		void _mergeIntersectingCL2( double distMergeCL, double snapOnVertexBorder);
		void _mergeIntersectingClWithGraph(double distMaxEdges, double snapProjCl2edge);
		
		bool _getCLToMerge(ign::feature::Feature fCL, double distMergeCL, std::map < std::string, ign::feature::Feature>& mCL2merge, std::set<std::string>& sCountryCode);


		void _getBorderFromEdge(ign::geometry::LineString& lsEdgeOnBorder, ign::geometry::LineString& lsBorder);

		//void _cleanEdgesOutOfCountry(std::string countryCC);

		//void _cleanAntennasOutOfCountry(std::string countryCC);

		bool _isNextEdgeInAntennas(ign::feature::Feature& fEdge, ign::geometry::Point& ptCurr, ign::feature::Feature&  edgeNext, ign::geometry::Point& ptNext);

		void _updateGeomCL(double snapOnVertex);

		void _getGeomProjClOnEdge(ign::geometry::LineString& lsCl, ign::geometry::LineString& lsEdge, ign::geometry::LineString& lsprojClOnEdg, double snapOnVertex);

		void _getClDoublonGeom();

		void _loadGraphCL(GraphType& graphCL);
		void _loadGraphEdges(std::string countryCodeSimple, GraphType& graphEdges);

		bool _isConnectedEdges(GraphType& graph, std::string idEdge1, std::string idEdge2);

		std::pair<std::string, std::string> _getClLinkedEdges( std::string const& linkedFeatIdName, GraphType& graphCL, GraphType::edge_descriptor eCl );

		bool _areParallelEdges( GraphType& graphCL, GraphType::edge_descriptor e1,  GraphType::edge_descriptor e2 );

		ign::geometry::Point _getLinkedEdgesConnectingPoint( GraphType const& graph, std::string const& idEdge1, std::string const& idEdge2 );

		void _setContinuityCl(GraphType& graphCL);
		//void _getClContinuity(std::map<std::string, std::vector<std::pair<std::string, bool>>>& mClConnect);
		//void _getEdgesConnectedOnPoint(ign::geometry::Point ptConnect, std::vector<std::pair<std::string, bool>>& vEdgesConnection);

		void _deleteCLUnderThreshold();

		void _getGeomCountry(std::string countryCodeSimple, ign::geometry::MultiPolygon& geomCountry);
 
		void _mergingClUnderThreshold(double threshold);
		//void _mergingClByLength(GraphType& graph, int threshold);

		//GraphType::edge_descriptor _mergeEdgesCl(GraphType & graph, GraphType::edge_descriptor d, GraphType::edges_path & path);
		//GraphType::edge_descriptor _switchEdge(GraphType& graph, GraphType::edge_descriptor oldEdge, GraphType::vertex_descriptor vSource, GraphType::vertex_descriptor vTarget, ign::geometry::LineString const& geom);
		
	private:
		ign::feature::sql::FeatureStorePostgis* _fsEdge;
		ign::feature::sql::FeatureStorePostgis* _fsBoundary;
		ign::feature::sql::FeatureStorePostgis* _fsLandmask;
		ign::feature::sql::FeatureStorePostgis* _fsCP;
		ign::feature::sql::FeatureStorePostgis* _fsCL;

		//ign::feature::FeatureFilter _filterEdges2generateCF;
		std::string _reqFilterEdges2generateCF;

		epg::log::EpgLogger*                               _logger;
		//--
		epg::log::ShapeLogger*                             _shapeLogger;
		//--
		bool                                               _verbose;

		std::string                                        _countryCodeDouble;

		std::vector<std::string>						   _vCountriesCodeName;

		epg::sql::tools::IdGeneratorInterfacePtr _idGeneratorCP;
		epg::sql::tools::IdGeneratorInterfacePtr _idGeneratorCL;

		std::set<std::string> _sAttrNameToConcat;
		std::set<std::string> _sAttrNameW;
		

	};

}
}

#endif