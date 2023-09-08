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

		CFeatGenerationOp(bool verbose = false);
		~CFeatGenerationOp();

		typedef ign::geometry::graph::GeometryGraph< ign::geometry::graph::PunctualVertexProperties, ign::geometry::graph::LinearEdgeProperties >  GraphType;
		typedef typename GraphType::edge_descriptor edge_descriptor;
		typedef typename GraphType::vertex_descriptor vertex_descriptor;

		//static void compute(std::string countryCode, bool verbose);
		void computeCL(std::string countryCodeDouble);
		void computeCP(std::string countryCodeDouble);


	private:


		void _init(bool verbose);

		void _getCLfromBorder(ign::geometry::LineString & lsBorder, ign::geometry::GeometryPtr& buffBorder,  double distBuffer, double thresholdNoCL, double angleMax, double ratioInBuff, double snapOnVertexBorder);

		double _getAngleEdgeWithBorder(ign::geometry::LineString& lsEdge, ign::geometry::LineString& lsBorder);

		void _getGeomCL(ign::geometry::LineString& lsCL, ign::geometry::LineString& lsBorder, ign::geometry::Point ptStartToProject, ign::geometry::Point ptEndToProject, double snapOnVertexBorder);
	

		void _addToUndershootNearBorder(ign::geometry::LineString & lsBorder, ign::geometry::GeometryPtr& buffBorder, double distUnderShoot);

		void _getCPfromIntersectBorder(ign::geometry::LineString & lsBorder, double distCLIntersected);



		bool _isEdgeIntersectedPtWithCL(ign::feature::Feature& fEdge, ign::geometry::Point ptIntersectBorder, double distCLIntersected);


		//void mergeCPNearBy(double distMergeCP, double snapOnVertexBorder);
		void _snapCPNearBy(std::string countryCodeDouble,double distMergeCP, double snapOnVertexBorder);

		bool _getNearestCP(ign::feature::Feature fCP,double distMergeCP, std::map < std::string, ign::feature::Feature>& mCPNear);

		void _addFeatAttributeMergingOnBorder(ign::feature::Feature& featMergedAttr, ign::feature::Feature& featAttrToAdd, std::string separator);

		void _deleteClByAngleEdges(std::string countryCodeDouble, double angleMax, double snapOnVertexBorder);

		//void mergeCL(double distMergeCL, double snapOnVertexBorder);
		void _mergeIntersectingCL(std::string countryCodeDouble, double distMergeCL,  double snapOnVertexBorder);

		bool _getCLToMerge(ign::feature::Feature fCL, double distMergeCL, std::map < std::string, ign::feature::Feature>& mCL2merge, std::set<std::string>& sCountryCode);


		void _getBorderFromEdge(ign::geometry::LineString& lsEdgeOnBorder, ign::geometry::LineString& lsBorder);

		//void _cleanEdgesOutOfCountry(std::string countryCC);

		//void _cleanAntennasOutOfCountry(std::string countryCC);

		bool _isNextEdgeInAntennas(ign::feature::Feature& fEdge, ign::geometry::Point& ptCurr, ign::feature::Feature&  edgeNext, ign::geometry::Point& ptNext);

		void _updateGeomCL(std::string countryCodeDouble, double snapOnVertex);

		void _getGeomProjClOnEdge(ign::geometry::LineString& lsCl, ign::geometry::LineString& lsEdge, ign::geometry::LineString& lsprojClOnEdg, double snapOnVertex);

		void _getClDoublonGeom(std::string countryCodeDouble);

		void _loadGraphCL(std::string countryCodeDouble, GraphType& graphCL);

		void _setContinuityCl(GraphType& graphCL);
		//void _getClContinuity(std::map<std::string, std::vector<std::pair<std::string, bool>>>& mClConnect);
		//void _getEdgesConnectedOnPoint(ign::geometry::Point ptConnect, std::vector<std::pair<std::string, bool>>& vEdgesConnection);

		void _deleteCLUnderThreshold(std::string countryCodeDouble);

		void _getGeomCountry(std::string countryCodeSimple, ign::geometry::MultiPolygon& geomCountry);
		
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

		//std::string                                        _countryCodeDouble;

		epg::sql::tools::IdGeneratorInterfacePtr _idGeneratorCP;
		epg::sql::tools::IdGeneratorInterfacePtr _idGeneratorCL;

		std::set<std::string> _sAttrNameToConcat;
		std::set<std::string> _sAttrNameW;
		

	};

}
}

#endif