#ifndef _APP_CALCUL_CFEATGENERATIONOP_H_
#define _APP_CALCUL_CFEATGENERATIONOP_H_


#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/sql/tools/IdGeneratorFactory.h>

namespace app{
namespace calcul{

	class CFeatGenerationOp {

	public:

		CFeatGenerationOp(bool verbose = false);
		~CFeatGenerationOp();

		//static void compute(std::string countryCode, bool verbose);

		void computeCL(std::string countryCodeDouble);
		void computeCP(std::string countryCodeDouble);


	private:


		void _init(bool verbose);

		void getCLfromBorder(ign::geometry::LineString & lsBorder, ign::geometry::GeometryPtr& buffBorder,  double distBuffer, double thresholdNoCL, double angleMax, double ratioInBuff, double snapOnVertexBorder);

		double getAngleEdgeWithBorder(ign::geometry::LineString& lsEdge, ign::geometry::LineString& lsBorder);

		void getGeomCL(ign::geometry::LineString& lsCL, ign::geometry::LineString& lsBorder, ign::geometry::Point ptStartToProject, ign::geometry::Point ptEndToProject, double snapOnVertexBorder);
	

		void addToUndershootNearBorder(ign::geometry::LineString & lsBorder, ign::geometry::GeometryPtr& buffBorder, double distUnderShoot);

		void getCPfromIntersectBorder(ign::geometry::LineString & lsBorder, double distCLIntersected);



		bool isEdgeIntersectedPtWithCL(ign::feature::Feature& fEdge, ign::geometry::Point ptIntersectBorder, double distCLIntersected);


		//void mergeCPNearBy(double distMergeCP, double snapOnVertexBorder);
		void snapCPNearBy(std::string countryCodeDouble,double distMergeCP, double snapOnVertexBorder);

		bool getNearestCP(ign::feature::Feature fCP,double distMergeCP, std::map < std::string, ign::feature::Feature>& mCPNear);

		void addFeatAttributeMergingOnBorder(ign::feature::Feature& featMergedAttr, ign::feature::Feature& featAttrToAdd, std::string separator);



		//void mergeCL(double distMergeCL, double snapOnVertexBorder);
		void mergeIntersectingCL(std::string countryCodeDouble, double distMergeCL, double snapOnVertexBorder);

		bool getCLToMerge(ign::feature::Feature fCL, double distMergeCL, std::map < std::string, ign::feature::Feature>& mCL2merge, std::set<std::string>& sCountryCode);


		void getBorderFromEdge(ign::geometry::LineString& lsEdgeOnBorder, ign::geometry::LineString& lsBorder);

		void cleanEdgesOutOfCountry(std::string countryCC);

		void cleanAntennasOutOfCountry(std::string countryCC);

		bool isNextEdgeInAntennas(ign::feature::Feature& fEdge, ign::geometry::Point& ptCurr, ign::feature::Feature&  edgeNext, ign::geometry::Point& ptNext);

		void updateGeomCL(std::string countryCodeDouble, double snapOnVertexBorder);

		void deleteCLUnderThreshold(std::string countryCodeDouble);

		void getGeomCountry(std::string countryCodeSimple, ign::geometry::MultiPolygon& geomCountry);
		
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