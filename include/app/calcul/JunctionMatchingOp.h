#ifndef _APP_CALCUL_JUNCTIONMATCHINGOP_H_
#define _APP_CALCUL_JUNCTIONMATCHINGOP_H_

#include <ign/geometry/graph/GeometryGraph.h>

#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/tools/MultiLineStringTool.h>


namespace app{
namespace calcul{

	class JunctionMatchingOp {

	public:

		JunctionMatchingOp(std::string countryCodeDouble, bool verbose = false);
		~JunctionMatchingOp();


		static void MatchJunctions(std::string countryCodeDouble, bool verbose = false);
		static void DisplaceJunctions(std::string countryCodeDouble, bool verbose = false);

		
	private:

		void _init(std::string countryCodeDouble, bool verbose);

		void _matchJunctions();
		void _displaceJunctions();


	private:
		ign::feature::sql::FeatureStorePostgis* _fsEdge;
		ign::feature::sql::FeatureStorePostgis* _fsBoundary;
		//ign::feature::sql::FeatureStorePostgis* _fsBoundarySmoothed;
		ign::feature::sql::FeatureStorePostgis* _fsLandmask;

		epg::log::EpgLogger*                               _logger;
		//--
		epg::log::ShapeLogger*                             _shapeLogger;
		//--
		bool                                               _verbose;

		std::string                                        _countryCodeDouble;

		std::vector<std::string>						   _vCountriesCodeName;

	};

}
}

#endif