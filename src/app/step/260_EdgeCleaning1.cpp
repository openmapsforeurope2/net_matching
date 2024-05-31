#include <app/step/260_EdgeCleaning1.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

//APP
#include <app/params/ThemeParameters.h>
#include <app/calcul/EdgeCleaningOp.h>


namespace app {
namespace step {

	///
	///
	///
	void EdgeCleaning1::init()
	{
		addWorkingEntity(EDGE_TABLE_INIT);
	}

	///
	///
	///
	void EdgeCleaning1::onCompute( bool verbose = false )
	{
		//--
		params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
		std::string const countryCodeW = themeParameters->getParameter(COUNTRY_CODE_W).getValue().toString();
		std::string const eclSqlFilter = themeParameters->getValue(ECL_SQL_FILTER).toString();
		std::string const cpRefTableName = themeParameters->getParameter(CP_TABLE).getValue().toString();

		//--
		_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
		ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);

		//--
		themeParameters->setParameter(CP_TABLE, ign::data::String(getLastWorkingTableName(CP_TABLE)));

		//--
		app::calcul::EdgeCleaningOp edgeCleaningOp(countryCodeW, verbose);
        edgeCleaningOp.cleanFaces();
        // edgeCleaningOp.cleanPathsOutOfCountry();
        edgeCleaningOp.cleanParalelleEdges();
        edgeCleaningOp.cleanFacesAndAntennaByCountry(eclSqlFilter);

		//--
		themeParameters->setParameter(CP_TABLE, ign::data::String(cpRefTableName));
	}

}
}
