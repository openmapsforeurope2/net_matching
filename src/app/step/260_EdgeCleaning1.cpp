#include <app/step/260_EdgeCleaning1.h>

// EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

// APP
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
		std::string const countryCodeW = _themeParams.getValue(COUNTRY_CODE_W).toString();
		std::string const eclSqlFilter = _themeParams.getValue(ECL_SQL_FILTER).toString();
		std::string const cpRefTableName = _themeParams.getValue(CP_TABLE).toString();

		//--
		_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
		ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true);

		//--
		_themeParams.setParameter(CP_TABLE, ign::data::String(getLastWorkingTableName(CP_TABLE)));

		//--
		app::calcul::EdgeCleaningOp edgeCleaningOp(countryCodeW, verbose);
		edgeCleaningOp.cleanParalelleEdges(); /*a faire avant le reste*/
		edgeCleaningOp.cleanFaces();
		edgeCleaningOp.cleanFacesAndAntennaByCountry(eclSqlFilter, false /*tagTreatedFeatures*/);

		//--
		_themeParams.setParameter(CP_TABLE, ign::data::String(cpRefTableName));
	}

}
}
