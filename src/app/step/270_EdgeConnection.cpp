#include <app/step/270_EdgeConnection.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

//APP
#include <app/params/ThemeParameters.h>
#include <app/calcul/EdgeConnectorOp.h>


namespace app {
namespace step {

	///
	///
	///
	void EdgeConnection::init()
	{
		addWorkingEntity(EDGE_TABLE_INIT);
	}

	///
	///
	///
	void EdgeConnection::onCompute( bool verbose = false )
	{
		//--
		std::string const countryCodeW = _themeParams.getParameter(COUNTRY_CODE_W).getValue().toString();

		//--
		_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
		ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);

		//--
		app::calcul::EdgeConnectorOp::Compute(countryCodeW, verbose);
	}

}
}
