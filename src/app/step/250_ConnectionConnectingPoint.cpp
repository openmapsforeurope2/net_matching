#include <app/step/250_ConnectionConnectingPoint.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <epg/utils/CopyTableUtils.h>

//APP
#include <app/calcul/CFeatConnectionOp.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

//APP
#include <app/calcul/CFeatConnectionOp.h>

namespace app {
	namespace step {

		///
		///
		///
		void ConnectionConnectingPoint::init()
		{
			addWorkingEntity(EDGE_TABLE_INIT);
		}

		///
		///
		///
		void ConnectionConnectingPoint::onCompute( bool verbose = false )
		{
			//--
			params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string countryCodeW = themeParameters->getParameter(COUNTRY_CODE_W).getValue().toString();
			std::string cpRefTableName = themeParameters->getParameter(CP_TABLE).getValue().toString();

			//--
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);

			//--
			themeParameters->setParameter(CP_TABLE, ign::data::String(getLastWorkingTableName(CP_TABLE)));

			//--
			app::calcul::CFeatConnectionOp::ComputeCp(countryCodeW, verbose);

			//--
			themeParameters->setParameter(CP_TABLE, ign::data::String(cpRefTableName));
		}

	}
}