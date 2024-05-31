#include <app/step/230_ImportConnectingLines.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

//APP
#include <app/calcul/CFeatGenerationOp.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <epg/utils/CopyTableUtils.h>

//APP
#include <app/calcul/CFeatConnectionOp.h>

namespace app {
	namespace step {

		///
		///
		///
		void ImportConnectingLines::init()
		{
			addWorkingEntity(EDGE_TABLE_INIT);
		}

		///
		///
		///
		void ImportConnectingLines::onCompute( bool verbose = false )
		{
			//--
			params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string countryCodeW = themeParameters->getParameter(COUNTRY_CODE_W).getValue().toString();
			std::string clRefTableName = themeParameters->getParameter(CL_TABLE).getValue().toString();

			//--
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);

			//--
			themeParameters->setParameter(CL_TABLE, ign::data::String(getLastWorkingTableName(CL_TABLE)));

			//--
			app::calcul::CFeatConnectionOp::ComputeClImport(countryCodeW, verbose);

			//--
			themeParameters->setParameter(CL_TABLE, ign::data::String(clRefTableName));
		}

	}
}