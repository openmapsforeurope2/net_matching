#include <app/step/201_CorrectCountryConnectivity.h>

// EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

// APP
 #include <app/calcul/ConnectivityCorrectorOp.h>

namespace app {
	namespace step {

		///
		///
		///
		void CorrectCountryConnectivity::init()
		{
			addWorkingEntity(EDGE_TABLE_INIT);
		}

		///
		///
		///
		void CorrectCountryConnectivity::onCompute(bool verbose = false)
		{
			//--
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);

			//--
			std::string countryCodeW = _themeParams.getParameter(COUNTRY_CODE_W).getValue().toString();

			//--
			app::calcul::ConnectivityCorrectorOp::Compute(countryCodeW, verbose);
		}

	}
}