#include <app/step/201_FillFictitiousField.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>


//APP
 #include <app/calcul/FillFictitiousFieldOp.h>

namespace app {
	namespace step {

		///
		///
		///
		void FillFictitiousField::init()
		{
			addWorkingEntity(EDGE_TABLE_INIT);
		}

		///
		///
		///
		void FillFictitiousField::onCompute(bool verbose = false)
		{
			//--
			std::string countryCodeW = _themeParams.getParameter(COUNTRY_CODE_W).getValue().toString();
			std::string clRefTableName = _themeParams.getParameter(CL_TABLE).getValue().toString();

			//--
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);

			app::calcul::FillFictitiousFieldOp::Compute(countryCodeW, verbose);

		}

	}
}