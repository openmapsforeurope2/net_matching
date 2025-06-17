#include <app/step/202_FillFictitiousField.h>

// EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

// APP
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
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);

			app::calcul::FillFictitiousFieldOp::Compute(verbose);
		}

	}
}