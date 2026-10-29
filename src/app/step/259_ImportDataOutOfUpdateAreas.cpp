#include <app/step/259_ImportDataOutOfUpdateAreas.h>

// EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

// APP
#include <app/calcul/ImportDataOutOfUpdateAreasOp.h>


namespace app {
	namespace step {

		///
		///
		///
		void ImportDataOutOfUpdateAreas::init()
		{
			addWorkingEntity(EDGE_TABLE_INIT);
		}

		///
		///
		///
		void ImportDataOutOfUpdateAreas::onCompute( bool verbose = false )
		{
			//--
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true);

			//--
			app::calcul::ImportDataOutOfUpdateAreasOp::Compute(verbose);
		}

	}
}