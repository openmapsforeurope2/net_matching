#include <app/step/255_GenerateCLinArea2.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

//APP
#include <app/calcul/CLInAreaGenerationOp.h>


namespace app {
	namespace step {

		///
		///
		///
		void GenerateCLinArea2::init()
		{
			addWorkingEntity(EDGE_TABLE_INIT);
		}

		///
		///
		///
		void GenerateCLinArea2::onCompute( bool verbose = false )
		{
			//--
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);

			//--
			app::calcul::CLInAreaGenerationOp::Compute(verbose);
		}

	}
}