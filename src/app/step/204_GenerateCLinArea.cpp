#include <app/step/204_GenerateCLinArea.h>

// EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

// APP
#include <app/calcul/CLInAreaGenerationOp.h>

namespace app {
	namespace step {

		///
		///
		///
		void GenerateCLinArea::init()
		{
			addWorkingEntity(EDGE_TABLE_INIT);
		}

		///
		///
		///
		void GenerateCLinArea::onCompute( bool verbose = false )
		{
			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const clMinLength = themeParameters->getValue( CLA_CL_LENGTH_THRESHOLD ).toDouble();
            double const clMinRatio = themeParameters->getValue( CLA_CL_MIN_RATIO_IN_AREA ).toDouble();

			//--
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);

			//--
			app::calcul::CLInAreaGenerationOp::Compute(clMinRatio, clMinLength, verbose);
		}

	}
}