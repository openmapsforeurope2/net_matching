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
            double const clMinLength = _themeParams.getValue( CLA_CL_LENGTH_THRESHOLD ).toDouble();
            double const clMinRatio = _themeParams.getValue( CLA_CL_MIN_RATIO_IN_AREA ).toDouble();
			std::string const countryCodeW = _themeParams.getValue(COUNTRY_CODE_W).toString();

			//--
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true);

			//--
			app::calcul::CLInAreaGenerationOp::Compute(countryCodeW, clMinRatio, clMinLength, verbose);
		}

	}
}