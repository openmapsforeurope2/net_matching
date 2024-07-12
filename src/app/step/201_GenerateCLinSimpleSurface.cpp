#include <app/step/201_GenerateCLinSimpleSurface.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

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
		void GenerateCLinSimpleSurface::init()
		{
			addWorkingEntity(EDGE_TABLE_INIT);
		}

		///
		///
		///
		void GenerateCLinSimpleSurface::onCompute( bool verbose = false )
		{
			//--
			std::string countryCodeW = _themeParams.getParameter(COUNTRY_CODE_W).getValue().toString();
			std::string clRefTableName = _themeParams.getParameter(CL_TABLE).getValue().toString();

			//--
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true, true);



			//--
			//_themeParams.setParameter(CL_TABLE, ign::data::String(getLastWorkingTableName(CL_TABLE)));

			//--
			//app::calcul::CFeatConnectionOp::ComputeCl(countryCodeW, verbose);

			//--
			//_themeParams.setParameter(CL_TABLE, ign::data::String(clRefTableName));
		}

	}
}