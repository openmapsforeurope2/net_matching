#include <app/step/214_UpdateGeomConnectingLines.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <epg/utils/CopyTableUtils.h>

//APP
#include <app/calcul/CFeatGenerationOp.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

//APP
#include <app/calcul/CFeatGenerationOp.h>

namespace app {
	namespace step {

		///
		///
		///
		void UpdateGeomConnectingLines::init()
		{
			addWorkingEntity(CL_TABLE);
		}

		///
		///
		///
		void UpdateGeomConnectingLines::onCompute( bool verbose = false )
		{
			//--
			epg::Context* context = epg::ContextS::getInstance();
			epg::params::EpgParameters & epgParams = context->getEpgParameters();
			std::string idName = epgParams.getValue( ID ).toString();
			std::string geomName = epgParams.getValue( GEOM ).toString();
			//--
			params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string countryCodeW = themeParameters->getParameter(COUNTRY_CODE_W).getValue().toString();
			std::string refTableName = themeParameters->getParameter(CL_TABLE).getValue().toString();

			//--
			epg::utils::CopyTableUtils::copyTable(
				getLastWorkingTableName(CL_TABLE),
				idName,
				geomName,
				ign::geometry::Geometry::GeometryTypeLineString,
				getCurrentWorkingTableName(CL_TABLE),
				"", false, true
			);

			//--
			themeParameters->setParameter(CL_TABLE, ign::data::String(getCurrentWorkingTableName(CL_TABLE)));
			epgParams.setParameter(EDGE_TABLE, ign::data::String(themeParameters->getParameter(EDGE_TABLE_INIT).getValue().toString()));

			//--
			app::calcul::CFeatGenerationOp::UpdateGeomConnectingLines(countryCodeW, verbose);

			//--
			themeParameters->setParameter(CL_TABLE, ign::data::String(refTableName));
			epgParams.setParameter(EDGE_TABLE, ign::data::String(""));
		}

	}
}