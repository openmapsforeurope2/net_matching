#include <app/step/215_UpdateGeomByContinuity.h>

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
		void UpdateGeomByContinuity::init()
		{
			addWorkingEntity(CL_TABLE);
		}

		///
		///
		///
		void UpdateGeomByContinuity::onCompute( bool verbose = false )
		{
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string idName = context->getEpgParameters().getValue( ID ).toString();
			std::string geomName = context->getEpgParameters().getValue( GEOM ).toString();
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
			themeParameters->setParameter(CL_TABLE, ign::data::String(getCurrentWorkingTableName(CL_TABLE)));

			//--
			app::calcul::CFeatGenerationOp::UpdateGeomByContinuity(countryCodeW, verbose);

			//--
			themeParameters->setParameter(CL_TABLE, ign::data::String(refTableName));
		}

	}
}