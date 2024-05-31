#include <app/step/240_GenerateConnectingPoint.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <ome2/utils/CopyTableUtils.h>

//APP
#include <app/calcul/CFeatGenerationOp.h>

//EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <epg/utils/CopyTableUtils.h>

//APP
#include <app/calcul/CFeatGenerationOp.h>

namespace app {
	namespace step {

		///
		///
		///
		void GenerateConnectingPoint::init()
		{
			addWorkingEntity(CP_TABLE);
		}

		///
		///
		///
		void GenerateConnectingPoint::onCompute( bool verbose = false )
		{
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string idName = context->getEpgParameters().getValue( ID ).toString();
			std::string geomName = context->getEpgParameters().getValue( GEOM ).toString();
			//--
			params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string countryCodeW = themeParameters->getParameter(COUNTRY_CODE_W).getValue().toString();
			std::string refTableName = themeParameters->getParameter(CP_TABLE).getValue().toString();

			//--
			epg::utils::CopyTableUtils::copyTable(
				getLastWorkingTableName(CP_TABLE),
				idName,
				geomName,
				ign::geometry::Geometry::GeometryTypeLineString,
				getCurrentWorkingTableName(CP_TABLE),
				"", false, true
			);
			themeParameters->setParameter(CP_TABLE, ign::data::String(getCurrentWorkingTableName(CP_TABLE)));

			//--
			app::calcul::CFeatGenerationOp::ComputeCP(countryCodeW, verbose);

			//--
			themeParameters->setParameter(CP_TABLE, ign::data::String(refTableName));
		}

	}
}