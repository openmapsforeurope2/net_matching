#include <app/step/211_MergeConnectingLinesOnBorder.h>

// EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <epg/utils/CopyTableUtils.h>

// APP
#include <app/calcul/CFeatGenerationOp.h>

namespace app {
	namespace step {

		///
		///
		///
		void MergeConnectingLinesOnBorder::init()
		{
			addWorkingEntity(CL_TABLE);
		}

		///
		///
		///
		void MergeConnectingLinesOnBorder::onCompute( bool verbose = false )
		{
			//--
			std::string idName = _epgParams.getValue( ID ).toString();
			std::string geomName = _epgParams.getValue( GEOM ).toString();
			//--
			params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string countryCodeW = _themeParams.getParameter(COUNTRY_CODE_W).getValue().toString();
			std::string clRefTableName = _themeParams.getParameter(CL_TABLE).getValue().toString();

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
			_themeParams.setParameter(CL_TABLE, ign::data::String(getCurrentWorkingTableName(CL_TABLE)));

			//--
			app::calcul::CFeatGenerationOp::MergeConnectingLinesOnBorder(countryCodeW, verbose);

			//--
			_themeParams.setParameter(CL_TABLE, ign::data::String(clRefTableName));
		}

	}
}