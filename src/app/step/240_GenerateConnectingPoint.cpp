#include <app/step/240_GenerateConnectingPoint.h>

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
			std::string idName = _epgParams.getValue( ID ).toString();
			std::string geomName = _epgParams.getValue( GEOM ).toString();
			std::string edgeRefTableName = _epgParams.getValue( EDGE_TABLE ).toString();
			//--
			params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string countryCodeW = _themeParams.getParameter(COUNTRY_CODE_W).getValue().toString();
			std::string cpRefTableName = _themeParams.getParameter(CP_TABLE).getValue().toString();

			//--
			epg::utils::CopyTableUtils::copyTable(
				getLastWorkingTableName(CP_TABLE),
				idName,
				geomName,
				ign::geometry::Geometry::GeometryTypePoint,
				getCurrentWorkingTableName(CP_TABLE),
				"", false, true
			);
			
			//--
			_themeParams.setParameter(CP_TABLE, ign::data::String(getCurrentWorkingTableName(CP_TABLE)));
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getLastWorkingTableName(EDGE_TABLE_INIT)));

			//--
			app::calcul::CFeatGenerationOp::ComputeCP(countryCodeW, verbose);

			//--
			_themeParams.setParameter(CP_TABLE, ign::data::String(cpRefTableName));
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(edgeRefTableName));
		}

	}
}