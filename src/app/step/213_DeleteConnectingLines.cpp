#include <app/step/213_DeleteConnectingLines.h>

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
		void DeleteConnectingLines::init()
		{
			addWorkingEntity(CL_TABLE);
		}

		///
		///
		///
		void DeleteConnectingLines::onCompute( bool verbose = false )
		{
			//--
			std::string idName = _epgParams.getValue( ID ).toString();
			std::string geomName = _epgParams.getValue( GEOM ).toString();
			std::string edgeRefTableName = _epgParams.getValue( EDGE_TABLE ).toString();
			//--
			std::string countryCodeW = _themeParams.getValue(COUNTRY_CODE_W).toString();
			std::string clRefTableName = _themeParams.getValue(CL_TABLE).toString();

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
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getLastWorkingTableName(EDGE_TABLE_INIT)));

			//--
			app::calcul::CFeatGenerationOp::DeleteConnectingLines(countryCodeW, verbose);

			//--
			_themeParams.setParameter(CL_TABLE, ign::data::String(clRefTableName));
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(edgeRefTableName));
		}

	}
}