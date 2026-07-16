#include <app/step/240_GenerateConnectingPoint.h>

// EPG
#include <epg/Context.h>
#include <epg/log/ScopeLogger.h>
#include <epg/utils/CopyTableUtils.h>
#include <ome2/utils/CopyTableUtils.h>

// APP
#include <app/calcul/CFeatGenerationOp.h>


namespace app {
	namespace step {

		///
		///
		///
		void GenerateConnectingPoint::init()
		{
			addWorkingEntity(CP_TABLE);
			addWorkingEntity(EDGE_TABLE_INIT);
		}

		///
		///
		///
		void GenerateConnectingPoint::onCompute( bool verbose = false )
		{
			//--
			std::string idName = _epgParams.getValue( ID ).toString();
			std::string geomName = _epgParams.getValue( GEOM ).toString();

			//--
			std::string countryCodeW = _themeParams.getValue(COUNTRY_CODE_W).toString();
			std::string cpRefTableName = _themeParams.getValue(CP_TABLE).toString();

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
			_epgParams.setParameter(EDGE_TABLE, ign::data::String(getCurrentWorkingTableName(EDGE_TABLE_INIT)));
			ome2::utils::CopyTableUtils::copyEdgeTable(getLastWorkingTableName(EDGE_TABLE_INIT), "", false, true);
			
			//--
			_themeParams.setParameter(CP_TABLE, ign::data::String(getCurrentWorkingTableName(CP_TABLE)));

			//--
			app::calcul::CFeatGenerationOp::ComputeCP(countryCodeW, verbose);

			//--
			_themeParams.setParameter(CP_TABLE, ign::data::String(cpRefTableName));
		}

	}
}