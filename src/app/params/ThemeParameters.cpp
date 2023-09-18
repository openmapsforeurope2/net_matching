
//APP
#include <app/params/ThemeParameters.h>

//SOCLE
#include <ign/Exception.h>


namespace app{
namespace params{


	///
	///
	///
	ThemeParameters::ThemeParameters()
	{
		_initParameter( LANDMASK_TABLE, "LANDMASK_TABLE" );
		_initParameter( LAND_COVER_TYPE, "LAND_COVER_TYPE" );
		_initParameter( TYPE_LAND_AREA, "TYPE_LAND_AREA" );
		_initParameter( EDGE_LINK, "EDGE_LINK" );
		_initParameter( SNAP_DIST, "SNAP_DIST" );
		_initParameter( CP_TABLE, "CP_TABLE");
		_initParameter( CL_TABLE, "CL_TABLE");
		_initParameter( CF_STATUS, "CF_STATUS");
		_initParameter( LIST_ATTR_TO_CONCAT, "LIST_ATTR_TO_CONCAT");
		_initParameter( LIST_ATTR_W, "LIST_ATTR_W");
		_initParameter( SQL_FILTER_EDGES_2_GENERATE_CF, "SQL_FILTER_EDGES_2_GENERATE_CF");
		_initParameter( SLIM_SURFACE_WIDTH, "SLIM_SURFACE_WIDTH");
		_initParameter( ANTENNA_RATIO_THRESHOLD, "ANTENNA_RATIO_THRESHOLD");
	}

	///
	///
	///
	ThemeParameters::~ThemeParameters()
	{
	}

	///
	///
	///
	std::string ThemeParameters::getClassName()const
	{
		return "app::params::ThemeParameters";
	}



}
}