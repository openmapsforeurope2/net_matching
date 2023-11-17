
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
		_initParameter( FORM_OF_WAY, "FORM_OF_WAY" );
		_initParameter( SNAP_DIST, "SNAP_DIST" );
		_initParameter( CP_TABLE, "CP_TABLE");
		_initParameter( CL_TABLE, "CL_TABLE");
		_initParameter( CF_STATUS, "CF_STATUS");

		_initParameter( CP_MERGE_DIST_CP, "CP_MERGE_DIST_CP");
		_initParameter( CP_MERGE_DIST_TRACTOR_CP, "CP_MERGE_DIST_TRACTOR_CP");
		_initParameter( CP_BUFFER_DIST, "CP_BUFFER_DIST");
		_initParameter( CP_INTERSECTED_CL_DIST, "CP_INTERSECTED_CL_DIST");
		_initParameter( CP_UNDERSHOOT_DIST, "CP_UNDERSHOOT_DIST");
		_initParameter( CP_CP_2_CL_SNAP_DIST, "CP_CP_2_CL_SNAP_DIST");
		_initParameter( CP_VERTEX_CL_SNAP_DIST, "CP_VERTEX_CL_SNAP_DIST");

		_initParameter( CL_BUFFER_DIST, "CL_BUFFER_DIST");
		_initParameter( CL_THRESHOLD_NO_CL, "CL_THRESHOLD_NO_CL");
		_initParameter( CL_RATIO_IN_BUFFER, "CL_RATIO_IN_BUFFER");
		_initParameter( CL_SNAP_ON_VERTEX_BORDER_DIST, "CL_SNAP_ON_VERTEX_BORDER_DIST");
		_initParameter( CL_CL_CLOSEST_MAX_DIST, "CL_CL_CLOSEST_MAX_DIST");
		_initParameter( CL_BORDER_MAX_ANGLE, "CL_BORDER_MAX_ANGLE");
		_initParameter( CL_EDGE_MAX_ANGLE, "CL_EDGE_MAX_ANGLE");
		_initParameter( CL_CL_INTERSECTED_DIST, "CL_CL_INTERSECTED_DIST");
		_initParameter( CL_MERGE_CL_DIST, "CL_MERGE_CL_DIST");
		_initParameter( CL_EDGE_MAX_DIST, "CL_EDGE_MAX_DIST");
		_initParameter( CL_SNAP_PROJ_CL_2_EDGE_DIST, "CL_SNAP_PROJ_CL_2_EDGE_DIST");
		_initParameter( CL_CL_2_MERGE_MIN_LENGTH, "CL_CL_2_MERGE_MIN_LENGTH");

		_initParameter( LIST_ATTR_TO_CONCAT, "LIST_ATTR_TO_CONCAT");
		_initParameter( LIST_ATTR_W, "LIST_ATTR_W");
		_initParameter( SQL_FILTER_EDGES_2_GENERATE_CF, "SQL_FILTER_EDGES_2_GENERATE_CF");
		_initParameter( SLIM_SURFACE_WIDTH, "SLIM_SURFACE_WIDTH");
		_initParameter( ANTENNA_RATIO_THRESHOLD, "ANTENNA_RATIO_THRESHOLD");
		_initParameter( ANTENNA_RATIO_WITH_BUFFER_THRESHOLD, "ANTENNA_RATIO_WITH_BUFFER_THRESHOLD");
		_initParameter( LANDMASK_BUFFER, "LANDMASK_BUFFER");
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