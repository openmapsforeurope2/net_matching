
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
		_initParameter( DB_CONF_FILE, "DB_CONF_FILE" );
		_initParameter( EDGE_TABLE_INIT, "EDGE_TABLE_INIT" );
		_initParameter( COUNTRY_CODE_W, "COUNTRY_CODE_W" );
		_initParameter( W_TAG, "W_TAG" );
		
		_initParameter(BOUNDARY_SMOOTHED_TABLE, "BOUNDARY_SMOOTHED_TABLE");
		_initParameter( LANDMASK_TABLE, "LANDMASK_TABLE" );
		_initParameter( LAND_COVER_TYPE, "LAND_COVER_TYPE" );
		_initParameter( TYPE_LAND_AREA, "TYPE_LAND_AREA" );
		_initParameter( EDGE_LINK, "EDGE_LINK" );
		_initParameter( FORM_OF_WAY, "FORM_OF_WAY" );
		_initParameter( SNAP_DIST, "SNAP_DIST");
		_initParameter( ANGLE_MAX_2_CUT_BORDER, "ANGLE_MAX_2_CUT_BORDER" );
		_initParameter( CP_TABLE, "CP_TABLE");
		_initParameter( CL_TABLE, "CL_TABLE");
		_initParameter( CP_TABLE_SUFFIX, "CP_TABLE_SUFFIX");
		_initParameter( CL_TABLE_SUFFIX, "CL_TABLE_SUFFIX");
		_initParameter( CF_STATUS, "CF_STATUS");

		_initParameter( CP_MERGE_DIST_CP, "CP_MERGE_DIST_CP");
		_initParameter(CP_MERGE_DIST_TRACTOR_CP, "CP_MERGE_DIST_TRACTOR_CP");
		_initParameter(CP_VALUE_FORMWAY_BIGDIST2MERGE, "CP_VALUE_FORMWAY_BIGDIST2MERGE");
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
		_initParameter( LIST_ATTR_JSON, "LIST_ATTR_JSON");
		_initParameter( LIST_ATTR_W, "LIST_ATTR_W");
		_initParameter( SQL_FILTER_EDGES_2_GENERATE_CF, "SQL_FILTER_EDGES_2_GENERATE_CF");

		_initParameter( EC_LANDMASK_BUFFER, "EC_LANDMASK_BUFFER");
		_initParameter( EC_SNAP_DIST, "EC_SNAP_DIST");
		_initParameter( EC_SNAP_2_EDGE_END_DIST, "EC_SNAP_2_EDGE_END_DIST");
		
		_initParameter( ECL_SQL_FILTER, "ECL_SQL_FILTER");
		_initParameter( ECL_SLIM_SURFACE_WIDTH, "ECL_SLIM_SURFACE_WIDTH");
		_initParameter( ECL_SLIM_SURFACE_MAX_AREA, "ECL_SLIM_SURFACE_MAX_AREA" );
		_initParameter( ECL_SLIM_SURFACE_MAX_NB_POINTS, "ECL_SLIM_SURFACE_MAX_NB_POINTS" );
		_initParameter( ECL_ANTENNA_RATIO_THRESHOLD, "ECL_ANTENNA_RATIO_THRESHOLD");
		_initParameter( ECL_ANTENNA_RATIO_WITH_BUFFER_THRESHOLD, "ECL_ANTENNA_RATIO_WITH_BUFFER_THRESHOLD");
		_initParameter( ECL_LANDMASK_BUFFER, "ECL_LANDMASK_BUFFER");
		_initParameter( ECL_ANTENNA_MIN_LENGTH, "ECL_ANTENNA_MIN_LENGTH");
		_initParameter( ECL_ANTENNA_MIN_LENGTH_IN_COUNTRY, "ECL_ANTENNA_MIN_LENGTH_IN_COUNTRY");
		_initParameter( ECL_PARALELLE_EDGE_MAX_DIST, "ECL_PARALELLE_EDGE_MAX_DIST");
		_initParameter( ECL_ANTENNA_MIN_DIST_2_NEIGHBOR, "ECL_ANTENNA_MIN_DIST_2_NEIGHBOR");
		_initParameter(ECL_TINY_EDGE_MAX_LENGTH, "ECL_TINY_EDGE_MAX_LENGTH");

		_initParameter(DIST_MAX_JUNCTIONS, "DIST_MAX_JUNCTIONS");
		_initParameter(ANGLE_MAX_ORIENTATION_EDGES, "ANGLE_MAX_ORIENTATION_EDGES");
	
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