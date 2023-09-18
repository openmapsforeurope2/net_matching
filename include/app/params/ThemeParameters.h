#ifndef _APP_PARAMS_THEMEPARAMETERS_H_
#define _APP_PARAMS_THEMEPARAMETERS_H_

//STL
#include <string>

//EPG
#include <epg/params/ParametersT.h>
#include <epg/SingletonT.h>



	enum TN_PARAMETERS{
		LANDMASK_TABLE,
		LAND_COVER_TYPE,
		TYPE_LAND_AREA,
		EDGE_LINK,
		SNAP_DIST,
		CP_TABLE,
		CL_TABLE,
		CF_STATUS,
		LIST_ATTR_TO_CONCAT,
		LIST_ATTR_W,
		SQL_FILTER_EDGES_2_GENERATE_CF,
		SLIM_SURFACE_WIDTH,
		ANTENNA_RATIO_THRESHOLD
	};

namespace app{
namespace params{

	class ThemeParameters : public epg::params::ParametersT< TN_PARAMETERS >
	{
		typedef  epg::params::ParametersT< TN_PARAMETERS > Base;

		public:

			/// \brief
			ThemeParameters();

			/// \brief
			~ThemeParameters();

			/// \brief
			virtual std::string getClassName()const;

	};

	typedef epg::Singleton< ThemeParameters >   ThemeParametersS;

}
}

#endif