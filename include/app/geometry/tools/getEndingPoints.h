#ifndef _APP_GEOMETRY_TOOLS_GETENDINGPOINTS_H_
#define _APP_GEOMETRY_TOOLS_GETENDINGPOINTS_H_

#include <ign/geometry/Geometry.h>


namespace app{
namespace geometry{
namespace tools{


	/// \brief Get the geometry of the main axis from a polygon
	std::pair<ign::geometry::Point, ign::geometry::Point>  getEndingPoints( 
		ign::geometry::LineString const& ls 
		);

}
}
}

#endif
