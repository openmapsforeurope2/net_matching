#include <app/geometry/tools/getEndingPoints.h>

//SOCLE
#include <ign/geometry/graph/tools/SnapRoundPlanarizer.h>
#include <ign/geometry/algorithm/TriangulationDelaunayCGAL.h>
#include <ign/geometry/algorithm/SimplifyOpGeos.h>


//EPG
#include <epg/calcul/selection/detail/MedialAxisGraph.h>
#include <epg/graph/algorithm/getLonguestPath.h>
#include <epg/graph/tools/convertPathToLineString.h>
#include <epg/tools/geometry/simplifyLineString.h>
#include <epg/graph/tools/reversePath.h>
#include <epg/detail/SimpleGraph.h>
#include <epg/log/ShapeLogger.h>


/////
/////

namespace app{
namespace geometry{
namespace tools{

	///
	///
	///
	std::pair<ign::geometry::Point, ign::geometry::Point> getEndingPoints( 
		ign::geometry::LineString const& ls
		)
	{
		ign::geometry::LineString mainAxis;

		ign::geometry::algorithm::TriangulationDelaunayCGAL triOp;
		ign::geometry::MultiLineString mls;
		triOp.medialAxis( ls, mls );

		
		typedef epg::calcul::selection::detail::MedialAxisGraphType             MedialAxisGraphType;
		typedef MedialAxisGraphType::edges_path                                 edges_path;
		typedef MedialAxisGraphType::edge_iterator                              edge_iterator;

		MedialAxisGraphType                                                     medialAxisGraph;
		ign::geometry::graph::tools::SnapRoundPlanarizer<MedialAxisGraphType>   builder( medialAxisGraph, 1e7 );


		for( size_t i = 0 ; i < mls.numGeometries() ; ++i ) builder.addEdge( mls.lineStringN(i) );
		
		builder.planarize();


		edge_iterator eit, eend;
		for( medialAxisGraph.edges( eit, eend ) ; eit != eend ; ++eit )
			medialAxisGraph[*eit].weight = medialAxisGraph.getGeometry(*eit).length();


		edges_path path;
		double dist = epg::graph::algorithm::getLonguestPath( medialAxisGraph, path );
        if ( dist <= 0 ) {
            return std::make_pair(ign::geometry::Point(), ign::geometry::Point());
        }

		// DEBUG
		ign::geometry::LineString pathLs = epg::graph::tools::convertPathToLineString( medialAxisGraph, path );
		ign::feature::Feature feat;
		feat.setGeometry(pathLs);
		epg::log::ShapeLoggerS::getInstance()->writeFeature("ecl_slim_surface_medial_axis", feat);
		
		
		ign::geometry::Point pStart = medialAxisGraph.getGeometry( medialAxisGraph.source( *path.begin() ) );
		ign::geometry::Point pTarget = medialAxisGraph.getGeometry( medialAxisGraph.target( *path.rbegin() ) );
		
		return std::make_pair(pStart, pTarget);
	}

}
}
}