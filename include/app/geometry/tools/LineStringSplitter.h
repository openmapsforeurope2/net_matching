#ifndef _APP_TOOLS_GEOMETRY_LINESTRINGSPLITTER_H_
#define _APP_TOOLS_GEOMETRY_LINESTRINGSPLITTER_H_

#include <set>

#include <ign/geometry/LineString.h>
#include <ign/geometry/algorithm/LineIntersectorOpGeos.h>
#include <ign/geometry/index/QuadTree.h>

namespace app{
namespace geometry{
namespace tools{

	class LineStringSplitter{

	public:
		/// \brief 
		LineStringSplitter( ign::geometry::LineString const& ls, double precision = 1e-5 );

		/// \brief
		~LineStringSplitter();

		/// \brief
		void addCuttingGeometry( ign::geometry::Geometry const& geom );

		/// \brief
		std::vector< ign::geometry::LineString > trimStart()const;

		/// \brief
		std::vector< ign::geometry::LineString > trimEnd()const;

		/// \brief
		ign::geometry::LineString truncAtEnds()const;

		/// \brief
		std::vector< ign::geometry::LineString > getSubLineStrings()const;

		/// \brief
		std::vector< ign::geometry::LineString > getSubLineStringsZ()const;

	private:
		//--
		void _addCuttingGeometry( ign::geometry::LineString const& ls );
		
		//--
		void _addCuttingGeometry( ign::geometry::Point const& pt );

		//--
		void _addCuttingGeometry( ign::geometry::Polygon const& poly );

		//--
		void _addCuttingGeometry( ign::geometry::MultiLineString const& mls );
		
		//--
		void _addCuttingGeometry( ign::geometry::MultiPoint const& mpt );

		//--
		void _addCuttingGeometry( ign::geometry::MultiPolygon const& mp );

		//--
		std::vector< ign::geometry::LineString > _cut( std::pair< int, double > const& cutAbs )const;


	private:

		//--
		ign::geometry::LineString const&                    _lsRef;
	
		//--
		std::vector< std::set< double > >                   _vCuttings;

		//--
		ign::geometry::algorithm::LineIntersectorOpGeos     _intersector;

		//--
		double                                              _precision;

		//--
		ign::geometry::index::QuadTree< int >              _qTreeSegment;

	
	};

}
}
}

#endif