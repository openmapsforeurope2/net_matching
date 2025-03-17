#ifndef _APP_TOOLS_GEOMETRY_LINESTRINGSPLITTER_H_
#define _APP_TOOLS_GEOMETRY_LINESTRINGSPLITTER_H_

#include <set>

#include <ign/geometry/LineString.h>
#include <ign/geometry/algorithm/LineIntersectorOpGeos.h>
#include <ign/geometry/index/QuadTree.h>
#include <epg/log/EpgLogger.h>

namespace app{
namespace geometry{
namespace tools{

	/// @brief Classe servant à découper les polylignes
	class LineStringSplitter{

	public:
		/// @brief Constructeur
		/// @param ls Polyligne à découper
		/// @param precision Précision
		LineStringSplitter( ign::geometry::LineString const& ls, double precision = 1e-5 );

		/// @brief Destructeur
		~LineStringSplitter();

		/// @brief Ajoute une géométrie de découpe. Si c'est une géométrie surfacique 
		/// celle-ci est considérée comme une collection de polyligne
		/// @param geom Géométrie de découpe
		void addCuttingGeometry( ign::geometry::Geometry const& geom );

		/// @brief Calcule la découpe en supprimant l'extrémité de départ
		/// @return Polylignes résultant de la découpe
		std::vector< ign::geometry::LineString > trimStart()const;

		/// @brief Calcule la découpe en supprimant l'extrémité de fin
		/// @return Polylignes résultant de la découpe
		std::vector< ign::geometry::LineString > trimEnd()const;

		/// @brief Calcule la découpe en supprimant les extremités
		/// @return Polylignes résultant de la découpe 
		ign::geometry::LineString truncAtEnds()const;

		/// @brief Calcule la découpe.
		/// @return Polylignes résultant de la découpe
		std::vector< ign::geometry::LineString > getSubLineStrings()const;

		/// @brief Calcule la découpe en conservant le Z au niveau des points de découpe
		/// @return Polylignes résultant de la découpe
		std::vector< ign::geometry::LineString > getSubLineStringsZ()const;

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
		//--
		epg::log::EpgLogger*                               _logger;

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
	};

}
}
}

#endif