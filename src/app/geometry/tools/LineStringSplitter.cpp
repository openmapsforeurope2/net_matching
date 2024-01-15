#include <app/geometry/tools/LineStringSplitter.h>
#include <app/geometry/tools/LengthIndexedLineString.h>

//SOCLE
#include <ign/geometry/Polygon.h>
#include <ign/geometry/MultiPoint.h>
#include <ign/data/all.h>

//EPG
#include <epg/tools/geometry/interpolate.h>
#include <epg/tools/geometry/project.h>



using namespace app;
using namespace app::geometry::tools ;


LineStringSplitter::LineStringSplitter( ign::geometry::LineString const& ls, double precision ):_lsRef( ls ),_precision(precision)
{
	_vCuttings.resize( _lsRef.numSegments() );

	_qTreeSegment.ensureExtent( _lsRef.getEnvelope() );
	for( size_t i = 0 ; i < _lsRef.numSegments() ; ++i )
	{
//		if( _lsRef.pointN(i).equals( _lsRef.pointN(i+1) ) )
//			IGN_THROW_EXCEPTION( "[epg::tools::geometry::LineStringSplitter::LineStringSplitter] null segment in Linestring '"+_lsRef.toString()+"' between points [index] "+ign::data::Integer(i).toString()+" and "+ign::data::Integer(i+1).toString()+"." );
		
		_qTreeSegment.insert( i, ign::geometry::Envelope( _lsRef.pointN(i), _lsRef.pointN(i+1) ) );
	}
}


LineStringSplitter::~LineStringSplitter()
{
}

///
///
///
void LineStringSplitter::addCuttingGeometry( ign::geometry::Geometry const& geom )
{
	switch( geom.getGeometryType() )
	{
		case ign::geometry::Geometry::GeometryTypePoint :{
			_addCuttingGeometry( geom.asPoint() );
			break;}
		case ign::geometry::Geometry::GeometryTypeLineString :{
			_addCuttingGeometry( geom.asLineString() );
			break;}
		case ign::geometry::Geometry::GeometryTypePolygon :{
			_addCuttingGeometry( geom.asPolygon() );
			break;}
		case ign::geometry::Geometry::GeometryTypeMultiPoint :{
			_addCuttingGeometry( geom.asMultiPoint() );
			break;}
		case ign::geometry::Geometry::GeometryTypeMultiLineString :{
			_addCuttingGeometry( geom.asMultiLineString() );
			break;}
		case ign::geometry::Geometry::GeometryTypeMultiPolygon :{
			_addCuttingGeometry( geom.asMultiPolygon() );
			break;}
		default :{
			IGN_THROW_EXCEPTION( "epg::tools::geometry::LineStringSplitter:addCuttingGeometry : geometry type "+geom.getGeometryTypeName()+" not allowed" );
			break;
		}
	}
}

///
///
///
std::vector< ign::geometry::LineString > LineStringSplitter::trimStart()const
{
	std::pair< int, double > firstCut(0,0.);
	for( size_t i = 0 ; i < _lsRef.numSegments() ; ++i )
	{
		if( _vCuttings[i].empty() ) continue;
		firstCut = std::make_pair( i, *_vCuttings[i].begin() );
		break;
	}

	return _cut( firstCut );
}

///
///
///
std::vector< ign::geometry::LineString > LineStringSplitter::trimEnd()const
{
	std::pair< int, double > lastCut(_lsRef.numSegments()-1,1.);
	for( size_t i = _lsRef.numSegments()-1 ; i > 0 ; --i )
	{
		if( _vCuttings[i].empty() ) continue;
		lastCut = std::make_pair( i, *_vCuttings[i].rbegin() );
		break;
	}

	return _cut( lastCut );
}

///
///
///
ign::geometry::LineString LineStringSplitter::truncAtEnds()const
{
	ign::geometry::LineString lsResult;

	LengthIndexedLineString indexedLs( _lsRef );
	double length = _lsRef.length();

	std::pair< int, double > firstCut(0,0.);
	for( size_t i = 0 ; i < _lsRef.numSegments() ; ++i )
	{
		if( _vCuttings[i].empty() ) continue;
		ign::geometry::Point firstCutPoint = epg::tools::geometry::interpolateZ( _lsRef, i, *_vCuttings[i].begin() );
		double test = indexedLs.project( firstCutPoint );
		if( indexedLs.project( firstCutPoint ) > .5*length ) break;
		firstCut = std::make_pair( i, *_vCuttings[i].begin() );
		break;
	}

	std::pair< int, double > lastCut(_lsRef.numSegments()-1,1.);
	for( size_t i = _lsRef.numSegments()-1 ; i > 0 ; --i )
	{
		if( _vCuttings[i].empty() ) continue;
		ign::geometry::Point lastCutPoint = epg::tools::geometry::interpolateZ( _lsRef, i, *_vCuttings[i].rbegin() );
		if( indexedLs.project( lastCutPoint ) < .5*length ) break;
		lastCut = std::make_pair( i, *_vCuttings[i].rbegin() );
		break;
	}

	
	for( int i = firstCut.first ; i < lastCut.first+1 ; ++i )
	{
		if( i == firstCut.first )
		{
			ign::geometry::Point firstCutPoint = epg::tools::geometry::interpolateZ( _lsRef, firstCut.first, firstCut.second );
			if( firstCut.second > 0.5 && firstCutPoint.distance2d( _lsRef.pointN( i+1 ) ) < _precision ) continue;

			lsResult.addPoint( firstCutPoint );
			continue;
		}

		lsResult.addPoint( _lsRef.pointN( i ) );
	}

	ign::geometry::Point lastCutPoint = epg::tools::geometry::interpolateZ( _lsRef, lastCut.first, lastCut.second );
	if( lastCutPoint.distance2d( lsResult.endPoint() ) > _precision )
	{
		lsResult.addPoint( lastCutPoint );
	}
	
	return lsResult;
}

///
///
///
std::vector< ign::geometry::LineString > LineStringSplitter::getSubLineStrings()const
{
	std::vector< ign::geometry::LineString > vSubLs;
	vSubLs.push_back( ign::geometry::LineString() );

	for( size_t i = 0 ; i < _lsRef.numSegments() ; ++i )
	{
		vSubLs.back().addPoint( _lsRef.pointN( i ) );

		std::set< double >::const_iterator cit;
		for( cit = _vCuttings[i].begin() ; cit != _vCuttings[i].end() ; ++cit )
		{
			double abs = *cit;

			//IGN_ASSERT( abs >= 0. && abs <= 1. );

			ign::geometry::Point pt = epg::tools::geometry::interpolate( _lsRef, i, abs );

			if( abs > 0.5 && pt.distance2d( _lsRef.pointN( i+1 ) ) < _precision )
			{
				vSubLs.back().addPoint( _lsRef.pointN( i+1 ) );
				vSubLs.push_back( ign::geometry::LineString() );
				break;
			}

			if( pt.distance2d( vSubLs.back().endPoint() ) < _precision )
			{
				pt = vSubLs.back().endPoint();
			}else{
				vSubLs.back().addPoint( pt );
			}

			if( vSubLs.back().numPoints() > 1 ) 
			{
				vSubLs.push_back( ign::geometry::LineString() );
				vSubLs.back().addPoint( pt );
			}
		}
	}

	if( vSubLs.back().isEmpty() ) vSubLs.pop_back();
	else vSubLs.back().addPoint( _lsRef.endPoint() );

	return vSubLs;
}

///
///
///
std::vector< ign::geometry::LineString > LineStringSplitter::getSubLineStringsZ()const
{
	std::vector< ign::geometry::LineString > vSubLs;
	vSubLs.push_back( ign::geometry::LineString() );

	for( size_t i = 0 ; i < _lsRef.numSegments() ; ++i )
	{
		vSubLs.back().addPoint( _lsRef.pointN( i ) );

		std::set< double >::const_iterator cit;
		for( cit = _vCuttings[i].begin() ; cit != _vCuttings[i].end() ; ++cit )
		{
			double abs = *cit;

			//IGN_ASSERT( abs >= 0. && abs <= 1. );

			ign::geometry::Point pt = epg::tools::geometry::interpolateZ( _lsRef, i, abs );

			if( abs > 0.5 && pt.distance2d( _lsRef.pointN( i+1 ) ) < _precision )
			{
				vSubLs.back().addPoint( _lsRef.pointN( i+1 ) );
				vSubLs.push_back( ign::geometry::LineString() );
				break;
			}

			if( pt.distance2d( vSubLs.back().endPoint() ) < _precision )
			{
				pt = vSubLs.back().endPoint();
			}else{
				vSubLs.back().addPoint( pt );
			}

			if( vSubLs.back().numPoints() > 1 ) 
			{
				vSubLs.push_back( ign::geometry::LineString() );
				vSubLs.back().addPoint( pt );
			}
		}
	}

	if( vSubLs.back().isEmpty() ) vSubLs.pop_back();
	else vSubLs.back().addPoint( _lsRef.endPoint() );

	return vSubLs;
}

///
///
///
void LineStringSplitter::_addCuttingGeometry( ign::geometry::LineString const& ls )
{
	for( size_t i = 0 ; i < ls.numSegments() ; ++i )
	{
		if( ls.pointN(i).equals( ls.pointN(i+1) ) ) {
			std::string mess = "[epg::tools::geometry::LineStringSplitter::_addCuttingGeometry] null segment in Linestring '"+ls.toString()+"' between points [index] "+ign::data::Integer(i).toString()+" and "+ign::data::Integer(i+1).toString()+".";
			_logger->log(epg::log::ERROR, mess);
			continue;
		}
		
		std::set< int > sSegments;
		_qTreeSegment.query( ign::geometry::Envelope( ls.pointN(i), ls.pointN(i+1) ), sSegments );

		std::set< int >::const_iterator sit;
		for( sit = sSegments.begin() ; sit != sSegments.end() ; ++sit )
		{
			_intersector.computeIntersection( _lsRef.pointN( *sit ), _lsRef.pointN( *sit+1 ), ls.pointN( i ), ls.pointN( i+1 ) );

			if( _intersector.getIntersectionNum() == 0 ) continue;

			for( size_t j = 0 ; j < _intersector.getIntersectionNum() ; ++j )
			{
				double abs = _intersector.abscisse( 0, j );

				_vCuttings[ *sit ].insert( abs );
			}
		}
	}
}

///
///
///
void LineStringSplitter::_addCuttingGeometry( ign::geometry::Point const& pt )
{
	std::pair< int, double > abs = epg::tools::geometry::projectAlong( _lsRef, pt, _precision/*1e-1*/ );

	_vCuttings[ abs.first ].insert( abs.second );
}

///
///
///
void LineStringSplitter::_addCuttingGeometry( ign::geometry::Polygon const& poly )
{
	for( size_t i = 0 ; i < poly.numRings() ; ++i )
		_addCuttingGeometry( poly.ringN( i ) );
}

///
///
///
void LineStringSplitter::_addCuttingGeometry( ign::geometry::MultiPoint const& mpt )
{
	for( size_t i = 0 ; i < mpt.numGeometries() ; ++i )
		_addCuttingGeometry( mpt.pointN( i ) );

}

///
///
///
void LineStringSplitter::_addCuttingGeometry( ign::geometry::MultiLineString const& mls )
{
	for( size_t i = 0 ; i < mls.numGeometries() ; ++i )
		_addCuttingGeometry( mls.lineStringN( i ) );
}

///
///
///
void LineStringSplitter::_addCuttingGeometry( ign::geometry::MultiPolygon const& mp )
{
	for( size_t i = 0 ; i < mp.numGeometries() ; ++i )
		_addCuttingGeometry( mp.polygonN( i ) );
}

///
///
///
std::vector< ign::geometry::LineString > LineStringSplitter::_cut( std::pair< int, double > const& cutAbs )const
{
	std::vector< ign::geometry::LineString > vResult( 1 );

	bool isFirstCut = true;
	for( size_t i = 0 ; i < _lsRef.numSegments() ; ++i )
	{
		vResult.back().addPoint( _lsRef.pointN( i ) );

		if( i != cutAbs.first ) continue;

		ign::geometry::Point firstCutPoint = epg::tools::geometry::interpolate( _lsRef, cutAbs.first, cutAbs.second );

		if( cutAbs.second < 0.5 && firstCutPoint.distance2d( _lsRef.pointN( i ) ) < _precision )
		{
			if( i != 0 ) 
			{
				vResult.push_back( ign::geometry::LineString() );
				vResult.back().addPoint( _lsRef.pointN( i ) );
			}
		}
		else if( cutAbs.second > 0.5 && firstCutPoint.distance2d( _lsRef.pointN( i+1 ) ) < _precision )
		{
			if( i != _lsRef.numSegments()-1 )
			{
				vResult.back().addPoint( _lsRef.pointN( i+1 ) );
				vResult.push_back( ign::geometry::LineString() );
			}
		}else{
			vResult.back().addPoint( firstCutPoint );
			vResult.push_back( ign::geometry::LineString() );
			vResult.back().addPoint( firstCutPoint );
		}
	}

	vResult.back().addPoint( _lsRef.endPoint() );

	IGN_ASSERT( vResult.size() <3 );
	
	return vResult;
}