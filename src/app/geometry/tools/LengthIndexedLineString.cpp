#include <app/geometry/tools/LengthIndexedLineString.h>
#include <ign/math/Line2T.h>

using namespace app;
using namespace app::geometry;
using namespace app::geometry::tools;

///
///
///
LengthIndexedLineString::LengthIndexedLineString():
	_lineString(),
	_index(),
    _is3d(false)
{

}

///
///
///
LengthIndexedLineString::LengthIndexedLineString( ign::geometry::LineString const& lineString ):
	_lineString(lineString),
	_index(),
    _is3d(lineString.is3D())
{
	_computeIndex();
}

///
///
///
LengthIndexedLineString::LengthIndexedLineString( LengthIndexedLineString const& other ):
	_lineString(other._lineString),
	_index(other._index),
    _is3d(other._lineString.is3D())
{

}

///
///
///
LengthIndexedLineString& LengthIndexedLineString::operator = ( LengthIndexedLineString const& other )
{
	_lineString = other._lineString;
	_index      = other._index;
    _is3d       = other._is3d;
	return *this;
}

///
///
///
LengthIndexedLineString::~LengthIndexedLineString()
{

}

///
///
///
ign::geometry::LineString const& LengthIndexedLineString::getLineString() const
{
	return _lineString;
}

///
///
///
void LengthIndexedLineString::setLineString( ign::geometry::LineString const& lineString )
{
	_lineString = lineString;
	_computeIndex();
}

///
///
///
double LengthIndexedLineString::getPointLocation( size_t numPoint ) const
{
	IGN_ASSERT( numPoint < _lineString.numPoints() );
	return _index[numPoint];
}

///
///
///
ign::geometry::LineString LengthIndexedLineString::getSubLineString( double const& sBegin, double const& sEnd )
{
	std::vector< double > abscisses;

	//point de depart
	abscisses.push_back( sBegin );

	//parcours de la fin au debut
	if ( sBegin <= sEnd ){
		for ( std::vector< double >::const_iterator it = _index.begin(); it != _index.end(); it++ ){
			if ( *it > sBegin && *it < sEnd ){
				abscisses.push_back( *it );
			}
		}
	}else{
		for ( std::vector< double >::const_reverse_iterator it = _index.rbegin(); it != _index.rend(); it++ ){
			if ( *it > sEnd && *it < sBegin ){
				abscisses.push_back( *it );
			}
		}
	}

	//point de fin
	abscisses.push_back( sEnd );

	ign::geometry::LineString res;
	res.reserve( abscisses.size() );
	for ( std::vector< double >::const_iterator it = abscisses.begin(); it != abscisses.end(); ++it )
    {
		res.addPoint( locateAlong( *it ) );
	}
	return res;
}

///
///
///
std::vector< double > const & LengthIndexedLineString::getPointAbscisses() const
{
	return _index;
}


///
///
///
double LengthIndexedLineString::length() const
{
	if ( _lineString.isEmpty() ){
		return ign::numeric::Numeric< double >::NaN();
	}
	return _index.back();
}

///
///
///
size_t LengthIndexedLineString::numLines() const
{
	size_t npt = _lineString.numPoints();
	if ( npt < 2 ){
		return 0;
	}else{
		return npt - 1 ;
	}
}

#include <iostream>

///
///
///
ign::geometry::Point LengthIndexedLineString::locateAlong( double const& s ) const
{
	IGN_ASSERT( _lineString.isValid() && ! _lineString.isEmpty() );
	size_t numSeg = _findLine( s );
	double ds = s - getPointLocation( numSeg );
	ign::math::Vector2T< double > vA = _lineString[ numSeg ].toVec2d();
	ign::math::Vector2T< double > vB = _lineString[ numSeg + 1 ].toVec2d();
	double zA = _is3d ? _lineString[numSeg].z():0;
	double zB = _is3d ? _lineString[numSeg+1].z():0;

    
	ign::math::Vector2T< double > AB = (vB - vA);
	double n = AB.norm();
	if ( ign::numeric::Numeric< double >::IsNULL(n) ){
		ign::geometry::Point res( vA[0], vA[1] );
		res.m() = s;
		return res;
	}
	ign::geometry::Point res(
		vA[0] + AB[0]*(ds/n),
		vA[1] + AB[1]*(ds/n)
	);
    if (_is3d) res.z() = zA + (zB-zA)*(ds/n);
	res.m() = s;
	return res;
}

///
///
///
double LengthIndexedLineString::project( ign::geometry::Point const& p )const
{
	double minDistance = std::numeric_limits<double>::max();

	ign::geometry::LineString::const_iterator it1 = _lineString.begin();
	ign::geometry::LineString::const_iterator it2 = it1;

	double distance;
	double AbsoluteAbs(0.);
	for( ++it2 ; it2 != _lineString.end() ; ++it1, ++it2 )
	{
        ign::math::Vec2d pt1 = it1->toVec2d();
        ign::math::Vec2d pt2 = it2->toVec2d();
        
        ign::math::Line2d line( pt1, pt2);
        ign::math::Vec2d v = p.toVec2d();
            
        double RelativeAbs = line.project( v, true );
        ign::math::Vec2d vProj = line.interpolate( RelativeAbs ) ;
        distance = v.distance2( vProj );
        
        if( distance < minDistance )
        {
            AbsoluteAbs = _index[ it1 - _lineString.begin() ] + line.a().distance( vProj );
            minDistance = distance;
        }
	}

	return AbsoluteAbs;
}

///
///
///
void LengthIndexedLineString::_computeIndex()
{
	if( ! _lineString.isValid() || _lineString.isEmpty() ){
		_index.resize(0);
		return;
	}

	_index.resize( _lineString.numPoints() );
	double L = 0.0;
	for ( size_t i = 0; i < _lineString.numPoints() - 1 ; i++ ){
		ign::geometry::Point const & a = _lineString.pointN(i);
		ign::geometry::Point const & b = _lineString.pointN(i+1);
		_index[i] = L;
#ifdef IGN_WITH_GEOS
		L += a.distance( b );
#else
		double dx = a.x() - b.x();
		double dy = a.y() - b.y();
        
		L += std::sqrt(dx*dx+dy*dy);
#endif
	}
	_index.back() = L;
}

///
///
///
size_t LengthIndexedLineString::_findLine( double const& s ) const
{

	if ( s <= 0.0 ){
		return 0;
	}
	for ( size_t i = 1; i < numLines(); i++ ){
		if ( _index[i] >= s ){
			return i-1;
		}
	}
	return numLines()-1;
}

