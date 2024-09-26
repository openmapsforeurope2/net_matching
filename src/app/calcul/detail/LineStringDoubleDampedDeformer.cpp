// #include <epg/calcul/matching/detail/LineStringAbsDampedDeformer.h>

#include <app/calcul/detail/LineStringAbsDampedDeformer.h>

using namespace epg::calcul::matching::detail;

///
///
///
LineStringAbsDampedDeformer::LineStringAbsDampedDeformer( 
	double absThreshold, 
	double influenceFactor,
	double snapDist
):
	_absThreshold( absThreshold ),
	_influenceFactor( influenceFactor ),
	_snapDist( snapDist )
{
}

///
///
///
LineStringAbsDampedDeformer::~LineStringAbsDampedDeformer()
{
}

///
///
///
void LineStringAbsDampedDeformer::deform( 
	ign::math::Vec2d const& startDeformation,
	ign::math::Vec2d const& endDeformation,
	ign::geometry::LineString& ls
	)const
{
	if(startDeformation.norm() == 0 && endDeformation.norm() == 0) return;
	
	bool hasZ = ls.is3D();
	double length = ls.length();

	ign::geometry::Point startPoint = ls.startPoint();
	ign::geometry::Point endPoint = ls.endPoint();

	double startAffectDist = std::min( std::min( startDeformation.norm() * _influenceFactor, length ), _absThreshold);
	double endAffectDist = std::min( std::min( endDeformation.norm() * _influenceFactor, length ), _absThreshold);

	if ( _snapDist >= 0 && startAffectDist != 0) {
		//ajouter un point de fin de déformation
		double abs = 0. ;
		double previousAbs = 0.;
		ign::geometry::Point previousPoint = startPoint;
		for( size_t i = 1 ; i < ls.numPoints()-1 ; ++i )
		{
			ign::geometry::Point& p = ls.pointN( i ) ;
			double lengthSegment = p.distance2d( previousPoint );
			abs += lengthSegment;
			
			if ( abs > startAffectDist ) {
				if ( (startAffectDist - previousAbs) < _snapDist ) break;
				if ( (abs - startAffectDist) < _snapDist ) break;

				double absSegment = abs - startAffectDist;
				double ratio = absSegment / lengthSegment;

				ign::geometry::Point newPoint = p;
				newPoint.x() += (previousPoint.x() - newPoint.x()) * ratio;
				newPoint.y() += (previousPoint.y() - newPoint.y()) * ratio;
				if ( hasZ ) newPoint.z() += (previousPoint.z() - newPoint.z()) * ratio;

				if( previousPoint.distance(newPoint) )

				ls.addPoint(newPoint, i);
				break;
			}
			previousAbs = abs;
			previousPoint = p;
		}
	}

	if ( _snapDist >= 0 && endAffectDist != 0) {
		//ajouter un point de fin de déformation
		double abs = 0. ;
		double previousAbs = 0.;
		ign::geometry::Point previousPoint = endPoint;
		for( int i = ls.numPoints()-2 ; i >= 0 ; --i )
		{
			ign::geometry::Point& p = ls.pointN( i ) ;
			double lengthSegment = p.distance2d( previousPoint );
			abs += lengthSegment;
			
			if ( abs > endAffectDist ) {
				if ( (endAffectDist - previousAbs) < _snapDist ) break;
				if ( (abs - endAffectDist) < _snapDist ) break;

				double absSegment = abs - endAffectDist;
				double ratio = absSegment / lengthSegment;

				ign::geometry::Point newPoint = p;
				newPoint.x() += (previousPoint.x() - newPoint.x()) * ratio;
				newPoint.y() += (previousPoint.y() - newPoint.y()) * ratio;
				if ( hasZ ) newPoint.z() += (previousPoint.z() - newPoint.z()) * ratio;

				ls.addPoint(newPoint, i+1);
				break;
			}
			previousAbs = abs;
			previousPoint = p;
		}
	}
	
	double abs = 0. ;
	ign::geometry::Point previousPoint = startPoint;
	for( size_t i = 1 ; i < ls.numPoints()-1 ; ++i )
	{
		ign::geometry::Point& p = ls.pointN( i ) ;
		abs += p.distance2d( previousPoint );
		previousPoint = p;

		double dampingFactorStart = abs >= startAffectDist ? 0 : (startAffectDist - abs) / startAffectDist ;
		double dampingFactorEnd = (length - abs) >= endAffectDist ? 0 : (endAffectDist - (length - abs)) / endAffectDist ;

		if (dampingFactorStart == 0 && dampingFactorEnd == 0) continue;

		ign::math::Vec2d translation = startDeformation * dampingFactorStart + endDeformation * dampingFactorEnd;
		p += translation;
	}

	ls.startPoint() += startDeformation;
	ls.endPoint() += endDeformation;
}