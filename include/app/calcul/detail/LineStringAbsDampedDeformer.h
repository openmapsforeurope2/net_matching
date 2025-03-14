#ifndef _EPG_CALCUL_MATCHING_DETAIL_LINESTRINGSABSDAMPEDDEFORMER_H_
#define _EPG_CALCUL_MATCHING_DETAIL_LINESTRINGSABSDAMPEDDEFORMER_H_

// EPG
#include <epg/calcul/matching/detail/LineStringDeformer.h>


namespace epg{
namespace calcul{
namespace matching{
namespace detail{

	/// \brief
	class LineStringAbsDampedDeformer : public LineStringDeformer{
	public:
		/// @brief 
		/// @param absThreshold 
		/// @param influenceFactor 
		/// @param snapDist 
		LineStringAbsDampedDeformer(double absThreshold, double influenceFactor, double snapDist);

		/// \brief
		virtual ~LineStringAbsDampedDeformer();

		/// @brief 
		/// @param startDeformation 
		/// @param endDeformation 
		/// @param ls 
		virtual void deform( 
			ign::math::Vec2d const& startDeformation,
			ign::math::Vec2d const& endDeformation,
			ign::geometry::LineString& ls
			)const;

        private:
			//--
		    double _absThreshold;
			//--
            double _influenceFactor;
			//--
            double _snapDist;
	};
}
}
}
}

#endif