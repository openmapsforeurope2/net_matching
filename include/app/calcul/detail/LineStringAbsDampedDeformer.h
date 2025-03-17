#ifndef _EPG_CALCUL_MATCHING_DETAIL_LINESTRINGSABSDAMPEDDEFORMER_H_
#define _EPG_CALCUL_MATCHING_DETAIL_LINESTRINGSABSDAMPEDDEFORMER_H_

// EPG
#include <epg/calcul/matching/detail/LineStringDeformer.h>


namespace epg{
namespace calcul{
namespace matching{
namespace detail{

	/// \brief Classe utilitaire pour la déformation amortie de linéaire
	class LineStringAbsDampedDeformer : public LineStringDeformer {
	public:
		/// @brief Constructeur
		/// @param absThreshold Abscisse curviligne définissant la limite de déformation
		/// @param influenceFactor Facteur permettant de définir la zone d'inluence d'un vecteur de déformation (cercle de rayon = norme du vecteur x influenceFactor)
		/// @param snapDist Distance de snapping permettant de prendre un point existant comme limite de déformation
		LineStringAbsDampedDeformer(double absThreshold, double influenceFactor, double snapDist);

		/// \brief Destructeur
		virtual ~LineStringAbsDampedDeformer();

		/// @brief Déforme la polyligne donnée en argument.
		/// @param startDeformation Vecteur de déformation appliqué au point de départ
		/// @param endDeformation Vecteur de déformation appliqué au dernier point
		/// @param ls Polyligne à déformer
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