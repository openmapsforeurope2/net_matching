#ifndef _APP_CALCUL_DETAIL_RATIOTOOLS_H_
#define _APP_CALCUL_DETAIL_RATIOTOOLS_H_

// SOCLE
#include <ign/feature/FeatureStore.h>

// APP
#include <epg/Context.h>

namespace app{
namespace calcul{
namespace detail{

    /// @brief Claase d'outil pour le calcul de ratio (proportion de polyligne incluse dans des surfaces)
    class RatioTools {
        public:

        /// @brief fonction calculant la proportion de polyligne incluse dans des surfaces
        /// @param ls polyligne dont on souhaite calculer le ratio
        /// @param country pays d'appartenance de la polyligne
        /// @param fs1 classe d'objets surfaciques
        /// @param fs2 [optionnel] classe supplémentaire d'objets surfaciques
        /// @return 
        static double GetRatio(
            ign::geometry::LineString const& ls,
            std::string const& country,
            ign::feature::FeatureStore * fs1,
            ign::feature::FeatureStore * fs2 = 0
        );

        private:

        /// @brief fonction calculant la longeur des composantes linéaires d'un objet géométrique
        /// @param geom 
        /// @return 
        static double _GetLength(
            ign::geometry::Geometry const& geom
        );
        
    };

}
}
}

#endif