#ifndef _APP_TOOLS_STRINGTOOLS_H_
#define _APP_TOOLS_STRINGTOOLS_H_

// SOCLE
#include <ign/geometry/graph/GeometryGraph.h>
#include <ign/geometry/graph/tools/SnapRoundPlanarizer.h>
#include <ign/geometry/graph/builder/SimpleGraphBuilder.h>

namespace app{
namespace tools{

    /// @brief Classe utilitaire pour les opérations sur les chaînes de charactères
    class StringTools {
        public:

		/// @brief Transforme une liste en chaîne de charactères
		/// @tparam Container 
		/// @param vStrings Liste de chaînes de charactères
		/// @param separator Séparateur
		/// @return Chaîne de charactères
		template<typename Container>
		static std::string ToString(Container const& vStrings, std::string separator = ",") {
            std::string result = "";
            typename Container::const_iterator vit;
            for ( vit = vStrings.begin() ; vit != vStrings.end() ; ++vit ) {
                if ( !result.empty() ) result += separator;
                result += *vit;
            }
            return result;
        };
    };

}
}

#endif



        