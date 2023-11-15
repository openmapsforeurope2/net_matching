// APP
#include <app/calcul/detail/graph/concept/EdgeCleaningGraphSpecializations.h>


namespace ign{
namespace graph{
namespace concept{
    //specialisation pour le concept direction
    template <>
    AllowedDirection direction( app::calcul::detail::graph::EdgeCleaningEdge const& properties )
    {
        return ALLOWED_BOTH;
    }

}
}
}

