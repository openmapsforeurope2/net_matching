#ifndef _APP_GEOMETRY_TOOLS_GETBUFFER_H_
#define _APP_GEOMETRY_TOOLS_GETBUFFER_H_


// SOCLE
#include <ign/geometry/algorithm/BufferOpGeos.h>


namespace app {
namespace geometry {
namespace tools {

    /// @brief Encapsulation de l'opérateur ign::geometry::algorithm::BufferOpGeos
    /// En cas d'échec du calcul, appelle la fonction PostGIS ST_Buffer.
    /// @param geom Géométrie à partir de laquelle le buffer est calculé
    /// @param distBuffer Rayon du buffer
    /// @param quadrantSegments Nombre de segments par quadrant
    /// @param endCapStyle Style aux extrémités
    /// @return Buffer
    ign::geometry::Geometry* getBuffer(
        ign::geometry::Geometry const& geom,
        double distBuffer,
        size_t quadrantSegments = 8,
		ign::geometry::algorithm::BufferOpGeos::EndCapStyle endCapStyle = ign::geometry::algorithm::BufferOpGeos::CAP_ROUND
    );
}
}
}

#endif