#ifndef _APP_UTILS_CREATECPCLTABLES_H_
#define _APP_UTILS_CREATECPCLTABLES_H_

#include <string>

namespace app{
namespace utils{

    /// @brief Fonction pour la création des tables pour le stockage
    /// des connecting lines et des connecting points
    /// @param edgeTableName Nom de la table du réseau
    void createCpClTables(std::string const& edgeTableName);
}
}

#endif