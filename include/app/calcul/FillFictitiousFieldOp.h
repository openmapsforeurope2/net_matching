#ifndef _APP_CALCUL_FILLFICTITIOUSFIELDOP_H_
#define _APP_CALCUL_FILLFICTITIOUSFIELDOP_H_

//SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>

//EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>


namespace app{
namespace calcul{

	/// @brief Classe pour le calcul du champ fictitious
	class FillFictitiousFieldOp {

        public:

            /// @brief Constructeur
            /// @param countryCode Code pays simple
            /// @param verbose Mode Verbeux
            FillFictitiousFieldOp( 
                std::string countryCode, 
                bool verbose 
            );

            /// @brief Destructeur
            ~FillFictitiousFieldOp();


            /// @brief Calcul le champ fictitious d'une classe d'objets linéaires en fonction de leur superposition ou non
            /// avec les objets de la (des) classe(s) de surfaciques associée(s). Un objet linéaire est considéré comme fictif
            /// si son ratio "longueur superposé/longueur totale" dépasse un certain seuil.
            /// @param countryCode Code pays simple
            /// @param verbose Mode verbeux
            static void Compute( 
                std::string countryCode, 
                bool verbose 
            );


        private:

            //--
            ign::feature::sql::FeatureStorePostgis*            _fsEdge;
            //--
            ign::feature::sql::FeatureStorePostgis*            _fsArea;
            //--
            ign::feature::sql::FeatureStorePostgis*            _fsStanding;
            //--
            epg::log::EpgLogger*                               _logger;
            //--
            epg::log::ShapeLogger*                             _shapeLogger;
            //--
            bool                                               _verbose;
            //--
            std::string                                        _countryCodeDouble;

        private:

            //--
            void _init();

            //--
            void _compute() const;

            //--
            double _getRatio(ign::geometry::LineString const& ls, std::string const& country) const;

            //--
            double _getLength(ign::geometry::Geometry const& geom) const;

    };

}
}

#endif