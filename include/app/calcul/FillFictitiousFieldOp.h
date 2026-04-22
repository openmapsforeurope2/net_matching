#ifndef _APP_CALCUL_FILLFICTITIOUSFIELDOP_H_
#define _APP_CALCUL_FILLFICTITIOUSFIELDOP_H_

// SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>
#include <ign/feature/ram/FeatureStoreRam.h>

// EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>


namespace app{
namespace calcul{

	/// @brief Classe pour le calcul du champ fictitious
	class FillFictitiousFieldOp {

        public:

            /// @brief Constructeur
            /// @param verbose Mode Verbeux
            FillFictitiousFieldOp(
                bool verbose 
            );

            /// @brief Destructeur
            ~FillFictitiousFieldOp();


            /// @brief Calcul le champ fictitious d'une classe d'objets linéaires en fonction de leur superposition ou non
            /// avec les objets de la (des) classe(s) de surfaciques associée(s). Un objet linéaire est considéré comme fictif
            /// si son ratio "longueur superposé/longueur totale" dépasse un certain seuil.
            /// @param verbose Mode verbeux
            static void Compute( 
                bool verbose 
            );


        private:

            //--
            ign::feature::sql::FeatureStorePostgis*            _fsEdge;
            //--
            ign::feature::ram::FeatureStoreRam*                _fsAreaRam;
            //--
            ign::feature::ram::FeatureStoreRam*                _fsStandingRam;
            //--
            epg::log::EpgLogger*                               _logger;
            //--
            epg::log::ShapeLogger*                             _shapeLogger;
            //--
            bool                                               _verbose;

        private:

            //--
            void _init();

            //--
            void _compute() const;

    };

}
}

#endif