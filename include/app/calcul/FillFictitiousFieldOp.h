#ifndef _APP_CALCUL_FILLFICTITIOUSFIELDOP_H_
#define _APP_CALCUL_FILLFICTITIOUSFIELDOP_H_

//SOCLE
#include <ign/feature/sql/FeatureStorePostgis.h>

//EPG
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>


namespace app{
namespace calcul{

	/// @brief 
	class FillFictitiousFieldOp {

        public:

            /// @brief 
            /// @param countryCode 
            /// @param verbose 
            FillFictitiousFieldOp( 
                std::string countryCode, 
                bool verbose 
            );

            /// @brief 
            ~FillFictitiousFieldOp();


            /// @brief 
            /// @param countryCode 
            /// @param verbose 
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