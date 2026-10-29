// APP
#include <app/calcul/ImportDataOutOfUpdateAreasOp.h>
#include <app/params/ThemeParameters.h>

// EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/sql/tools/numFeatures.h>
#include <epg/tools/FilterTools.h>

// OME2
#include <ome2/feature/sql/NotDestroyedTools.h>

// BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>


namespace app
{
    namespace calcul
    {

		///
        ///
        ///
        ImportDataOutOfUpdateAreasOp::ImportDataOutOfUpdateAreasOp( 
            bool verbose 
        ): 
            _verbose(verbose)
        {
            _init();
        }


        ///
        ///
        ///
        ImportDataOutOfUpdateAreasOp::~ImportDataOutOfUpdateAreasOp() {

        }


        ///
        ///
        ///
        void ImportDataOutOfUpdateAreasOp::Compute( 
            bool verbose 
        ) {
            ImportDataOutOfUpdateAreasOp op(verbose);
            op._compute();
        }


        ///
        ///
        ///
        void ImportDataOutOfUpdateAreasOp::_init() {
            //--
            _logger = epg::log::EpgLoggerS::getInstance();
            _logger->log(epg::log::INFO, "[START] initialization: " + epg::tools::TimeTools::getTime());

            //--
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const idName = epgParams.getValue(ID).toString();
            std::string const geomName = epgParams.getValue(GEOM).toString();
            std::string const edgeTableName = epgParams.getValue(EDGE_TABLE).toString();
            
            // app parameters
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            std::string const ofuaTableName = themeParameters->getValue(OFUA_TABLE).toString();

            //--
            _fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);

            //--
            _fsOfua = context->getDataBaseManager().getFeatureStore(ofuaTableName, idName, geomName);

            //--
            _logger->log(epg::log::INFO, "[END] initialization: " + epg::tools::TimeTools::getTime());
        }

        ///
        ///
        ///
        void ImportDataOutOfUpdateAreasOp::_compute() const {
			epg::Context* context = epg::ContextS::getInstance();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

            ign::feature::FeatureFilter filter;
            ign::feature::FeatureIteratorPtr itOfua = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsOfua, filter);

			while (itOfua->hasNext())
            {
                ign::feature::Feature fOfua = itOfua->next();
				_fsEdge->createFeature(fOfua);
			}
        }
    }
}
