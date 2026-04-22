// APP
#include <app/calcul/FillFictitiousFieldOp.h>
#include <app/params/ThemeParameters.h>
#include <app/calcul/detail/LineStringAbsDampedDeformer.h>
#include <app/calcul/detail/RatioTools.h>

// BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

// EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/sql/tools/numFeatures.h>
#include <epg/tools/FilterTools.h>

// OME2
#include <ome2/feature/sql/NotDestroyedTools.h>



namespace app
{
    namespace calcul
    {

		///
        ///
        ///
        FillFictitiousFieldOp::FillFictitiousFieldOp( 
            bool verbose 
        ): 
            _verbose(verbose)
        {
            _init();
        }


        ///
        ///
        ///
        FillFictitiousFieldOp::~FillFictitiousFieldOp() {

        }


        ///
        ///
        ///
        void FillFictitiousFieldOp::Compute( 
            bool verbose 
        ) {
            FillFictitiousFieldOp op(verbose);
            op._compute();
        }


        ///
        ///
        ///
        void FillFictitiousFieldOp::_init() {
            //--
            _logger = epg::log::EpgLoggerS::getInstance();
            _logger->log(epg::log::INFO, "[START] initialization: " + epg::tools::TimeTools::getTime());

            //--
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const idName = epgParams.getValue(ID).toString();
            std::string const geomName = epgParams.getValue(GEOM).toString();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            std::string const edgeTableName = epgParams.getValue(EDGE_TABLE).toString();
            
            // app parameters
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            std::string areaTableName = themeParameters->getValue(WATERCOURSE_AREA_TABLE).toString();
            std::string standingTableName = themeParameters->getValue(STANDING_WATER_TABLE).toString();

            //--
            _fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);

            //--
            _fsAreaRam = context->getDataBaseManager().getFeatureStoreRam(areaTableName, idName, geomName); // veiller a ce que la table soit bien issu d'une extraction (et non table complete)
            ome2::feature::sql::NotDestroyedTools::RemoveDestroyed(_fsAreaRam);

            //--
            _fsStandingRam = context->getDataBaseManager().getFeatureStoreRam(standingTableName, idName, geomName); // veiller a ce que la table soit bien issu d'une extraction (et non table complete)
            ome2::feature::sql::NotDestroyedTools::RemoveDestroyed(_fsStandingRam);

            //--
            _logger->log(epg::log::INFO, "[END] initialization: " + epg::tools::TimeTools::getTime());
        }

        ///
        ///
        ///
        void FillFictitiousFieldOp::_compute() const {
            //--
            epg::Context *context = epg::ContextS::getInstance();

            //--
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            //--
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            double const minRatio = themeParameters->getValue(FFF_RATIO).toDouble();
            std::string const fictitiousFieldName = themeParameters->getValue(EDGE_FICTITIOUS_NAME).toString();

            //DEBUG
            std::string const wTagName = themeParameters->getParameter(W_TAG_NAME).getValue().toString();

            //--
            ign::feature::FeatureFilter filter;
            size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filter);
            boost::progress_display display(numFeatures, std::cout, "[ filling fictitious field % complete ]\n");

            ign::feature::FeatureIteratorPtr itEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filter);
            while (itEdge->hasNext())
            {
                ++display;
                ign::feature::Feature fEdge = itEdge->next();
                ign::geometry::LineString const& ls = fEdge.getGeometry().asLineString();
                std::string edgeId = fEdge.getId();
                std::string country = fEdge.getAttribute(countryCodeName).toString();
                std::string fictitious = fEdge.getAttribute(fictitiousFieldName).toString();

                double ratio = detail::RatioTools::GetRatio(ls, country, _fsAreaRam, _fsStandingRam);

                if (ratio >= minRatio && fictitious != "true") {
                    ign::feature::Feature fEdge_ = fEdge;
                    fEdge_.setAttribute(fictitiousFieldName, ign::data::String("true"));

                    //DEBUG
                    fEdge_.setAttribute(wTagName, ign::data::String("debug_201_to_true"));
                    
                    _fsEdge->modifyFeature(fEdge_);
                } else if (ratio < minRatio && fictitious != "false") {
                    ign::feature::Feature fEdge_ = fEdge;
                    fEdge_.setAttribute(fictitiousFieldName, ign::data::String("false"));

                    //DEBUG
                    fEdge_.setAttribute(wTagName, ign::data::String("debug_201_to_false"));
                    
                    _fsEdge->modifyFeature(fEdge_);
                }
            }
        }

    }
}
