#include <app/calcul/FillFictitiousFieldOp.h>

//APP
#include <app/params/ThemeParameters.h>

//BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

//EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/sql/tools/numFeatures.h>


namespace app
{
    namespace calcul
    {

		///
        ///
        ///
        FillFictitiousFieldOp::FillFictitiousFieldOp( 
            std::string countryCode, 
            bool verbose 
        ) {
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
            std::string countryCode, 
            bool verbose 
        ) {
            FillFictitiousFieldOp op(countryCode, verbose);
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
            // _shapeLogger = epg::log::ShapeLoggerS::getInstance();
            // _shapeLogger->addShape("ecl_deleted_edges", epg::log::ShapeLogger::LINESTRING);

            //--
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            // std::string const boundaryTableName = epgParams.getValue(TARGET_BOUNDARY_TABLE).toString();
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
            _fsArea = context->getDataBaseManager().getFeatureStore(areaTableName, idName, geomName);

            //--
            _fsStanding = context->getDataBaseManager().getFeatureStore(standingTableName, idName, geomName);

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
            std::string const fictitiousFieldName = themeParameters->getValue(EDGE_FICTITIOUS).toString();

            //--
            ign::feature::FeatureFilter filter;
            int numFeatures = epg::sql::tools::numFeatures(*_fsEdge, filter);
            boost::progress_display display(numFeatures, std::cout, "[ filling fictitious field  % complete ]\n");

            ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(filter);
            while (itEdge->hasNext())
            {
                ++display;
                ign::feature::Feature const& fEdge = itEdge->next();
                ign::geometry::LineString const& ls = fEdge.getGeometry().asLineString();
                std::string edgeId = fEdge.getId();
                std::string country = fEdge.getAttribute(countryCodeName).toString();
                std::string fictitious = fEdge.getAttribute(fictitiousFieldName).toString();

                double ratio = _getRatio(ls);

                if (ratio > minRatio && fictitious == "false") {
                    ign::feature::Feature fEdge_ = fEdge;
                    fEdge_.setAttribute(fictitiousFieldName, ign::data::String("true"));
                    _fsEdge->modifyFeature(fEdge_);
                }
            }
        }

        ///
        ///
        ///
        double FillFictitiousFieldOp::_getRatio(ign::geometry::LineString const& ls) const {

        }

        ///
        ///
        ///
        std::pair<double, double> FillFictitiousFieldOp::_addLengths(
            std::string country, 
            ign::geometry::LineString const& ls,
            double & lengthInCountry,
            double & length
        ) const {
            // std::map<std::string, ign::geometry::GeometryPtr>::const_iterator mit = _mCountryGeomPtr.find(country);
            // if (mit == _mCountryGeomPtr.end()) {
            //     _logger->log(epg::log::ERROR, "Unknown country [country code] " + country);
            //     return std::make_pair(0, 0);
            // }

            // ign::geometry::GeometryPtr resultPtr(mit->second->Intersection(ls));

            // std::pair<double, double> lengths = _getLengths(*resultPtr, &ls.startPoint());

            // lengthInCountry += lengths.first;
            // length += ls.length();

            // return std::make_pair(lengths.second, lengths.second/ls.length());
        }

        ///
        std::pair<double, double> FillFictitiousFieldOp::_getLengths( ign::geometry::Geometry const& geom, ign::geometry::Point const* startPoint ) const
        {
            double length = 0;
            double lengthFirstPart = 0;

            ign::geometry::Geometry::GeometryType geomType = geom.getGeometryType();
            switch( geomType )
            {
                case ign::geometry::Geometry::GeometryTypeNull :
                case ign::geometry::Geometry::GeometryTypePoint :
                case ign::geometry::Geometry::GeometryTypeMultiPoint :
                case ign::geometry::Geometry::GeometryTypeTriangle :
                case ign::geometry::Geometry::GeometryTypeTriangulatedSurface :
                case ign::geometry::Geometry::GeometryTypePolyhedralSurface :
                case ign::geometry::Geometry::GeometryTypePolygon :
                case ign::geometry::Geometry::GeometryTypeMultiPolygon :
                    return std::make_pair(0, 0);
                case ign::geometry::Geometry::GeometryTypeLineString :
                    {
                        ign::geometry::LineString const& ls = geom.asLineString();
                        if (ls.isEmpty()) return std::make_pair(0, 0);
                        double length = ls.length();
                        return std::make_pair(length, length);
                    }
                    
                case ign::geometry::Geometry::GeometryTypeMultiLineString : 
                    {
                        ign::geometry::MultiLineString const& mls = geom.asMultiLineString();
                        for( size_t i = 0 ; i < mls.numGeometries() ; ++i ) {
                            length += mls.lineStringN(i).length();
                            if ( startPoint != 0 && (startPoint->distance(mls.lineStringN(i).startPoint()) < 1e-5 || startPoint->distance(mls.lineStringN(i).endPoint()) < 1e-5))
                                lengthFirstPart += length;
                        }
                        return std::make_pair(length, lengthFirstPart);
                    }
                
                case ign::geometry::Geometry::GeometryTypeGeometryCollection :
                    {
                        ign::geometry::GeometryCollection const& collection = geom.asGeometryCollection();
                        for( size_t i = 0 ; i < collection.numGeometries() ; ++i ) {
                            std::pair<double, double> lengths = _getLengths( collection.geometryN(i), startPoint);
                            length += lengths.first;
                            lengthFirstPart += lengths.second;
                        }
                    
                        return std::make_pair(length, lengthFirstPart);
                    }
                default :
                    IGN_THROW_EXCEPTION( "Geometry type not allowed" );
            }
        }


    }
}
