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
#include <epg/tools/FilterTools.h>


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

            //DEBUG
            std::string const wTagName = themeParameters->getParameter(W_TAG).getValue().toString();

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

                double ratio = _getRatio(ls, country);

                if (ratio > minRatio && fictitious == "false") {
                    ign::feature::Feature fEdge_ = fEdge;
                    fEdge_.setAttribute(fictitiousFieldName, ign::data::String("true"));

                    //DEBUG
                    fEdge_.setAttribute(wTagName, ign::data::String("debug_201"));
                    
                    _fsEdge->modifyFeature(fEdge_);
                }
            }
        }

        ///
        ///
        ///
        double FillFictitiousFieldOp::_getRatio(ign::geometry::LineString const& ls, std::string const& country) const {
            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const geomName = epgParams.getValue(GEOM).toString();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            ign::geometry::GeometryPtr areaUnionPtr(new ign::geometry::Polygon());

            ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + ls.toString() + "'),3035))");
            epg::tools::FilterTools::addAndConditions(filter, countryCodeName +" = '"+country+"'");
            ign::feature::FeatureIteratorPtr itArea = _fsArea->getFeatures(filter);
            while (itArea->hasNext())
            {
                ign::feature::Feature const& fArea = itArea->next();
                ign::geometry::MultiPolygon const& areaGeom = fArea.getGeometry().asMultiPolygon();

                areaUnionPtr.reset(areaUnionPtr->Union(areaGeom));
            }
            ign::feature::FeatureIteratorPtr itStand = _fsStanding->getFeatures(filter);
            while (itStand->hasNext())
            {
                ign::feature::Feature const& fStand = itStand->next();
                ign::geometry::MultiPolygon const& standGeom = fStand.getGeometry().asMultiPolygon();
                
                areaUnionPtr.reset(areaUnionPtr->Union(standGeom));
            }

            if(areaUnionPtr->isEmpty() || areaUnionPtr->isNull()) return 0;

            ign::geometry::GeometryPtr resultPtr(areaUnionPtr->Intersection(ls));

            double lengthInter = _getLength(*resultPtr);

            return lengthInter / ls.length();

        }

        ///
        double FillFictitiousFieldOp::_getLength( ign::geometry::Geometry const& geom ) const
        {
            double length = 0;

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
                    return 0;
                case ign::geometry::Geometry::GeometryTypeLineString :
                    {
                        return geom.asLineString().length();
                    }
                    
                case ign::geometry::Geometry::GeometryTypeMultiLineString : 
                    {
                        ign::geometry::MultiLineString const& mls = geom.asMultiLineString();
                        for( size_t i = 0 ; i < mls.numGeometries() ; ++i ) {
                            length += mls.lineStringN(i).length();
                        }
                        return length;
                    }
                
                case ign::geometry::Geometry::GeometryTypeGeometryCollection :
                    {
                        ign::geometry::GeometryCollection const& collection = geom.asGeometryCollection();
                        for( size_t i = 0 ; i < collection.numGeometries() ; ++i ) {
                            length += _getLength(collection.geometryN(i));
                        }
                        return length;
                    }
                default :
                    IGN_THROW_EXCEPTION( "Geometry type not allowed" );
            }
        }


    }
}
