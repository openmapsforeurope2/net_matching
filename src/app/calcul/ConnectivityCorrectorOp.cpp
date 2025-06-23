// APP
#include <app/calcul/ConnectivityCorrectorOp.h>
#include <app/params/ThemeParameters.h>
#include <app/calcul/detail/graph/concept/EdgeCleaningGraphSpecializations.h>
#include <app/geometry/tools/LineStringSplitter.h>

// BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

// EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/sql/tools/numFeatures.h>
#include <epg/tools/geometry/project.h>
#include <epg/tools/StringTools.h>

// OME2
#include <ome2/feature/sql/NotDestroyedTools.h>

// STL
#include <sstream>


namespace app
{
    namespace calcul
    {

		///
        ///
        ///
        void ConnectivityCorrectorOp::Compute(
            std::string const& borderCode,
            bool verbose
        ) {
            ConnectivityCorrectorOp op(verbose);

            std::vector<std::string> vCountry;
            epg::tools::StringTools::Split(borderCode, "#", vCountry);
            for (size_t i = 0 ; i < vCountry.size() ; ++i) {
                op._compute(vCountry[i]);
            }
        }

		///
        ///
        ///
        ConnectivityCorrectorOp::ConnectivityCorrectorOp(
            bool verbose
        ) : _verbose(verbose)
        {
            _init();
        }

		///
		///
		///
		ConnectivityCorrectorOp::~ConnectivityCorrectorOp()
		{
            _shapeLogger->closeShape("disconnection");
            _shapeLogger->closeShape("missing_connection");
		}

		///
        ///
        ///
        void ConnectivityCorrectorOp::_init()
        {
			//--
            _logger = epg::log::EpgLoggerS::getInstance();
            _logger->log(epg::log::INFO, "[START] initialization: " + epg::tools::TimeTools::getTime());

            //--
            _shapeLogger = epg::log::ShapeLoggerS::getInstance();
            _shapeLogger->addShape("disconnection", epg::log::ShapeLogger::POINT);
            _shapeLogger->addShape("missing_connection", epg::log::ShapeLogger::POINT);

            //--
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const edgeTableName = epgParams.getValue(EDGE_TABLE).toString();
            std::string const idName = epgParams.getValue(ID).toString();
            std::string const geomName = epgParams.getValue(GEOM).toString();

			//--
            _fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);

            //--
            _logger->log(epg::log::INFO, "[END] initialization: " + epg::tools::TimeTools::getTime());
		}

        ///
        ///
        ///
        void ConnectivityCorrectorOp::_compute(
            std::string const& countryCode
        ) const {
            _connectAtEndings(countryCode);
            _createMissingConnections(countryCode);
        }

        ///
        ///
        ///
        void ConnectivityCorrectorOp::_createMissingConnections(
            std::string const& countryCode
        ) const {
            //--
            epg::Context* context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            //--
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const distThreshold = themeParameters->getValue(CC_DIST_THRESHOLD).toDouble();
            std::string const natIdName = themeParameters->getParameter(NATIONAL_IDENTIFIER_NAME).getValue().toString();

            //--
            ign::geometry::index::QuadTree< std::string > qTreeEdges;
            std::map<std::string, ign::geometry::LineString> mEdges;
            std::map<ign::geometry::Point, std::list<std::string>, PointCompare> mVertices;

            //--
            ign::feature::FeatureFilter filterEdge(countryCodeName + " = '" + countryCode + "'");

            size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterEdge);
            boost::progress_display display(numFeatures, std::cout, "[ edge loading % complete ]\n");

            ign::feature::FeatureIteratorPtr itEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterEdge);
            while (itEdge->hasNext())
            {
                ++display;
                ign::feature::Feature const& fEdge = itEdge->next();
                ign::geometry::LineString const& ls = fEdge.getGeometry().asLineString();
                std::string edgeId = fEdge.getId();

                qTreeEdges.insert(edgeId, ls.getEnvelope());
                mEdges.insert(std::make_pair(edgeId, ls));

                if ( mVertices.find(ls.startPoint()) == mVertices.end() )
                    mVertices.insert(std::make_pair(ls.startPoint(), std::list<std::string>()));
                mVertices[ls.startPoint()].push_back(edgeId);

                if ( mVertices.find(ls.endPoint()) == mVertices.end() )
                    mVertices.insert(std::make_pair(ls.endPoint(), std::list<std::string>()));
                mVertices[ls.endPoint()].push_back(edgeId);
            }

            std::map<std::string, std::vector<ign::geometry::Point>> mEdgeSplittingPoints;
            for( std::map<ign::geometry::Point, std::list<std::string>>::const_iterator mit = mVertices.begin() ; mit != mVertices.end() ; ++mit) {
                
                if( mit->second.size() > 1 ) //not dangle
                    continue;
                
                std::set<std::string> sEdges;
		        qTreeEdges.query( mit->first.getEnvelope().expandBy(distThreshold), sEdges );

                for( std::set<std::string>::const_iterator sit = sEdges.begin() ; sit != sEdges.end() ; ++sit ) {
                    if( *sit == *mit->second.begin() ) 
                        continue;

                    if( mEdges[*sit].distance(mit->first) < distThreshold) {
                        if ( mEdgeSplittingPoints.find(*sit) == mEdgeSplittingPoints.end() )
                            mEdgeSplittingPoints.insert(std::make_pair(*sit, std::vector<ign::geometry::Point>()));
                        mEdgeSplittingPoints[*sit].push_back(mit->first);

                        ign::feature::Feature feat;
                        feat.setGeometry(mit->first);
                        _shapeLogger->writeFeature("missing_connection", feat);
                    }
                }
            }

            std::map<std::string, std::vector<ign::geometry::Point>>::const_iterator mit;
            for( mit = mEdgeSplittingPoints.begin() ; mit != mEdgeSplittingPoints.end() ; ++mit ) {
                ign::feature::Feature fEdge;
                _fsEdge->getFeatureById(mit->first, fEdge);
                ign::geometry::LineString const& ls = fEdge.getGeometry().asLineString();
                std::string const natId = fEdge.getAttribute(natIdName).toString();

                app::geometry::tools::LineStringSplitter lsSplitter(ls);
                for (size_t i = 0 ; i < mit->second.size() ; ++i) {
                    ign::geometry::Point proj = epg::tools::geometry::project(ls,  mit->second[i]);
                    lsSplitter.addCuttingGeometry(proj);
                }
                std::vector<ign::geometry::LineString> subEdges = lsSplitter.getSubLineStringsZ();

                size_t count = 1;
                for (size_t i = 0 ; i < subEdges.size() ; ++i) {
                    if ( i != 0 ) {
                        subEdges[i].startPoint() = _getClosest(subEdges[i].startPoint(), mit->second);
                    }
                    if ( i != subEdges.size()-1 ) {
                        subEdges[i].endPoint() = _getClosest(subEdges[i].endPoint(), mit->second);
                    }

                    std::ostringstream ss;
                    ss << natId << "_" << count++;

                    fEdge.setAttribute(natIdName, ign::data::String(ss.str()) );
                    fEdge.setGeometry(subEdges[i]);
                    _fsEdge->createFeature(fEdge);
                }
                _fsEdge->deleteFeature(mit->first);
            }
        }

        ///
        ///
        ///
        void ConnectivityCorrectorOp::_connectAtEndings(
            std::string const& countryCode
        ) const {
            //--
            epg::Context* context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            //--
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const distThreshold = themeParameters->getValue(CC_DIST_THRESHOLD).toDouble();

            //--
            std::vector<std::pair<std::string, ENDING>> vVertices;
            ign::geometry::index::QuadTree< int > qTreeVertices;
            std::map<std::string, ign::geometry::LineString> mEdges;
            std::set<std::string> sModifiedEdges;
            std::set<int> sTreatedVertices;

            //--
            ign::feature::FeatureFilter filterEdge(countryCodeName + " = '" + countryCode + "'");

            size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterEdge);
            boost::progress_display display(numFeatures, std::cout, "[ vertices loading % complete ]\n");

            ign::feature::FeatureIteratorPtr itEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterEdge);
            while (itEdge->hasNext())
            {
                ++display;
                ign::feature::Feature const& fEdge = itEdge->next();
                ign::geometry::LineString const& ls = fEdge.getGeometry().asLineString();
                std::string edgeId = fEdge.getId();

                qTreeVertices.insert(vVertices.size(), ls.startPoint().getEnvelope());
                vVertices.push_back( std::make_pair(edgeId, START));

                qTreeVertices.insert(vVertices.size(), ls.endPoint().getEnvelope());
                vVertices.push_back( std::make_pair(edgeId, END));
                
                mEdges.insert( std::make_pair(edgeId, ls) );
            }

            for ( size_t i = 0 ; i < vVertices.size() ; ++i )
            {
                if ( sTreatedVertices.find(i) != sTreatedVertices.end() )
                    continue;

                ign::geometry::Point const& refVertexGeom = vVertices[i].second == START ? mEdges[vVertices[i].first].startPoint() : mEdges[vVertices[i].first].endPoint();

                std::set<int> sVertices;
		        qTreeVertices.query( refVertexGeom.getEnvelope().expandBy(distThreshold), sVertices );

                for (std::set<int>::const_iterator sit = sVertices.begin() ; sit != sVertices.end() ; ++sit )
                {
                    if( *sit == i )
                        continue;

                    if ( sTreatedVertices.find(*sit) != sTreatedVertices.end() )
                        continue;

                    ign::geometry::Point const& vertexGeom = vVertices[*sit].second == START ? mEdges[vVertices[*sit].first].startPoint() : mEdges[vVertices[*sit].first].endPoint();

                    double dist = refVertexGeom.distance(vertexGeom);
                    if (dist > 0 && dist < distThreshold) {
                        sModifiedEdges.insert(vVertices[*sit].first);
                        sTreatedVertices.insert(*sit);

                        ign::geometry::Point* endingPointToModifyPtr = vVertices[*sit].second == START ? &mEdges[vVertices[*sit].first].startPoint() : &mEdges[vVertices[*sit].first].endPoint();

                        ign::feature::Feature feat;
                        feat.setGeometry(*endingPointToModifyPtr);
                        _shapeLogger->writeFeature("disconnection", feat);

                        *endingPointToModifyPtr = refVertexGeom;
                    }
                }
            }

            for (std::set<std::string>::const_iterator sit = sModifiedEdges.begin(); sit != sModifiedEdges.end() ; ++sit) {
                ign::feature::Feature fEdge;
                _fsEdge->getFeatureById( *sit, fEdge );
                fEdge.setGeometry(mEdges[*sit]);
                _fsEdge->modifyFeature( fEdge );
            }
        }

        ///
        ///
        ///
        ign::geometry::Point ConnectivityCorrectorOp::_getClosest(
            ign::geometry::Point const& pt, 
            std::vector<ign::geometry::Point> const& vPt
        ) const {
            double minDist = std::numeric_limits<double>::max();
            int index = -1;
            for ( size_t i = 0 ; i < vPt.size() ; ++i ) {
                double dist = pt.distance(vPt[i]);
                if( dist < minDist) {
                    minDist = dist;
                    index = i;
                }
            }
            return vPt[index];
        }
    }
}