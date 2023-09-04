// APP
#include <app/calcul/CFeatConnectionOp.h>
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LineStringSplitter.h>

// BOOST
#include <boost/progress.hpp>
#include <boost/optional.hpp>

// EPG
#include <epg/Context.h>
#include <epg/params/EpgParameters.h>
#include <epg/sql/tools/numFeatures.h>
#include <epg/sql/DataBaseManager.h>
#include <epg/tools/StringTools.h>
#include <epg/tools/TimeTools.h>
#include <epg/tools/geometry/project.h>
#include <epg/calcul/matching/detail/LineStringSimpleDampedDeformer.h>

// SOCLE
#include <ign/geometry/graph/builder/SimpleGraphBuilder.h>
#include <ign/tools/stringtools.h>

namespace app
{
    namespace calcul
    {

        ///
        ///
        ///
        void CFeatConnectionOp::computeCp(
            std::string edgeTable,
            std::string cpTable,
            std::string countryCode,
            bool verbose)
        {
            CFeatConnectionOp CFeatConnectionOp(edgeTable, cpTable, "", countryCode, verbose);
            CFeatConnectionOp._computeCp();
        }

        ///
        ///
        ///
        void CFeatConnectionOp::computeCl(
            std::string edgeTable,
            std::string clTable,
            std::string countryCode,
            bool verbose)
        {
            CFeatConnectionOp CFeatConnectionOp(edgeTable, "", clTable, countryCode, verbose);
            CFeatConnectionOp._computeCl();
        }

        ///
        ///
        ///
        void CFeatConnectionOp::computeCpCl(
            std::string edgeTable,
            std::string cpTable,
            std::string clTable,
            std::string countryCode,
            bool verbose)
        {
            CFeatConnectionOp CFeatConnectionOp(edgeTable, cpTable, clTable, countryCode, verbose);
            CFeatConnectionOp._computeCpCl();
        }

        ///
        ///
        ///
        CFeatConnectionOp::CFeatConnectionOp(
            std::string edgeTable,
            std::string cpTable,
            std::string clTable,
            std::string countryCode,
            bool verbose
        ) : _countryCode(countryCode),
            _verbose(verbose)
        {
            _init(edgeTable, cpTable, clTable);
        }

        ///
        ///
        ///
        CFeatConnectionOp::~CFeatConnectionOp()
        {
            _shapeLogger->closeShape("resulting_edges");
            _shapeLogger->closeShape("projected_cp");
            _shapeLogger->closeShape("multiple_result");
            _shapeLogger->closeShape("cp_displacements");
            _shapeLogger->closeShape("created_edges");
            _shapeLogger->closeShape("cl_displacements");
            _shapeLogger->closeShape("cl_created_features");
            _shapeLogger->closeShape("cl_deleted_features");

            epg::log::ShapeLoggerS::kill();
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_init(
            std::string edgeTable,
            std::string cpTable,
            std::string clTable)
        {
            //--
            _logger = epg::log::EpgLoggerS::getInstance();
            _logger->log(epg::log::INFO, "[START] initialization: " + epg::tools::TimeTools::getTime());

            //--
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const &epgParams = context->getEpgParameters();

            std::string const idName = epgParams.getValue(ID).toString();
            std::string const geomName = epgParams.getValue(GEOM).toString();

            // app parameters
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            std::string const landmaskTableName = themeParameters->getValue(LANDMASK_TABLE).toString();

            //--
            _fsLandmask = context->getDataBaseManager().getFeatureStore(landmaskTableName, idName, geomName);

            //--
            if (cpTable != "")
            {
                _fsCp = context->getDataBaseManager().getFeatureStore(cpTable, idName, geomName);
            }

            //--
            if (clTable != "")
            {
                _fsCl = context->getDataBaseManager().getFeatureStore(clTable, idName, geomName);
            }

            //--
            _fsEdge = context->getDataBaseManager().getFeatureStore(edgeTable, idName, geomName);

            //--
            _shapeLogger = epg::log::ShapeLoggerS::getInstance();
            _shapeLogger->setDataDirectory(context->getLogDirectory());
            _shapeLogger->addShape("resulting_edges", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("projected_cp", epg::log::ShapeLogger::POINT);
            _shapeLogger->addShape("multiple_result", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("cp_displacements", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("created_edges", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("cl_displacements", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("cl_created_features", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("cl_deleted_features", epg::log::ShapeLogger::LINESTRING);

            //--
            _logger->log(epg::log::INFO, "[END] initialization: " + epg::tools::TimeTools::getTime());
        };

        ///
        ///
        ///
        void CFeatConnectionOp::_computeClDisplacements(std::map<ign::geometry::Point, ign::math::Vec2d> & mDisplacements) const
        {
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const &epgParams = context->getEpgParameters();

            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            // app params
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();

            std::string const edgeLinkName = themeParameters->getValue(EDGE_LINK).toString();
            std::string const landCoverTypeName = themeParameters->getValue(LAND_COVER_TYPE).toString();
            std::string const landAreaValue = themeParameters->getValue(TYPE_LAND_AREA).toString();
            double snapDistance = themeParameters->getValue(SNAP_DIST).toDouble();

            ign::feature::FeatureFilter filterCl(countryCodeName + " LIKE '%" + _countryCode + "%'");
            std::set<std::string> sEdge2Delete;

            // patience
            int numFeatures = epg::sql::tools::numFeatures(*_fsCl, filterCl);
            boost::progress_display display(numFeatures, std::cout, "[ cl_connection  % complete ]\n");

            ign::feature::FeatureIteratorPtr itCl = _fsCl->getFeatures(filterCl);
            while (itCl->hasNext())
            {
                ign::feature::Feature const &fCl = itCl->next();
                ign::geometry::LineString const &clGeom = fCl.getGeometry().asLineString();
                std::string edgeLink = fCl.getAttribute(edgeLinkName).toString();
                std::string countryCode = fCl.getAttribute(countryCodeName).toString();

                if (_verbose) _logger->log(epg::log::DEBUG, fCl.getId());

                // std::string idDebug = fCl.getId();
                // // if ( idDebug != "CONNECTINGLINE179" && idDebug != "CONNECTINGLINE178" ) continue;
                // if ( idDebug != "CONNECTINGLINE177" ) continue;

                std::pair<bool, std::string> foundEdgeLink = _getCountryEdgeLink(edgeLink, countryCode, _countryCode);
                if (!foundEdgeLink.first)
                {
                    _logger->log(epg::log::ERROR, "Edge link not found [connecting edge id] " + fCl.getId());
                    continue;
                }

                std::pair<bool, std::string> foundEdge = _getNearestEdge(clGeom, countryCodeName, edgeLinkName, foundEdgeLink.second, sEdge2Delete);
                if (!foundEdge.first)
                {
                    _logger->log(epg::log::ERROR, "No candidate edge found [" + edgeLinkName + "] " + foundEdgeLink.second);
                    continue;
                }

                // if (foundEdge.second ==  "50ab969e-2186-4c19-a0bb-ddc2d8ada2b8" || foundEdge.second ==  "38cb0d6f-50cf-4ce3-ae74-9d6f8a5915d4") {
                //     bool test =true;
                // }

                ign::feature::Feature fEdge;
                _fsEdge->getFeatureById(foundEdge.second, fEdge);

                if (_verbose) _logger->log(epg::log::DEBUG, "  "+foundEdge.second);

                ign::geometry::LineString const &edgeGeom = fEdge.getGeometry().asLineString();

                // peut-on connecter l'edge aux extremités du connecting edge ?
                std::vector<double> vDist;
                vDist.push_back(clGeom.startPoint().distance(edgeGeom.startPoint()));
                vDist.push_back(clGeom.startPoint().distance(edgeGeom.endPoint()));
                vDist.push_back(clGeom.endPoint().distance(edgeGeom.startPoint()));
                vDist.push_back(clGeom.endPoint().distance(edgeGeom.endPoint()));

                size_t minId = 0;
                for (size_t i = 1; i < vDist.size(); ++i)
                {
                    if (vDist[i] < vDist[minId])
                        minId = i;
                }

                ign::geometry::MultiPoint mpCuttingPoints;
                boost::optional< ign::geometry::Point > endingPoint;
                if (minId < 2)
                {
                    ign::geometry::Point startPoint;
                    if (vDist[minId] <= snapDistance)
                    {
                        startPoint = (minId == 0) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                        endingPoint = startPoint;
                    }
                    else
                    {
                        epg::tools::geometry::projectZ(edgeGeom, clGeom.startPoint(), startPoint);
                        mpCuttingPoints.addGeometry(startPoint);
                    }
                    mDisplacements.insert(std::make_pair(startPoint, clGeom.startPoint().toVec2d() - startPoint.toVec2d()));

                    // size_t minId2 = vDist[2] < vDist[3] ? 2 : 3;
                    size_t minId2 = minId == 0  ? 3 : 2;
                    if (vDist[minId2] <= snapDistance)
                    {
                        startPoint = (minId2 == 2) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                        endingPoint = startPoint;
                    }
                    else
                    {
                        epg::tools::geometry::projectZ(edgeGeom, clGeom.endPoint(), startPoint);
                        mpCuttingPoints.addGeometry(startPoint);
                    }
                    mDisplacements.insert(std::make_pair(startPoint, clGeom.endPoint().toVec2d() - startPoint.toVec2d()));
                }
                else
                {
                    ign::geometry::Point startPoint;
                    if (vDist[minId] <= snapDistance)
                    {
                        startPoint = (minId == 2) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                        endingPoint = startPoint;
                    }
                    else
                    {
                        epg::tools::geometry::projectZ(edgeGeom, clGeom.endPoint(), startPoint);
                        mpCuttingPoints.addGeometry(startPoint);
                    }
                    mDisplacements.insert(std::make_pair(startPoint, clGeom.endPoint().toVec2d() - startPoint.toVec2d()));

                    // size_t minId2 = vDist[0] < vDist[1] ? 0 : 1;
                    size_t minId2 = minId == 2  ? 1 : 0;
                    if (vDist[minId2] <= snapDistance)
                    {
                        startPoint = (minId2 == 0) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                        endingPoint = startPoint;
                    }
                    else
                    {
                        epg::tools::geometry::projectZ(edgeGeom, clGeom.startPoint(), startPoint);
                        mpCuttingPoints.addGeometry(startPoint);
                    }
                    mDisplacements.insert(std::make_pair(startPoint, clGeom.startPoint().toVec2d() - startPoint.toVec2d()));
                }

                std::vector< ign::geometry::LineString > vNewGeom;
                if (!mpCuttingPoints.isEmpty())
                {
                    app::geometry::tools::LineStringSplitter lsSplitter(edgeGeom, 1e-5);
                    lsSplitter.addCuttingGeometry(mpCuttingPoints);

                    vNewGeom = lsSplitter.getSubLineStringsZ();

                    if ( endingPoint ) mpCuttingPoints.addGeometry(endingPoint.get());
                    for (int i = vNewGeom.size()-1 ; i >= 0 ; --i) {
                        if(vNewGeom[i].startPoint().distance(mpCuttingPoints) < 1e-5 && vNewGeom[i].endPoint().distance(mpCuttingPoints) < 1e-5) {
                            vNewGeom.erase(vNewGeom.begin()+i);
                        }
                    }
                }
                if (vNewGeom.empty()) {
                    sEdge2Delete.insert(foundEdge.second);
                } else if (vNewGeom.size() == 1) {
                    fEdge.setGeometry(vNewGeom.front());
                    _fsEdge->modifyFeature(fEdge);
                } else {
                    sEdge2Delete.insert(foundEdge.second);

                    for (size_t i = 0 ; i < vNewGeom.size() ; ++i) {
                        fEdge.setGeometry(vNewGeom[i]);
                        _fsEdge->createFeature(fEdge);

                        vNewGeom[i].clearZ();
                        fEdge.setGeometry(vNewGeom[i]);
                        _shapeLogger->writeFeature("cl_created_features", fEdge);
                    }
                }
                ++display;
            }

            std::set< std::string >::const_iterator sit;
            for( sit = sEdge2Delete.begin() ; sit != sEdge2Delete.end() ; ++sit ) {
                ign::feature::Feature dFeat;
                _fsEdge->getFeatureById(*sit, dFeat);

                ign::geometry::LineString tempLs = dFeat.getGeometry().asLineString();
                tempLs.clearZ();
                dFeat.setGeometry(tempLs);
                _shapeLogger->writeFeature("cl_deleted_features", dFeat);

                _fsEdge->deleteFeature(*sit);
            }

            std::map<ign::geometry::Point, ign::math::Vec2d>::const_iterator mit;
            for (mit = mDisplacements.begin(); mit != mDisplacements.end(); ++mit)
            {
                ign::feature::Feature feat;
                ign::geometry::Point p = mit->first;
                p.clearZ();
                feat.setGeometry(ign::geometry::LineString(p, ign::geometry::Point(mit->first.x() + mit->second.x(), mit->first.y() + mit->second.y())));
                _shapeLogger->writeFeature("cl_displacements", feat);
            }
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_computeCl()
        {
            std::map<ign::geometry::Point, ign::math::Vec2d> mDisplacements;
            _computeClDisplacements(mDisplacements);

            // on charge le graph
            GraphType graph;
            _loadEdgeGraph(graph);

            // On calcul les déplacements
            std::set<std::string> sCollapsedEdges;
            std::vector<edge_descriptor> vDeformedEdges;
            _applyEdgeDisplacement(graph, mDisplacements, vDeformedEdges, sCollapsedEdges);

            // on enregistre les modifications
            _persistEdgeDisplacement(graph, vDeformedEdges);
        };

        ///
        ///
        ///
        void CFeatConnectionOp::_computeCpDisplacements(std::map<ign::geometry::Point, ign::math::Vec2d> & mDisplacements) const 
        {
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const &epgParams = context->getEpgParameters();

            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            // app params
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();

            std::string const edgeLinkName = themeParameters->getValue(EDGE_LINK).toString();
            std::string const landCoverTypeName = themeParameters->getValue(LAND_COVER_TYPE).toString();
            std::string const landAreaValue = themeParameters->getValue(TYPE_LAND_AREA).toString();

            ign::geometry::MultiPolygon mpLandmask;
            ign::feature::FeatureIteratorPtr itLandmask = _fsLandmask->getFeatures(ign::feature::FeatureFilter(landCoverTypeName + " = '" + landAreaValue + "' AND " + countryCodeName + " = '" + _countryCode + "'"));
            while (itLandmask->hasNext())
            {
                ign::feature::Feature const &fLandmask = itLandmask->next();
                ign::geometry::MultiPolygon const &mp = fLandmask.getGeometry().asMultiPolygon();
                for (int i = 0; i < mp.numGeometries(); ++i)
                {
                    mpLandmask.addGeometry(mp.polygonN(i));
                }
            }

            ign::feature::FeatureFilter filterCp(countryCodeName + " LIKE '%" + _countryCode + "%'");
            std::set<std::string> sEdge2Delete;

            // patience
            int numFeatures = epg::sql::tools::numFeatures(*_fsCp, filterCp);
            boost::progress_display display(numFeatures, std::cout, "[ cp_connection  % complete ]\n");

            ign::feature::FeatureIteratorPtr itCp = _fsCp->getFeatures(filterCp);
            while (itCp->hasNext())
            {
                ign::feature::Feature const &fCp = itCp->next();
                ign::geometry::Point const &cpGeom = fCp.getGeometry().asPoint();
                std::string edgeLink = fCp.getAttribute(edgeLinkName).toString();

                // if (fCp.getId() == "CONNECTINGPOINT48" || fCp.getId() == "CONNECTINGPOINT49") {
                //     bool test = true;
                // } else {
                //     continue;
                // }

                if (_verbose)
                    _logger->log(epg::log::DEBUG, fCp.getId());

                std::pair<bool, std::string> foundEdge = _getNearestEdge(cpGeom, countryCodeName, edgeLinkName, edgeLink, sEdge2Delete);

                if (!foundEdge.first)
                {
                    _logger->log(epg::log::ERROR, "No candidate edge found [" + edgeLinkName + "] " + edgeLink);
                    continue;
                }

                ign::feature::Feature fEdge;
                _fsEdge->getFeatureById(foundEdge.second, fEdge);

                if (_verbose)
                    _logger->log(epg::log::DEBUG, foundEdge.second);

                // if (foundEdge.second == "41f0d5bf-d002-4c8e-a828-779aeed290ae") {
                //     bool test = true;
                // } else {
                //     continue;
                // }

                ign::geometry::LineString const &edgeGeom = fEdge.getGeometry().asLineString();

                // On recupere la (les) partie(s) qui touche(nt) le pays
                ign::geometry::GeometryPtr intersectionGeom(mpLandmask.Intersection(edgeGeom));
                ign::geometry::LineString intersectedEdgeGeom;
                size_t nbResultingEdges = 0;

                if (intersectionGeom->isLineString())
                {
                    intersectedEdgeGeom = intersectionGeom->asLineString();
                }
                else if (intersectionGeom->isGeometryCollection())
                {
                    double distanceMax = std::numeric_limits<double>::infinity();
                    ign::geometry::GeometryCollection const &collectGeom = intersectionGeom->asGeometryCollection();
                    for (size_t i = 0; i < collectGeom.numGeometries(); ++i)
                    {
                        if (collectGeom.geometryN(i).isLineString())
                        {
                            ++nbResultingEdges;
                            ign::geometry::LineString const &ls = collectGeom.geometryN(i).asLineString();

                            ign::feature::Feature feat = fCp;
                            ign::geometry::LineString lsNoZ = ls;
                            lsNoZ.clearZ();
                            feat.setGeometry(lsNoZ);
                            _shapeLogger->writeFeature("multiple_result", feat);

                            double distance = ls.distance(cpGeom);
                            if (distance < distanceMax)
                            {
                                distanceMax = distance;
                                intersectedEdgeGeom = ls;
                            }
                        }
                    }
                }

                if (intersectedEdgeGeom.isEmpty())
                {
                    _logger->log(epg::log::ERROR, "The edge doesn't intersect the landmask [id] " + foundEdge.second);
                    continue;
                }

                bool bCreateNewFeature = false;
                if (nbResultingEdges > 1)
                {
                    sEdge2Delete.insert(foundEdge.second);
                    bCreateNewFeature = true;
                }

                app::geometry::tools::LineStringSplitter lsSplitter(intersectedEdgeGeom, 1e-5);

                ign::geometry::Point projectedPoint;
                epg::tools::geometry::projectZ(intersectedEdgeGeom, cpGeom, projectedPoint);

                ign::feature::Feature feat = fEdge;
                ign::geometry::Point projectedPointNoZ = projectedPoint;
                projectedPointNoZ.clearZ();
                feat.setGeometry(projectedPointNoZ);
                _shapeLogger->writeFeature("projected_cp", feat);

                lsSplitter.addCuttingGeometry(projectedPoint);
                ign::geometry::LineString newEdgeGeom = lsSplitter.truncAtEnds();

                newEdgeGeom.clearZ();

                ign::geometry::Point displacementStart = newEdgeGeom.startPoint().distance(cpGeom) < newEdgeGeom.endPoint().distance(cpGeom) ? newEdgeGeom.startPoint() : newEdgeGeom.endPoint();
                mDisplacements.insert(std::make_pair(displacementStart, cpGeom.toVec2d() - displacementStart.toVec2d()));

                fEdge.setGeometry(newEdgeGeom);
                if (bCreateNewFeature)
                {
                    fEdge.setId("");
                    _fsEdge->createFeature(fEdge);
                    _shapeLogger->writeFeature("created_edges", fEdge);
                }
                else
                {
                    _fsEdge->modifyFeature(fEdge);
                }

                ign::feature::Feature featLog;
                featLog.setGeometry(newEdgeGeom);
                _shapeLogger->writeFeature("resulting_edges", featLog);

                ++display;
            }

            std::set< std::string >::const_iterator sit;
            for( sit = sEdge2Delete.begin() ; sit != sEdge2Delete.end() ; ++sit ) {
                _fsEdge->deleteFeature(*sit);
            }

            std::map<ign::geometry::Point, ign::math::Vec2d>::const_iterator mit;
            for (mit = mDisplacements.begin(); mit != mDisplacements.end(); ++mit)
            {
                ign::feature::Feature feat;
                feat.setGeometry(ign::geometry::LineString(mit->first, ign::geometry::Point(mit->first.x() + mit->second.x(), mit->first.y() + mit->second.y())));
                _shapeLogger->writeFeature("cp_displacements", feat);
            }

        }

        ///
        ///
        ///
        void CFeatConnectionOp::_computeCp()
        {
            std::map<ign::geometry::Point, ign::math::Vec2d> mDisplacements;
            _computeCpDisplacements(mDisplacements);

            // on charge le graph
            GraphType graph;
            _loadEdgeGraph(graph);

            // On calcul les déplacements
            std::set<std::string> sCollapsedEdges;
            std::vector<edge_descriptor> vDeformedEdges;
            _applyEdgeDisplacement(graph, mDisplacements, vDeformedEdges, sCollapsedEdges);

            // on enregistre les modifications
            _persistEdgeDisplacement(graph, vDeformedEdges);
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_computeCpCl()
        {
            std::map<ign::geometry::Point, ign::math::Vec2d> mDisplacements;
            _computeCpDisplacements(mDisplacements);
            _computeClDisplacements(mDisplacements);

            // on charge le graph
            GraphType graph;
            _loadEdgeGraph(graph);

            // On calcul les déplacements
            std::set<std::string> sCollapsedEdges;
            std::vector<edge_descriptor> vDeformedEdges;
            _applyEdgeDisplacement(graph, mDisplacements, vDeformedEdges, sCollapsedEdges);

            // on enregistre les modifications
            _persistEdgeDisplacement(graph, vDeformedEdges);
        }

        ///
        ///
        ///
        std::pair<bool, std::string> CFeatConnectionOp::_getNearestEdge(
            ign::geometry::Geometry const& refGeom,
            std::string countryCodeName,
            std::string edgeLinkName,
            std::string edgeLink,
            std::set<std::string> const& sEdge2Delete ) const
        {

            ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(countryCodeName + " LIKE '%" + _countryCode + "%' AND " + edgeLinkName + " = '" + edgeLink + "'");

            double dMax = std::numeric_limits<double>::infinity();
            std::string idMax = "";
            while (itEdge->hasNext())
            {
                ign::feature::Feature const& fEdge = itEdge->next();
                if ( sEdge2Delete.find(fEdge.getId()) != sEdge2Delete.end()  ) continue;

                ign::geometry::LineString const &edgeGeom = fEdge.getGeometry().asLineString();
                double distance = edgeGeom.distance(refGeom);

                if (distance < dMax)
                {
                    dMax = distance;
                    idMax = fEdge.getId();
                }
            }

            bool found = (dMax != std::numeric_limits<double>::infinity());

            return std::make_pair(found, idMax);
        }

        ///
        ///
        ///
        std::pair<bool, std::string> CFeatConnectionOp::_getCountryEdgeLink(
            std::string edgeLinks,
            std::string countryCodes,
            std::string country) const
        {
            std::vector<std::string> vEdgesLink;
            ign::tools::StringManip::Split(edgeLinks, "#", vEdgesLink);

            std::vector<std::string> vCountryCode;
            ign::tools::StringManip::Split(countryCodes, "#", vCountryCode);

            int index = -1;
            for (size_t i = 0; i < vCountryCode.size(); ++i)
            {
                if (vCountryCode[i] == country)
                {
                    index = i;
                    break;
                }
            }

            if (index < 0 || index > vEdgesLink.size() - 1)
                return std::make_pair(false, "");

            return std::make_pair(true, vEdgesLink[index]);
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_loadEdgeGraph(GraphType & graph) const {
            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const &epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            

            ign::geometry::graph::builder::SimpleGraphBuilder<GraphType> graphBuilder(graph, 1e-5);

            ign::feature::FeatureFilter filterEdge(countryCodeName + " LIKE '%" + _countryCode + "%'");
            // ign::feature::FeatureFilter filterEdge(countryCodeName + " LIKE '%" + _countryCode + "%' AND inspireid='02d8e754-05bc-44bf-96c1-3bcc75ff740e'");
            ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(filterEdge);

            // patience
            int numFeatures2 = epg::sql::tools::numFeatures(*_fsEdge, filterEdge);
            boost::progress_display display2(numFeatures2, std::cout, "[ graph loading  % complete ]\n");

            while (itEdge->hasNext())
            {
                ign::feature::Feature const &fEdge = itEdge->next();
                ign::geometry::LineString const &edgeGeom = fEdge.getGeometry().asLineString();

                graphBuilder.addEdge(edgeGeom, fEdge.getId());

                ++display2;
            }
        }

        ///
        ///
        ///
        std::pair<bool, app::calcul::CFeatConnectionOp::vertex_descriptor> CFeatConnectionOp::_getNearestVertex(
            GraphType const &graph,
            ign::geometry::Point const &pt,
            double searchDistance) const
        {

            std::vector<vertex_descriptor> vVertices;
            graph.verticesIntersectingBox(pt.getEnvelope().expandBy(searchDistance), vVertices);

            vertex_descriptor v = GraphType::nullVertex();

            double minDistance = std::numeric_limits<double>::max();
            for (size_t i = 0; i < vVertices.size(); ++i)
            {
                double distance = pt.distance(graph.getGeometry(vVertices[i]));
                if (distance < minDistance)
                {
                    minDistance = distance;
                    v = vVertices[i];
                }
            }
            bool found = (v != GraphType::nullVertex());
            return std::make_pair(found, v);
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_applyEdgeDisplacement(
            GraphType & graph,
            std::map<ign::geometry::Point, ign::math::Vec2d> const & mReferences,
            std::vector<edge_descriptor> & vDeformedEdges,
            std::set<std::string> & sCollapsedEdges
        ) const {
            double influenceDist = 2;
            // double mergingDist = 1e-5;
            epg::calcul::matching::detail::LineStringSimpleDampedDeformer deformer;

            _applyDisplacement(graph, mReferences, deformer, vDeformedEdges, sCollapsedEdges, influenceDist /*, mergingDist*/);
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_applyDisplacement(
            GraphType & graph,
            std::map<ign::geometry::Point, ign::math::Vec2d> const & mReferences,
            epg::calcul::matching::detail::LineStringDeformer const & lineStringDeformer,
            std::vector<edge_descriptor> & vDeformedEdges,
            std::set<std::string> & sCollapsedEdges,
            double influenceDist /*,
             double mergingDist*/
        ) const
        {
            epg::log::ShapeLogger *shapeLogger = epg::log::ShapeLoggerS::getInstance();
            epg::Context *context = epg::ContextS::getInstance();

            if (_verbose)
                shapeLogger->addShape("applyDisplacements_deformedEdges", epg::log::ShapeLogger::LINESTRING);
            if (_verbose)
                shapeLogger->addShape("applyDisplacements_edgesBeforeDeformation", epg::log::ShapeLogger::LINESTRING);


            typedef typename GraphType::vertex_descriptor vertex_descriptor;
            typedef typename GraphType::edge_descriptor edge_descriptor;
            typedef typename GraphType::edge_iterator edge_iterator;

            // on construit la map des deplacements
            std::map<vertex_descriptor, std::vector<ign::math::Vec2d>> mDisplacements;

            // et on enregistre les vertices de reference
            std::map<vertex_descriptor, ign::math::Vec2d> mVertexRefs;

            // deplacement et enregistrement des references
            std::map<ign::geometry::Point, ign::math::Vec2d>::const_iterator rit;
            for (rit = mReferences.begin(); rit != mReferences.end(); ++rit)
            {
                std::pair<bool, vertex_descriptor> foundVertex = _getNearestVertex(graph, rit->first, 1e-5);
                if (!foundVertex.first)
                {
                    _logger->log(epg::log::WARN, "No node found at location : " + rit->first.toString());
                    continue;
                }

                typename std::map<vertex_descriptor, std::vector<ign::math::Vec2d>>::iterator
                    rit_bis = mDisplacements.insert(std::make_pair(foundVertex.second, std::vector<ign::math::Vec2d>())).first;

                rit_bis->second.push_back(rit->second);
                mVertexRefs.insert(std::make_pair(foundVertex.second, rit->second));
            }

            // patience
            boost::progress_display display1(mReferences.size(), std::cout, "[ applyDisplacement (1/2) % complete ]\n");

            std::map<vertex_descriptor, ign::math::Vec2d>::const_iterator mit;
            for (mit = mVertexRefs.begin(); mit != mVertexRefs.end(); ++mit, ++display1)
            {
                ign::math::Vec2d const &refVector = mit->second;

                ign::geometry::Point refPoint = graph.getGeometry(mit->first);
                double affectDist = refVector.norm() * influenceDist;

                std::set<vertex_descriptor> sVertices;
                graph.getVertices(refPoint.getEnvelope().expandBy(affectDist), sVertices);

                typename std::set<vertex_descriptor>::const_iterator vit;
                for (vit = sVertices.begin(); vit != sVertices.end(); ++vit)
                {
                    if (mVertexRefs.find(*vit) != mVertexRefs.end())
                        continue; // si c est un point de reference on passe

                    ign::geometry::Point vertexGeom = graph.getGeometry(*vit);

                    double distance = vertexGeom.distance2d(refPoint);
                    if (distance >= affectDist)
                        continue;

                    ign::math::Vec2d displacement = ((affectDist - distance) / affectDist) * refVector;

                    typename std::map<vertex_descriptor, std::vector<ign::math::Vec2d>>::iterator mit = mDisplacements.find(*vit);

                    if (mit == mDisplacements.end())
                        mit = mDisplacements.insert(std::make_pair(*vit, std::vector<ign::math::Vec2d>())).first;

                    mit->second.push_back(displacement);
                }
            }

            // on applique les deplacements
            // sur les arcs
            int numDeformedEdges = 0;
            std::vector<edge_descriptor> edgesToRemove;

            // patience
            boost::progress_display display2(graph.numEdges(), std::cout, "[ applyDisplacement (2/2) % complete ]\n");

            edge_iterator eit, eend;
            for (graph.edges(eit, eend); eit != eend; ++eit)
            {
                ++display2;

                typename std::map<vertex_descriptor, std::vector<ign::math::Vec2d>>::iterator dit_source, dit_target;

                dit_source = mDisplacements.find(graph.source(*eit));
                dit_target = mDisplacements.find(graph.target(*eit));

                if (dit_source == mDisplacements.end() && dit_target == mDisplacements.end())
                    continue;

                ign::geometry::LineString ls = graph.getGeometry(*eit);

                ign::math::Vec2d vectSource = (dit_source != mDisplacements.end()) ? _computeDisplacement(dit_source->second) : ign::math::Vec2d();

                ign::math::Vec2d vectTarget = (dit_target != mDisplacements.end()) ? _computeDisplacement(dit_target->second) : ign::math::Vec2d();

                // deformation
                lineStringDeformer.deform(vectSource, vectTarget, ls);
                // DEBUG
                // double distance = ls.startPoint().distance2d( ls.endPoint());
                /// \todo en cas de geometrie effondree assurer la continuite topologique
                /// supprimer le sommet deplace, mettre a jour les arcs lies au sommet deplacer (geometrie et sommet ini ou fin)
                if ((ls.startPoint().distance2d(ls.endPoint()) < 1e-5) && (graph.source(*eit) != graph.target(*eit)))
                {
                    edgesToRemove.push_back(*eit);
                    sCollapsedEdges.insert(graph.origins(*eit)[0]);
                    _logger->log(epg::log::INFO, "resulting geometry is collapsed for the edge [id] : " + graph.origins(*eit)[0]);
                }
                else
                {

                    if (_verbose)
                    {
                        ign::feature::Feature feat;
                        ign::geometry::LineString lsNoZ = graph.getGeometry(*eit);
                        lsNoZ.clearZ();
                        feat.setGeometry(lsNoZ);
                        shapeLogger->writeFeature("applyDisplacements_edgesBeforeDeformation", feat);
                    }

                    // 1) patch pour by-passer les IGN_ASSERT de setGeometry(edge)
                    ign::geometry::Point oldSourceGeom = graph.getGeometry(graph.source(*eit));
                    ign::geometry::Point oldTargetGeom = graph.getGeometry(graph.target(*eit));
                    graph.setGeometry(graph.source(*eit), ls.startPoint());
                    graph.setGeometry(graph.target(*eit), ls.endPoint());

                    graph.setGeometry(*eit, ls);

                    // 2) patch pour by-passer les IGN_ASSERT de setGeometry(edge)
                    graph.setGeometry(graph.source(*eit), oldSourceGeom);
                    graph.setGeometry(graph.target(*eit), oldTargetGeom);

                    if (_verbose)
                    {
                        ign::feature::Feature feat;
                        ls.clearZ();
                        feat.setGeometry(ls);
                        shapeLogger->writeFeature("applyDisplacements_deformedEdges", feat);
                    }

                    vDeformedEdges.push_back(*eit);
                    ++numDeformedEdges;
                }
            }

            for (size_t i = 0; i < edgesToRemove.size(); ++i)
                graph.removeEdge(edgesToRemove[i]);

            // sur les sommets
            typename std::map<vertex_descriptor, std::vector<ign::math::Vec2d>>::const_iterator dit;
            for (dit = mDisplacements.begin(); dit != mDisplacements.end(); ++dit)
            {
                ign::math::Vec2d displacement = _computeDisplacement(dit->second);
                // epg::graph::tools::translateVertex( graph, dit->first, displacement, false/*with merging*/ );
                ign::geometry::Point oldGeom = graph.getGeometry(dit->first);
                graph.setGeometry(dit->first, ign::geometry::Point(oldGeom.x() + displacement.x(), oldGeom.y() + displacement.y(), oldGeom.z()));
            }

            if (_verbose)
                shapeLogger->closeShape("applyDisplacements_deformedEdges");
            if (_verbose)
                shapeLogger->closeShape("applyDisplacements_edgesBeforeDeformation");
            _logger->log(epg::log::INFO, "nombre d arcs deformes :" + ign::data::Integer(numDeformedEdges).toString());
        }

        ///
        ///
        ///
        ign::math::Vec2d CFeatConnectionOp::_computeDisplacement(std::vector<ign::math::Vec2d> const &vVectors) const
        {
            ign::math::Vec2d displacement;

            for (size_t i = 0; i < vVectors.size(); ++i)
                displacement += vVectors[i];

            return (displacement / vVectors.size());
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_persistEdgeDisplacement(
            GraphType & graph,
            std::vector<edge_descriptor> & vDeformedEdges
        ) const
        {
            std::vector<edge_descriptor>::const_iterator vit;
            for ( vit = vDeformedEdges.begin() ; vit != vDeformedEdges.end() ; ++vit ) {
                std::string edgeId = graph.origins(*vit)[0];
                ign::geometry::LineString edgeGeom = graph.getGeometry(*vit);
                
                ign::feature::Feature fEdge;
                _fsEdge->getFeatureById(edgeId, fEdge);

                fEdge.setGeometry(edgeGeom);

                _fsEdge->modifyFeature(fEdge);
            }

        }
    }
}