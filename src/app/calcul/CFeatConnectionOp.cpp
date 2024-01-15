// APP
#include <app/calcul/CFeatConnectionOp.h>
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LineStringSplitter.h>
#include <app/tools/StringTools.h>
#include <app/tools/translateVertex.h>
#include <app/calcul/detail/ClMerger.h>

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
#include <ign/geometry/algorithm/LineMergerOpGeos.h>
#include <ign/geometry/algorithm/HausdorffDistanceOp.h>

namespace app
{
    namespace calcul
    {

        ///
        ///
        ///
        CFeatConnectionOp::CFeatConnectionOp(
            std::string countryCode,
            bool verbose
        ) : _countryCode(countryCode),
            _verbose(verbose)
        {
            _init();
        }

        ///
        ///
        ///
        CFeatConnectionOp::~CFeatConnectionOp()
        {
        }

        ///
        ///
        ///
        void CFeatConnectionOp::computeCp()
        {
            std::vector<std::string> vCountriesCodeName;
		    epg::tools::StringTools::Split(_countryCode, "#", vCountriesCodeName);

			for (std::vector<std::string>::iterator vit = vCountriesCodeName.begin(); vit != vCountriesCodeName.end(); ++vit) {
                _computeCp(*vit);
            }
        }

        ///
        ///
        ///
        void CFeatConnectionOp::computeCl()
        {
            std::vector<std::string> vCountriesCodeName;
		    epg::tools::StringTools::Split(_countryCode, "#", vCountriesCodeName);
            
            for (std::vector<std::string>::iterator vit = vCountriesCodeName.begin(); vit != vCountriesCodeName.end(); ++vit) {
                _computeCl(*vit);
            }
        }

        ///
        ///
        ///
        void CFeatConnectionOp::computeCpCl()
        {
            std::vector<std::string> vCountriesCodeName;
		    epg::tools::StringTools::Split(_countryCode, "#", vCountriesCodeName);
            
            for (std::vector<std::string>::iterator vit = vCountriesCodeName.begin(); vit != vCountriesCodeName.end(); ++vit) {
                _computeCpCl(*vit);
            }
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_init()
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
            std::string const edgeTableName = epgParams.getValue(EDGE_TABLE).toString();

            // app parameters
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            std::string const landmaskTableName = themeParameters->getValue(LANDMASK_TABLE).toString();
            std::string clTableName = themeParameters->getValue(CL_TABLE).toString();
            if ( clTableName == "" ) {
                std::string const clTableSuffix = themeParameters->getValue(CL_TABLE_SUFFIX).toString();
                clTableName = edgeTableName + clTableSuffix;
            }
		    std::string cpTableName = themeParameters->getValue(CP_TABLE).toString();
            if ( cpTableName == "" ) {
                std::string const cpTableSuffix = themeParameters->getValue(CP_TABLE_SUFFIX).toString();
                cpTableName = edgeTableName + cpTableSuffix;
            }

            //--
            _fsLandmask = context->getDataBaseManager().getFeatureStore(landmaskTableName, idName, geomName);

            //--
            if (cpTableName != "")
            {
                _fsCp = context->getDataBaseManager().getFeatureStore(cpTableName, idName, geomName);
            }

            //--
            if (clTableName != "")
            {
                _fsCl = context->getDataBaseManager().getFeatureStore(clTableName, idName, geomName);
            }

            //--
            _fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);

            //--
            _logger->log(epg::log::INFO, "[END] initialization: " + epg::tools::TimeTools::getTime());
        };

        ///
        ///
        ///
        void CFeatConnectionOp::_computeClDisplacements(std::map<ign::geometry::Point, ign::math::Vec2d> & mDisplacements, std::string const& country) const
        {
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const &epgParams = context->getEpgParameters();

            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            std::string const linkedFeatureIdName = epgParams.getValue(LINKED_FEATURE_ID).toString();

            // app params
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();

            std::string const edgeLinkName = themeParameters->getValue(EDGE_LINK).toString();
            std::string const landCoverTypeName = themeParameters->getValue(LAND_COVER_TYPE).toString();
            std::string const landAreaValue = themeParameters->getValue(TYPE_LAND_AREA).toString();
            double snapDistance = themeParameters->getValue(SNAP_DIST).toDouble();

            ign::feature::FeatureFilter filterCl(countryCodeName + " LIKE '%" + country + "%'");

            // patience
            int numFeatures = epg::sql::tools::numFeatures(*_fsCl, filterCl);
            boost::progress_display display(numFeatures, std::cout, "[ cl_connection  % complete ]\n");

            // pour garder le lien entre les CL et les edges nouvellement créés
            bimap_t mParentChilds;

            std::set<std::string> sTreatedCl;
            ign::feature::FeatureIteratorPtr itCl = _fsCl->getFeatures(filterCl);
            while (itCl->hasNext())
            {
                ++display;
                ign::feature::Feature const& fCl = itCl->next();
                ign::geometry::LineString const& clGeom = fCl.getGeometry().asLineString();
                std::string const linkedFeatureId = fCl.getAttribute(linkedFeatureIdName).toString();
                std::string const countryCode = fCl.getAttribute(countryCodeName).toString();

                // if (_verbose) _logger->log(epg::log::DEBUG, fCl.getId());

                if ( sTreatedCl.find(fCl.getId()) != sTreatedCl.end() ) continue;

                //DEBUG
                // std::string idDebug = fCl.getId();
                // if ( idDebug != "CONNECTINGLINE2096" && idDebug != "CONNECTINGLINE2093" && idDebug != "CONNECTINGLINE2095" && idDebug != "CONNECTINGLINE2094" && idDebug != "CONNECTINGLINE2086" ) continue;
                // if ( idDebug == "CONNECTINGLINE1817" || idDebug == "CONNECTINGLINE1816") {
                //     bool test = true;
                // }
                // if ( idDebug != "CONNECTINGLINE2218" && idDebug != "CONNECTINGLINE2230") continue;

                std::pair<bool, std::string> foundFeatureId = _getSingleValue(linkedFeatureId, countryCode, country);
                if (!foundFeatureId.first)
                {
                    _logger->log(epg::log::ERROR, "Feature id not found [connecting edge id] " + fCl.getId());
                    continue;
                }

                //DEBUG
                // std::map<std::string, std::vector<std::string>>::const_iterator mit = mParentChilds.find("fd8c7927-7e77-422e-a50e-67a4989550fc");
                // if (mit !=  mParentChilds.end()) {
                //     std::vector<std::string>>::const_iterator vit;
                //     for( vit = mit->second.begin() ; vit != mit->second.end() ; ++vit) {
                //         if (vit == foundFeatureId.second)
                //     }
                // }
                // if (foundFeatureId.second == "fd8c7927-7e77-422e-a50e-67a4989550fc") {
                //     bool test =true;
                // }

                //fusionner les cl adjacentes avec le même edgeLink, récuperer le géométrie fusionnée et lister les cl traitées pour ne pas les traiter de nouveau
                ign::geometry::LineString mergedClGeom = detail::ClMerger::merge(_fsCl, fCl, foundFeatureId.second, sTreatedCl);
                
                std::pair<bool, ign::feature::Feature> foundEdge = _getNearestChild(mergedClGeom, foundFeatureId.second, mParentChilds);
                std::string edgeId = foundEdge.second.getId();
                if (!foundEdge.first)
                {
                    _logger->log(epg::log::ERROR, "No candidate edge found [" + linkedFeatureIdName + "] " + foundFeatureId.second);
                    continue;
                }

                // if (_verbose) _logger->log(epg::log::DEBUG, "  "+edgeId);

                ign::geometry::LineString const& edgeGeom = foundEdge.second.getGeometry().asLineString();

                // peut-on connecter l'edge aux extremités du connecting edge ?
                std::vector<double> vDist;
                vDist.push_back(mergedClGeom.startPoint().distance(edgeGeom.startPoint()));
                vDist.push_back(mergedClGeom.startPoint().distance(edgeGeom.endPoint()));
                vDist.push_back(mergedClGeom.endPoint().distance(edgeGeom.startPoint()));
                vDist.push_back(mergedClGeom.endPoint().distance(edgeGeom.endPoint()));

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
                        epg::tools::geometry::projectZ(edgeGeom, mergedClGeom.startPoint(), startPoint);
                        if (startPoint.distance((minId == 0) ? edgeGeom.startPoint() : edgeGeom.endPoint()) <= snapDistance) {
                            startPoint = (minId == 0) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                            endingPoint = startPoint;
                        } else {
                            mpCuttingPoints.addGeometry(startPoint);
                        }
                    }
                    ign::math::Vec2d v1 = mergedClGeom.startPoint().toVec2d() - startPoint.toVec2d();
                    mDisplacements.insert(std::make_pair(startPoint, v1));
                    double n1 = v1.norm();
                    if (n1 > 10) {
                        _logger->log(epg::log::WARN, "Big displacement : "+std::to_string(n1)+" [" + linkedFeatureIdName + "] " + foundFeatureId.second);
                    }

                    // size_t minId2 = vDist[2] < vDist[3] ? 2 : 3;
                    size_t minId2 = minId == 0  ? 3 : 2;
                    if (vDist[minId2] <= snapDistance)
                    {
                        startPoint = (minId2 == 2) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                        endingPoint = startPoint;
                    }
                    else
                    {
                        epg::tools::geometry::projectZ(edgeGeom, mergedClGeom.endPoint(), startPoint);
                        if (startPoint.distance((minId2 == 2) ? edgeGeom.startPoint() : edgeGeom.endPoint()) <= snapDistance) {
                            startPoint = (minId2 == 2) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                            endingPoint = startPoint;
                        } else {
                            mpCuttingPoints.addGeometry(startPoint);
                        }
                    }
                    ign::math::Vec2d v2 = mergedClGeom.endPoint().toVec2d() - startPoint.toVec2d();
                    mDisplacements.insert(std::make_pair(startPoint, v2));
                    double n2 = v2.norm();
                    if (n2 > 10) {
                        _logger->log(epg::log::WARN, "Big displacement : "+std::to_string(n2)+" [" + linkedFeatureIdName + "] " + foundFeatureId.second);
                    }
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
                        epg::tools::geometry::projectZ(edgeGeom, mergedClGeom.endPoint(), startPoint);
                        if (startPoint.distance((minId == 2) ? edgeGeom.startPoint() : edgeGeom.endPoint()) <= snapDistance) {
                            startPoint = (minId == 2) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                            endingPoint = startPoint;
                        } else {
                            mpCuttingPoints.addGeometry(startPoint);
                        }
                    }
                    ign::math::Vec2d v1 = mergedClGeom.endPoint().toVec2d() - startPoint.toVec2d();
                    mDisplacements.insert(std::make_pair(startPoint, v1));
                    double n1 = v1.norm();
                    if (n1 > 10) {
                        _logger->log(epg::log::WARN, "Big displacement : "+std::to_string(n1)+" [" + linkedFeatureIdName + "] " + foundFeatureId.second);
                    }

                    // size_t minId2 = vDist[0] < vDist[1] ? 0 : 1;
                    size_t minId2 = minId == 2  ? 1 : 0;
                    if (vDist[minId2] <= snapDistance)
                    {
                        startPoint = (minId2 == 0) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                        endingPoint = startPoint;
                    }
                    else
                    {
                        epg::tools::geometry::projectZ(edgeGeom, mergedClGeom.startPoint(), startPoint);
                        if (startPoint.distance((minId2 == 0) ? edgeGeom.startPoint() : edgeGeom.endPoint()) <= snapDistance) {
                            startPoint = (minId2 == 0) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                            endingPoint = startPoint;
                        } else {
                            mpCuttingPoints.addGeometry(startPoint);
                        }
                    }
                    ign::math::Vec2d v2 = mergedClGeom.startPoint().toVec2d() - startPoint.toVec2d();
                    mDisplacements.insert(std::make_pair(startPoint, v2));
                    double n2 = v2.norm();
                    if (n2 > 10) {
                        _logger->log(epg::log::WARN, "Big displacement : "+std::to_string(n2)+" [" + linkedFeatureIdName + "] " + foundFeatureId.second);
                    }
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

                //on supprime l'edge
                _shapeLogger->writeFeature("cl_deleted_features_"+country, foundEdge.second);
                _fsEdge->deleteFeature(edgeId);

                auto r_mit = mParentChilds.right.find(edgeId);
                std::string parentId = r_mit != mParentChilds.right.end() ? r_mit->second : edgeId;
                if ( vNewGeom.size() > 0 && r_mit != mParentChilds.right.end()) mParentChilds.right.erase(r_mit);

                for (size_t i = 0 ; i < vNewGeom.size() ; ++i) {
                    foundEdge.second.setGeometry(vNewGeom[i]);
                    
                    _fsEdge->createFeature(foundEdge.second);
                    mParentChilds.insert(value_type(parentId, foundEdge.second.getId()));

                    _shapeLogger->writeFeature("cl_created_features_"+country, foundEdge.second);
                }
            }

            std::map<ign::geometry::Point, ign::math::Vec2d>::const_iterator mit;
            for (mit = mDisplacements.begin(); mit != mDisplacements.end(); ++mit)
            {
                ign::feature::Feature feat;
                ign::geometry::Point p = mit->first;
                feat.setGeometry(ign::geometry::LineString(p, ign::geometry::Point(mit->first.x() + mit->second.x(), mit->first.y() + mit->second.y())));
                _shapeLogger->writeFeature("cl_displacements_"+country, feat);
            }
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_computeCl(std::string const& country)
        {
            _shapeLogger = epg::log::ShapeLoggerS::getInstance();
            _shapeLogger->addShape("cl_displacements_"+country, epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("cl_created_features_"+country, epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("cl_deleted_features_"+country, epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("cl_superposed_edges_"+country, epg::log::ShapeLogger::LINESTRING);

            std::map<ign::geometry::Point, ign::math::Vec2d> mDisplacements;
            _computeClDisplacements(mDisplacements, country);

            // on charge le graph
            GraphType graph;
            _loadEdgeGraph(graph, country);

            // On calcul les déplacements
            std::set<std::string> sCollapsedEdges;
            std::vector<edge_descriptor> vDeformedEdges;
            _applyEdgeDisplacement(graph, mDisplacements, vDeformedEdges, sCollapsedEdges);


            //DEBUG
            // std::set<edge_descriptor> sVisitedEdge;
            // edge_iterator eit, eend;
            // for (graph.edges(eit, eend); eit != eend; ++eit)
            // {
            //     if ( sVisitedEdge.find(*eit) != sVisitedEdge.end() ) continue;
            //     std::vector< oriented_edge_descriptor > vParallelEdges;
            //     graph.edges( graph.source(*eit), graph.target(*eit), vParallelEdges );
            //     if(vParallelEdges.size() < 2) continue;

            //     std::vector< oriented_edge_descriptor >::const_iterator vit;
            //     std::vector< oriented_edge_descriptor >::const_iterator vit_last = --vParallelEdges.end();
            //     std::set<edge_descriptor> sVisited;
            //     for (vit = vParallelEdges.begin() ; vit != vit_last ; ++vit ) {
            //         if ( sVisitedEdge.find(vit->descriptor) != sVisitedEdge.end() ) continue;
            //         sVisitedEdge.insert(vit->descriptor);
            //         ign::geometry::LineString lsRef = graph.getGeometry(*vit);
            //         std::vector< oriented_edge_descriptor >::const_iterator vit2 = vit;
            //         for ( ++vit2 ; vit2 != vParallelEdges.end() ; ++vit2 ) {
            //             if( vit2->descriptor == vit->descriptor ) continue; /*gestion des boucle*/
            //             ign::geometry::LineString ls = graph.getGeometry(vit2->descriptor);
            //             if ( ign::geometry::algorithm::HausdorffDistanceOp::distance(lsRef, ls) < 0.1 /*todo a ajuster*/ ) {
            //                 sVisitedEdge.insert(vit2->descriptor);
            //                 ign::feature::Feature eFeat;
            //                 eFeat.setGeometry(ls);
            //                 _shapeLogger->writeFeature("cl_superposed_edges_"+country, eFeat);
            //             }
            //         }
            //     }
            // }

            // on enregistre les modifications
            _persistEdgeDisplacement(graph, vDeformedEdges);

            //remove collapsed edges
            for (std::set<std::string>::const_iterator sit = sCollapsedEdges.begin() ; sit != sCollapsedEdges.end() ; ++sit) {
                _fsEdge->deleteFeature(*sit);
            }

            _shapeLogger->closeShape("cl_displacements_"+country);
            _shapeLogger->closeShape("cl_created_features_"+country);
            _shapeLogger->closeShape("cl_deleted_features_"+country);
            _shapeLogger->closeShape("cl_superposed_edges_"+country);
        };

        ///
        ///
        ///
        void CFeatConnectionOp::_computeCpDisplacements(std::map<ign::geometry::Point, ign::math::Vec2d> & mDisplacements, std::string const& country) const 
        {
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const &epgParams = context->getEpgParameters();

            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            std::string const linkedFeatureIdName = epgParams.getValue(LINKED_FEATURE_ID).toString();

            // app params
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            
            // std::string const edgeLinkName = themeParameters->getValue(EDGE_LINK).toString();
            // std::string const landCoverTypeName = themeParameters->getValue(LAND_COVER_TYPE).toString();
            // std::string const landAreaValue = themeParameters->getValue(TYPE_LAND_AREA).toString();
            double const snapDistance = themeParameters->getValue(SNAP_DIST).toDouble();

            // ign::geometry::MultiPolygon mpLandmask;
            // ign::feature::FeatureIteratorPtr itLandmask = _fsLandmask->getFeatures(ign::feature::FeatureFilter(landCoverTypeName + " = '" + landAreaValue + "' AND " + countryCodeName + " = '" + _countryCode + "'"));
            // while (itLandmask->hasNext())
            // {
            //     ign::feature::Feature const &fLandmask = itLandmask->next();
            //     ign::geometry::MultiPolygon const &mp = fLandmask.getGeometry().asMultiPolygon();
            //     for (int i = 0; i < mp.numGeometries(); ++i)
            //     {
            //         mpLandmask.addGeometry(mp.polygonN(i));
            //     }
            // }

            ign::feature::FeatureFilter filterCp(countryCodeName + " LIKE '%" + country + "%'");

            // patience
            int numFeatures = epg::sql::tools::numFeatures(*_fsCp, filterCp);
            boost::progress_display display(numFeatures, std::cout, "[ cp_connection  % complete ]\n");

            // pour garder le lien entre les CP et les edges nouvellement créés
            bimap_t mParentChilds;

            ign::feature::FeatureIteratorPtr itCp = _fsCp->getFeatures(filterCp);
            while (itCp->hasNext())
            {
                ++display;
                ign::feature::Feature const& fCp = itCp->next();

                if (_verbose) _logger->log(epg::log::DEBUG, fCp.getId());

                ign::geometry::Point const& cpGeom = fCp.getGeometry().asPoint();
                std::string const linkedFeatureId = fCp.getAttribute(linkedFeatureIdName).toString();
                std::string const countryCode = fCp.getAttribute(countryCodeName).toString();

                // if (fCp.getId() == "CONNECTINGPOINT48" || fCp.getId() == "CONNECTINGPOINT49") {
                //     bool test = true;
                // } else {
                //     continue;
                // }

                std::pair<bool, std::string> foundFeatureId = _getSingleValue(linkedFeatureId, countryCode, country);
                if (!foundFeatureId.first)
                {
                    _logger->log(epg::log::ERROR, "Feature id not found [connecting point id] " + fCp.getId());
                    continue;
                }

                // ign::feature::Feature fEdge;
                // _fsEdge->getFeatureById(foundFeatureId.second, fEdge);

                // if (_verbose)
                //     _logger->log(epg::log::DEBUG, foundFeatureId.second);

                // if (foundEdge.second == "41f0d5bf-d002-4c8e-a828-779aeed290ae") {
                //     bool test = true;
                // } else {
                //     continue;
                // }

                // ign::geometry::LineString const& edgeGeom = fEdge.getGeometry().asLineString();

                std::pair<bool, ign::feature::Feature> foundEdge = _getNearestChild(cpGeom, foundFeatureId.second, mParentChilds);
                std::string edgeId = foundEdge.second.getId();
                if (!foundEdge.first)
                {
                    _logger->log(epg::log::ERROR, "No candidate edge found [" + linkedFeatureIdName + "] " + foundFeatureId.second);
                    continue;
                }

                // if (_verbose) _logger->log(epg::log::DEBUG, "  "+edgeId);

                ign::geometry::LineString const& edgeGeom = foundEdge.second.getGeometry().asLineString();


                ign::geometry::Point startPoint;
                // peut-on connecter l'edge aux extremités du connecting edge ?
                double dCpStart = edgeGeom.startPoint().distance(cpGeom);
                double dCpEnd = edgeGeom.endPoint().distance(cpGeom);
                if (std::min(dCpStart, dCpEnd) > snapDistance) {
                    epg::tools::geometry::projectZ(edgeGeom, cpGeom, startPoint);
                    dCpStart = edgeGeom.startPoint().distance(startPoint);
                    dCpEnd = edgeGeom.endPoint().distance(startPoint);
                }
                    
                if (std::min(dCpStart, dCpEnd) <= snapDistance) {
                    startPoint = dCpStart < dCpEnd ? edgeGeom.startPoint() : edgeGeom.endPoint();
                } else {
                    app::geometry::tools::LineStringSplitter lsSplitter(edgeGeom, 1e-5);
                    lsSplitter.addCuttingGeometry(startPoint);
                    std::vector<ign::geometry::LineString> vNewEdgeGeom = lsSplitter.getSubLineStringsZ();

                    if (vNewEdgeGeom.size() == 2) {
                        startPoint = vNewEdgeGeom.front().endPoint();
                    } else {
                        startPoint = edgeGeom.startPoint().distance(cpGeom) < edgeGeom.endPoint().distance(cpGeom) ? edgeGeom.startPoint() : edgeGeom.endPoint();
                    }

                    if (vNewEdgeGeom.size() > 1) {
                        _fsEdge->deleteFeature(edgeId);

                        auto r_mit = mParentChilds.right.find(edgeId);
                        std::string parentId = r_mit != mParentChilds.right.end() ? r_mit->second : edgeId;
                        if (r_mit != mParentChilds.right.end()) mParentChilds.right.erase(r_mit);
                        for (size_t i = 0 ; i < vNewEdgeGeom.size() ; ++i) {
                            foundEdge.second.setGeometry(vNewEdgeGeom[i]);
                            
                            _fsEdge->createFeature(foundEdge.second);
                            mParentChilds.insert(value_type(parentId, foundEdge.second.getId()));
                        }
                    }
                }
                mDisplacements.insert(std::make_pair(startPoint, cpGeom.toVec2d() - startPoint.toVec2d()));
            }

            std::map<ign::geometry::Point, ign::math::Vec2d>::const_iterator mit;
            for (mit = mDisplacements.begin(); mit != mDisplacements.end(); ++mit)
            {
                ign::feature::Feature feat;
                ign::geometry::Point p = mit->first;
                feat.setGeometry(ign::geometry::LineString(p, ign::geometry::Point(mit->first.x() + mit->second.x(), mit->first.y() + mit->second.y())));
                _shapeLogger->writeFeature("cp_displacements_"+country, feat);
            }
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_computeCp(std::string const& country)
        {
            _shapeLogger = epg::log::ShapeLoggerS::getInstance();
            _shapeLogger->addShape("cp_displacements_"+country, epg::log::ShapeLogger::LINESTRING);

            std::map<ign::geometry::Point, ign::math::Vec2d> mDisplacements;
            _computeCpDisplacements(mDisplacements, country);

            // on charge le graph
            GraphType graph;
            _loadEdgeGraph(graph, country);

            // On calcul les déplacements
            std::set<std::string> sCollapsedEdges;
            std::vector<edge_descriptor> vDeformedEdges;
            _applyEdgeDisplacement(graph, mDisplacements, vDeformedEdges, sCollapsedEdges);

            // on enregistre les modifications
            _persistEdgeDisplacement(graph, vDeformedEdges);

            //remove collapsed edges
            for (std::set<std::string>::const_iterator sit = sCollapsedEdges.begin() ; sit != sCollapsedEdges.end() ; ++sit) {
                _fsEdge->deleteFeature(*sit);
            }

            _shapeLogger->closeShape("cp_displacements_"+country);
        }

        ///
        ///
        ///
        void CFeatConnectionOp::_computeCpCl(std::string const& country)
        {
            std::map<ign::geometry::Point, ign::math::Vec2d> mDisplacements;
            _computeCpDisplacements(mDisplacements, country);
            _computeClDisplacements(mDisplacements, country);

            // on charge le graph
            GraphType graph;
            _loadEdgeGraph(graph, country);

            // On calcul les déplacements
            std::set<std::string> sCollapsedEdges;
            std::vector<edge_descriptor> vDeformedEdges;
            _applyEdgeDisplacement(graph, mDisplacements, vDeformedEdges, sCollapsedEdges);

            // on enregistre les modifications
            _persistEdgeDisplacement(graph, vDeformedEdges);

            //remove collapsed edges
            for (std::set<std::string>::const_iterator sit = sCollapsedEdges.begin() ; sit != sCollapsedEdges.end() ; ++sit) {
                _fsEdge->deleteFeature(*sit);
            }
        }

        ///
        ///
        ///
        std::pair<bool, ign::feature::Feature> CFeatConnectionOp::_getNearestChild(
            ign::geometry::Geometry const& refGeom,
            std::string const& parentFeatureId,
            bimap_t const& mParentChilds ) const
        {
            epg::Context* context = epg::ContextS::getInstance();
            epg::params::EpgParameters const &epgParams = context->getEpgParameters();
            std::string const edgeIdName = epgParams.getValue(ID).toString();

            ign::feature::Feature fEdgeResult;
            bool found = false;

            std::vector<std::string> vCandidates;

            // bimap_t::const_iterator mit = mParentChilds.find(parentFeatureId);
            // if ( mit == mParentChilds.end() ) vCandidates.push_back(parentFeatureId);
            // else vCandidates = mit->second;

            auto range = mParentChilds.left.equal_range(parentFeatureId);
            for (auto l_it = range.first; l_it != range.second; ++l_it) {
                vCandidates.push_back(l_it->second);
            }
            if (vCandidates.size() == 0) {
                vCandidates.push_back(parentFeatureId);
            }

            if (vCandidates.size() == 1) {
                _fsEdge->getFeatureById(vCandidates.front(), fEdgeResult);
                found = !fEdgeResult.getId().empty();
            } else {
                double dMax = std::numeric_limits<double>::infinity();
                ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(edgeIdName + " IN ('" + tools::StringTools::toString(vCandidates, "','") + "')");

                while (itEdge->hasNext())
                {
                    ign::feature::Feature const& fEdge = itEdge->next();

                    ign::geometry::LineString const &edgeGeom = fEdge.getGeometry().asLineString();
                    double distance = edgeGeom.distance(refGeom);

                    if (distance < dMax)
                    {
                        dMax = distance;
                        fEdgeResult = fEdge;
                    }
                }
                found = (dMax != std::numeric_limits<double>::infinity());
            }
            
            return std::make_pair(found, fEdgeResult);
        }

        ///
        ///
        ///
        std::pair<bool, std::string> CFeatConnectionOp::_getSingleValue(
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
        void CFeatConnectionOp::_loadEdgeGraph(GraphType & graph, std::string const& country) const {
            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const &epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            

            ign::geometry::graph::builder::SimpleGraphBuilder<GraphType> graphBuilder(graph, 1e-5);

            ign::feature::FeatureFilter filterEdge(countryCodeName + " LIKE '%" + country + "%'");
            ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(filterEdge);

            // patience
            int numFeatures2 = epg::sql::tools::numFeatures(*_fsEdge, filterEdge);
            boost::progress_display display2(numFeatures2, std::cout, "[ graph loading  % complete ]\n");

            while (itEdge->hasNext())
            {
                ign::feature::Feature const& fEdge = itEdge->next();
                ign::geometry::LineString const& edgeGeom = fEdge.getGeometry().asLineString();

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
            double influenceDist = 0;
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
            boost::progress_display display1(mVertexRefs.size(), std::cout, "[ applyDisplacement (1/2) % complete ]\n");

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
                        ign::geometry::LineString const& ls = graph.getGeometry(*eit);
                        feat.setGeometry(ls);
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
            std::map<edge_descriptor, edge_descriptor> mOldNewEdges;
            std::set<edge_descriptor> sEdges2remove;
		    std::set<vertex_descriptor> sVertices2remove;
            typename std::map<vertex_descriptor, std::vector<ign::math::Vec2d>>::const_iterator dit;
            for (dit = mDisplacements.begin(); dit != mDisplacements.end(); ++dit)
            {
                ign::math::Vec2d displacement = _computeDisplacement(dit->second);
                tools::translateVertex( graph, dit->first, displacement, mOldNewEdges, sEdges2remove, sVertices2remove, true/*with merging*/, 1e-5 );
                // ign::geometry::Point oldGeom = graph.getGeometry(dit->first);
                // graph.setGeometry(dit->first, ign::geometry::Point(oldGeom.x() + displacement.x(), oldGeom.y() + displacement.y(), oldGeom.z()));
            }
            for ( std::set<edge_descriptor>::const_iterator sit = sEdges2remove.begin() ; sit != sEdges2remove.end() ; ++sit )
                graph.removeEdge(*sit);
            for ( std::set<vertex_descriptor>::const_iterator sit = sVertices2remove.begin() ; sit != sVertices2remove.end() ; ++sit )
                graph.removeVertex(*sit);

            for ( std::vector<edge_descriptor>::iterator vit = vDeformedEdges.begin() ; vit != vDeformedEdges.end() ; ++vit) {
                std::map<edge_descriptor, edge_descriptor>::const_iterator mit = mOldNewEdges.find(*vit);
                if (mit != mOldNewEdges.end()) {
                    *vit = mit->second;

                    // on fait une 2eme passe car l edge peut avoir ete deforme a la source et a la target
                    // ce qui fait qu il a 2 parents dans la map mOldNewEdges
                    mit = mOldNewEdges.find(*vit);
                    if (mit != mOldNewEdges.end()) {
                        *vit = mit->second;
                    }
                }
                
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
                ign::geometry::LineString edgeGeom = graph.getGeometry(*vit);

                std::string edgeId = graph.origins(*vit)[0];

                // _logger->log(epg::log::DEBUG, edgeId);
                
                ign::feature::Feature fEdge;
                _fsEdge->getFeatureById(edgeId, fEdge);

                fEdge.setGeometry(edgeGeom);

                _fsEdge->modifyFeature(fEdge);
            }

        }

		///
		///
		///
		void CFeatConnectionOp::computeClImport()
		{

			epg::Context* context = epg::ContextS::getInstance();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

			ign::feature::FeatureIteratorPtr itCL = _fsCl->getFeatures(ign::feature::FeatureFilter(countryCodeName + " = '" + _countryCode + "'"));

			while (itCL->hasNext())
			{
				ign::feature::Feature fCL = itCL->next();
				ign::geometry::LineString lsCL = fCL.getGeometry().asLineString();

				ign::feature::Feature fNewEdge = fCL;
				lsCL.setFillZ(0);
				fNewEdge.setGeometry(lsCL);
				_fsEdge->createFeature(fNewEdge);
			}
		}
    }
}

