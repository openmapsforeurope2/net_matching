// APP
#include <app/calcul/EdgeConnectorOp.h>
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LineStringSplitter.h>
#include <app/tools/StringTools.h>
// #include <app/tools/translateVertex.h>


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
#include <app/calcul/detail/LineStringAbsDampedDeformer.h>
#include <epg/tools/FilterTools.h>
#include <epg/tools/geometry/LineIntersector.h>

//#include <ome2/calcul/detail/ClMerger.h>

// SOCLE
#include <ign/geometry/graph/builder/SimpleGraphBuilder.h>
#include <ign/tools/stringtools.h>
#include <ign/geometry/algorithm/HausdorffDistanceOp.h>
#include <ign/math/Line2T.h>
#include <ign/numeric/Numeric.h>

namespace app
{
    namespace calcul
    {

        ///
        ///
        ///
        void EdgeConnectorOp::compute(
            std::string borderCode,
            bool verbose)
        {
            EdgeConnectorOp EdgeConnectorOp(borderCode, verbose);
            EdgeConnectorOp._compute();
        }

        ///
        ///
        ///
        EdgeConnectorOp::EdgeConnectorOp(
            std::string borderCode,
            bool verbose
        ) : _borderCode(borderCode),
            _verbose(verbose)
        {
            _init();
        }

        ///
        ///
        ///
        EdgeConnectorOp::~EdgeConnectorOp()
        {
            _shapeLogger->closeShape("ec_projected_antennas");
            _shapeLogger->closeShape("ec_split_edges");
            _shapeLogger->closeShape("ec_trimed_edges_parts");
        }

        ///
        ///
        ///
        void EdgeConnectorOp::_init()
        {
            //--
            _logger = epg::log::EpgLoggerS::getInstance();
            _logger->log(epg::log::INFO, "[START] initialization: " + epg::tools::TimeTools::getTime());

            //--
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();

            std::string const edgeTableName = epgParams.getValue(EDGE_TABLE).toString();
            std::string const idName = epgParams.getValue(ID).toString();
            std::string const geomName = epgParams.getValue(GEOM).toString();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            std::string const boundaryTableName = epgParams.getValue(TARGET_BOUNDARY_TABLE).toString();

            // app parameters
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            std::string const landCoverTypeName = themeParameters->getValue(LAND_COVER_TYPE).toString();
            std::string const landAreaValue = themeParameters->getValue(TYPE_LAND_AREA).toString();
            std::string const landmaskTableName = themeParameters->getValue(LANDMASK_TABLE).toString();
            double const landmaskBuffer = themeParameters->getValue(EC_LANDMASK_BUFFER).toDouble();

            //--
            _fsLandmask = context->getDataBaseManager().getFeatureStore(landmaskTableName, idName, geomName);

            // on recupere un buffer autour de la frontiere
            ign::geometry::GeometryPtr boundBuffPtr(new ign::geometry::Polygon());
            ign::feature::sql::FeatureStorePostgis* fsBoundary = context->getDataBaseManager().getFeatureStore(boundaryTableName, idName, geomName);
            ign::feature::FeatureIteratorPtr itBoundary = fsBoundary->getFeatures(ign::feature::FeatureFilter(countryCodeName +" = '"+_borderCode+"'"));
            while (itBoundary->hasNext())
            {
                ign::feature::Feature const& fBoundary = itBoundary->next();
                ign::geometry::LineString const& ls = fBoundary.getGeometry().asLineString();

                ign::geometry::GeometryPtr tmpBuffPtr(ls.buffer(3000));

                boundBuffPtr.reset(boundBuffPtr->Union(*tmpBuffPtr));
            }

            //on recupere la geometry des pays
            std::vector<std::string> vCountry;
		    epg::tools::StringTools::Split(_borderCode, "#", vCountry);

            for (std::vector<std::string>::iterator vit = vCountry.begin() ; vit != vCountry.end() ; ++vit) {
                ign::geometry::MultiPolygon mpLandmask;
                ign::feature::sql::FeatureStorePostgis* fsLandmask = context->getDataBaseManager().getFeatureStore(landmaskTableName, idName, geomName);
                ign::feature::FeatureIteratorPtr itLandmask = fsLandmask->getFeatures(ign::feature::FeatureFilter(landCoverTypeName + " = '" + landAreaValue + "' AND " + countryCodeName + " = '" + *vit + "'"));
                while (itLandmask->hasNext())
                {
                    ign::feature::Feature const& fLandmask = itLandmask->next();
                    ign::geometry::MultiPolygon const& mp = fLandmask.getGeometry().asMultiPolygon();
                    for (int i = 0; i < mp.numGeometries(); ++i)
                    {
                        mpLandmask.addGeometry(mp.polygonN(i));
                    }
                }

                //on calcul la geometry de travail
                _mCountryGeomPtr.insert(std::make_pair(*vit, ign::geometry::GeometryPtr(boundBuffPtr->Intersection(mpLandmask)) ));
                // _mCountryGeomWithBuffPtr.insert(std::make_pair(*vit, ign::geometry::GeometryPtr(_mCountryGeomPtr[*vit]->buffer(-1*landmaskBuffer))));
            }

            //--
            _fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);

            //--
            _shapeLogger = epg::log::ShapeLoggerS::getInstance();
            _shapeLogger->addShape("ec_projected_antennas", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("ec_split_edges", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("ec_trimed_edges_parts", epg::log::ShapeLogger::LINESTRING);

            //--
            _logger->log(epg::log::INFO, "[END] initialization: " + epg::tools::TimeTools::getTime());
        };

        ///
        ///
        ///
        void EdgeConnectorOp::_compute()
        {
            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const geomName = epgParams.getValue(GEOM).toString();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            double const snapDist = themeParameters->getValue(EC_SNAP_DIST).toDouble();
            double const snap2EdgeEndingsDist = themeParameters->getValue(EC_SNAP_2_EDGE_END_DIST).toDouble();
            double const antennaMinLength = themeParameters->getValue( ECL_ANTENNA_MIN_LENGTH ).toDouble();
            double const antennaMinDist2Neighbor = themeParameters->getValue( ECL_ANTENNA_MIN_DIST_2_NEIGHBOR ).toDouble();
            double const landmaskBuffer = themeParameters->getValue(EC_LANDMASK_BUFFER).toDouble();

            std::vector<std::string> vCountry;
		    epg::tools::StringTools::Split(_borderCode, "#", vCountry);

            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager, false);
            GraphType const& graph = graphManager.getGraph();
            
            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> vpAntennas;
            std::set< vertex_descriptor > visitedVertices;

            boost::progress_display display(graph.numVertices(), std::cout, "[ search antennas  % complete ]\n");
            vertex_iterator vit, vend;
            for( graph.vertices( vit, vend ) ; vit != vend ; ++vit )
            {
                ++display;

                if( visitedVertices.find( *vit ) != visitedVertices.end() ) continue;
                if( graph.degree( *vit ) != 1 ) continue;

                // seulement les vertex qui touchent une CL ?

                std::list<oriented_edge_descriptor> lAntennaEdges;

                std::vector< oriented_edge_descriptor > vEdges;
                graph.incidentEdges( *vit, vEdges );

                oriented_edge_descriptor nextEdge = vEdges.front(); // si nextEdge n'est pas une CL ?
                if (graphManager.isCl(nextEdge.descriptor)) {
                    _logger->log(epg::log::WARN, "Antenna is connecting line [cl id] "+graph.origins(nextEdge.descriptor)[0]);
                    continue;
                }
                std::string country = graphManager.getCountry(nextEdge.descriptor);

                vertex_descriptor vTarget = GraphType::nullVertex();
                
                while( true )
                {
                    if (graphManager.getCountry(nextEdge.descriptor) != country) break;

                    lAntennaEdges.push_back(nextEdge);

                    vTarget = graph.target( nextEdge );

                    if (graphManager.isCl(nextEdge.descriptor) && graph.degree(graph.source( nextEdge )) == 2 ) {
                        _logger->log(epg::log::WARN, "Antenna connected to connecting line [cl id] "+graph.origins(nextEdge.descriptor)[0]);
                        break;
                    }

                    if( graph.degree( vTarget ) != 2 ) { // ou si nextEdge est une CL ?
                        // if( graph.degree( vTarget ) == 1 /*antenne isolee*/)
                        // {
                        //     visitedVertices.insert( vTarget );
                        // }
                        break;
                    }

                    std::vector< oriented_edge_descriptor > vIncEdges;
                    graph.incidentEdges( vTarget, vIncEdges );

                    nextEdge = ( vIncEdges.front().descriptor == nextEdge.descriptor )? vIncEdges.back():vIncEdges.front();
                }

                if (!lAntennaEdges.empty()) vpAntennas.push_back(std::make_pair(country, lAntennaEdges));
            }

            double absThreshold = 30;
            double influenceFactor = 5;
            double snapDist2 = 1;
            epg::calcul::matching::detail::LineStringAbsDampedDeformer deformer(absThreshold, influenceFactor, snapDist2);

            boost::progress_display display2(vpAntennas.size(), std::cout, "[ displace antennas  % complete ]\n");

            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>>::const_iterator vpit;
            for( vpit = vpAntennas.begin() ; vpit != vpAntennas.end() ; ++vpit ) {
                ++display2;

                std::string country = vpit->first;
                std::string otherCountry = vCountry.front() == country ? vCountry.back() : vCountry.front();
                ign::geometry::LineString antennaGeom = graph.getGeometry(vpit->second.begin()->descriptor);
                ign::geometry::Point dangleEndPoint = vpit->second.begin()->direction == ign::graph::DIRECT? antennaGeom.startPoint() : antennaGeom.endPoint();
                ign::geometry::Point dangleNextPoint = vpit->second.begin()->direction == ign::graph::DIRECT? antennaGeom.pointN(1) : antennaGeom.pointN(antennaGeom.numPoints()-2);


                //DEBUG
                // if( dangleEndPoint.distance(ign::geometry::Point(3845556.01,3072979.21)) < 0.5 ) {
                //     bool test =true;
                // }

                // std::map<std::string, ign::geometry::GeometryPtr>::const_iterator mit = _mCountryGeomWithBuffPtr.find(vpit->first);
                // if (mit == _mCountryGeomWithBuffPtr.end()) {
                //     _logger->log(epg::log::ERROR, "Unknown country [country code] " + vpit->first);
                //     continue;
                // }

                // On ne snap que les edges "étrangers" vers les edges "locaux"
                // if( mit->second->intersects(dangleEndPoint)) continue;

                ign::geometry::MultiLineString mlsEdgesAround;

                ign::geometry::Polygon bbox = dangleEndPoint.getEnvelope().expandBy(2*snapDist).toPolygon();

                ign::feature::FeatureFilter filter(countryCodeName +" = '"+otherCountry+"'");
                epg::tools::FilterTools::addAndConditions(filter, "ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + bbox.toString() + "'),3035))");
                ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(filter);
                while (itEdge->hasNext())
                {
                    ign::feature::Feature const& fEdge = itEdge->next();
                    ign::geometry::LineString const& edgeGeom = fEdge.getGeometry().asLineString();
                    mlsEdgesAround.addGeometry(edgeGeom);
                }

                // regarder d'abord si un overshoot est possible
                app::geometry::tools::LineStringSplitter splitter(antennaGeom);
                splitter.addCuttingGeometry(mlsEdgesAround);
                if (dangleEndPoint == antennaGeom.startPoint()) {
                    std::vector<ign::geometry::LineString> vsubLs = splitter.trimStart();
                    if (vsubLs.size() == 2 && vsubLs.front().length() < snapDist) {
                        antennaGeom = vsubLs.back();
                        continue;
                    }
                } else {
                    std::vector<ign::geometry::LineString> vsubLs = splitter.trimEnd();
                    if (vsubLs.size() == 2 && vsubLs.back().length() < snapDist) {
                        antennaGeom = vsubLs.front();
                        continue;
                    }
                }

                // sinon calculer la projection
                double maxDist = snapDist;
                double maxOrthoDist = snapDist/2;
                ign::geometry::Point maxPt;
                ign::geometry::Point maxOrthoPt;
                
                for ( size_t i = 0 ; i < mlsEdgesAround.numGeometries() ; ++i ) {     
                    ign::geometry::LineString const& edgeGeom = mlsEdgesAround.lineStringN(i);

                    // projection ortho
                    ign::geometry::Point projOrthoPt;
                    epg::tools::geometry::projectZ(edgeGeom, dangleEndPoint, projOrthoPt);
                    double distanceOrtho = dangleEndPoint.distance(projOrthoPt);
                    if (distanceOrtho < maxOrthoDist) {
                        maxOrthoDist = distanceOrtho;
                        maxOrthoPt = projOrthoPt;
                        
                        // snap
                        double startDist = maxOrthoPt.distance(edgeGeom.startPoint());
                        double endDist = maxOrthoPt.distance(edgeGeom.endPoint());
                        if (std::min(startDist, endDist) < snap2EdgeEndingsDist) {
                            maxOrthoPt = startDist < endDist ? edgeGeom.startPoint() : edgeGeom.endPoint();
                        }
                    }

                    // projection axiale
                    bool foundNewPt = false;
                    std::vector< ign::geometry::Point > vPtIntersect = epg::tools::geometry::LineIntersector::compute(dangleEndPoint, dangleNextPoint, edgeGeom);
                    for (std::vector< ign::geometry::Point >::iterator vit2 = vPtIntersect.begin(); vit2 != vPtIntersect.end(); ++vit2) {
                        ign::math::Line2d line(dangleNextPoint.toVec2d(), dangleEndPoint.toVec2d());
                        double abs = line.project(vit2->toVec2d());
                        if (abs < 0) continue;
                        double distance = dangleEndPoint.distance(*vit2);
                        if (distance < maxDist) {
                            maxDist = distance;
                            maxPt = *vit2;
                            foundNewPt = true;
                        }
                    }
                    // snap
                    if (foundNewPt) {
                        double startDist = maxPt.distance(edgeGeom.startPoint());
                        double endDist = maxPt.distance(edgeGeom.endPoint());
                        if (std::min(startDist, endDist) < snap2EdgeEndingsDist) {
                            maxPt = startDist < endDist ? edgeGeom.startPoint() : edgeGeom.endPoint();
                        }
                    }
                }
                if (maxPt.isEmpty())
                    maxPt = maxOrthoPt;

                if (maxPt.isEmpty()) continue;

                if( _inCountry(country, dangleEndPoint) && _inCountry(otherCountry, maxPt) ) continue;

                // deformation  
                std::string edgeId = graph.origins(vpit->second.begin()->descriptor)[0];
                ign::feature::Feature fEdge;
                _fsEdge->getFeatureById(edgeId, fEdge);
                // ign::geometry::LineString edgeGeom = fEdge.getGeometry().asLineString();

                ign::math::Vec2d vectDeform = maxPt.toVec2d() - dangleEndPoint.toVec2d();

                bool deformAtSource = antennaGeom.startPoint().distance(dangleEndPoint) < 1e-5;
                
                ign::math::Vec2d vectSource = deformAtSource ? vectDeform : ign::math::Vec2d();
                ign::math::Vec2d vectTarget = !deformAtSource ? vectDeform : ign::math::Vec2d();
                deformer.deform(vectSource, vectTarget, antennaGeom);

                // remplacement du point extremité
                if (deformAtSource) {
                    double z = antennaGeom.startPoint().z();
                    antennaGeom.startPoint() = maxPt;
                    antennaGeom.startPoint().setZ(z);
                } else {
                    double z = antennaGeom.endPoint().z();
                    antennaGeom.endPoint() = maxPt;
                    antennaGeom.endPoint().setZ(z);
                }

                fEdge.setGeometry(antennaGeom);
                _fsEdge->modifyFeature(fEdge);

                _shapeLogger->writeFeature("ec_projected_antennas", fEdge);
            }

            _loadGraphAndPlanarize(graphManager/*, vCountry[i]*/);
            // std::string otherCountry = vCountry.front() == vCountry[i] ? vCountry.back() : vCountry.front();

            boost::progress_display display3(graph.numLinearOrigins(), std::cout, "[ split edges  % complete ]\n");

            linear_origin_iterator oit, oend; 
            graph.origins( oit, oend );
            for( graph.origins( oit, oend ) ; oit != oend ; ++oit )
            {
                ++display3;

                std::pair< bool, std::vector< oriented_edge_descriptor >> foundInducedEdges = graph.getInducedEdges( *oit );
                if (!foundInducedEdges.first) {
                    _logger->log(epg::log::ERROR, "No induced edges for edge [id] " + *oit);
                    continue;
                }

                //DEBUG
                // if( *oit == "9b39ec5b-4f1a-4706-b07e-06249a4f105e") {
                //     bool testt = true;
                //     ign::geometry::LineString ls1 = graph.getGeometry(foundInducedEdges.second.front());
                //     ign::geometry::LineString ls2 = graph.getGeometry(foundInducedEdges.second.back());
                // }

                if( foundInducedEdges.second.size() == 0 ) continue;

                //TODO a revoir _isCuttingPoint si les edges de pays differents se touchent à leur extremité
                bool sourceIsCuttingPoint = _isCuttingPoint(graphManager, graph.source(foundInducedEdges.second.front()));
                bool targetIsCuttingPoint = _isCuttingPoint(graphManager, graph.target(foundInducedEdges.second.back()));
                if( foundInducedEdges.second.size() == 1 ) {
                    // pour la robustesse : si l'extrémité est un point de coupure il faut mettre la geometrie en cohérence avec la coupure
                    // dont on traite l edge si l une de ses extremites est un point de coupure
                    if( sourceIsCuttingPoint || targetIsCuttingPoint ) {
                        ign::feature::Feature featOrigin;
                        _fsEdge->getFeatureById(*oit, featOrigin);
                        ign::geometry::LineString originGeom = featOrigin.getGeometry().asLineString();

                        if( sourceIsCuttingPoint ){
                            double z = originGeom.startPoint().z();
                            originGeom.startPoint() = graph.getGeometry(graph.source(foundInducedEdges.second.front()));
                            originGeom.startPoint().z() = z;
                        }
                            
                        if( targetIsCuttingPoint ) {
                            double z = originGeom.endPoint().z();
                            originGeom.endPoint() = graph.getGeometry(graph.target(foundInducedEdges.second.back()));
                            originGeom.endPoint().z() = z;
                        }
                            
                        featOrigin.setGeometry(originGeom);
                        _fsEdge->modifyFeature(featOrigin);
                    }
                    continue;
                }

                ign::feature::Feature featOrigin;
                _fsEdge->getFeatureById(*oit, featOrigin);
                std::string country = featOrigin.getAttribute(countryCodeName).toString();
                std::string otherCountry = vCountry.front() == country ? vCountry.back() : vCountry.front();

                ign::geometry::LineString originGeom = featOrigin.getGeometry().asLineString();
                geometry::tools::LengthIndexedLineString originLengthIndexedGeom(originGeom);

                // on fusionne les subEdges si les coupures sont causées par des edges du meme pays
                std::vector< std::pair<vertex_descriptor, vertex_descriptor >> vpSourceTarget;
                // vLs.push_back(graph.getGeometry(foundInducedEdges.second.front()));
                vertex_descriptor currentSource = graph.source(foundInducedEdges.second.front());
                for (size_t i = 1 ; i < foundInducedEdges.second.size() ; ++i) {
                    if( !_isCuttingPoint(graphManager, graph.source(foundInducedEdges.second[i]) ) ) continue;

                    vertex_descriptor source = graph.source(foundInducedEdges.second[i]);
                    vpSourceTarget.push_back(std::make_pair(currentSource, source));
                    currentSource = source;
                }
                vpSourceTarget.push_back(std::make_pair(currentSource, graph.target(foundInducedEdges.second.back())));

                if (vpSourceTarget.size() == 1 && !sourceIsCuttingPoint && !targetIsCuttingPoint) continue;

                std::vector< ign::geometry::LineString > vLs;
                
                double startAbs = 0;
                for ( std::vector< std::pair<vertex_descriptor, vertex_descriptor >>::const_iterator vpit = vpSourceTarget.begin() ; vpit != vpSourceTarget.end() ; ++vpit ) {
                    ign::geometry::Point endPoint = graph.getGeometry(vpit->second);
                    double endAbs = originLengthIndexedGeom.project(endPoint);
                    vLs.push_back(originLengthIndexedGeom.getSubLineString(startAbs, endAbs));

                    if( vpit == vpSourceTarget.begin() ) {
                        ign::geometry::Point startPoint = graph.getGeometry(vpSourceTarget.front().first);

                        double z = vLs.back().startPoint().z();
                        vLs.back().startPoint() = startPoint;
                        vLs.back().startPoint().z() = z;
                    } else {
                        vLs.back().startPoint() = vLs[vLs.size()-2].endPoint();
                    }

                    double z = vLs.back().endPoint().z();
                    vLs.back().endPoint() = endPoint;
                    vLs.back().endPoint().z() = z;
                    
                    //nettoyage
                    _removeMcoordinate(vLs.back());
                    if (vLs.back().numPoints() > 2) {
                        if ( vLs.back().startPoint().distance(vLs.back().pointN(1)) < 0.1 /*TODO: a confirmer*/) {
                            vLs.back().removePointN(1);
                        }
                    }
                    if (vLs.back().numPoints() > 2) {
                        if ( vLs.back().endPoint().distance(vLs.back().pointN(vLs.back().numPoints()-2)) < 0.1 /*TODO: a confirmer*/) {
                            vLs.back().removePointN(vLs.back().numPoints()-2);
                        }
                    }

                    startAbs = endAbs;
                }
                // pour la robustesse
                if (!sourceIsCuttingPoint)
                    vLs.front().startPoint() = originGeom.startPoint();
                if (!targetIsCuttingPoint)
                    vLs.back().endPoint() = originGeom.endPoint();

                _fsEdge->deleteFeature(*oit);

                // double previousZ = originGeom.startPoint().z();
                for (size_t i = 0 ; i < vLs.size() ; ++i) {
                    //on enleve tous les dangles aux extremités
                    if ( 
                        i == 0 && graph.degree(vpSourceTarget[i].first) == 1 ||
                        i == vLs.size()-1 && graph.degree(vpSourceTarget[i].second) == 1
                    ) {
                        ign::geometry::Point dangleEndPoint = i == 0 ? 
                            graph.getGeometry(vpSourceTarget[i].first):
                            graph.getGeometry(vpSourceTarget[i].second);

                        std::map<std::string, ign::geometry::GeometryPtr>::const_iterator mit = _mCountryGeomPtr.find(country);
                        if (mit == _mCountryGeomPtr.end()) {
                            _logger->log(epg::log::ERROR, "Unknown country [country code] " + country);
                        } else {
                            if( (mit->second->distance(dangleEndPoint) > landmaskBuffer) || vLs[i].length() < antennaMinLength) {
							//if (!mit->second->intersects(dangleEndPoint) || vLs[i].length() < antennaMinLength) {
                                ign::feature::Feature feat;
                                feat.setGeometry(vLs[i]);
                                _shapeLogger->writeFeature("ec_trimed_edges_parts", feat);

                                continue;
                            }
                        }
                        // ou si hausdorff avec edges alentours < 5 ?

                        // on recupere les edges voisins appartenant à l'autre pays
                        ign::geometry::MultiLineString mlsOtherCountryEdgesAround;
                        ign::geometry::Envelope bbox = vLs[i].getEnvelope().expandBy(2*antennaMinDist2Neighbor);
                        std::set< edge_descriptor > sEdgesAround;
                        graph.getEdges(bbox, sEdgesAround);
                        for ( std::set< edge_descriptor >::const_iterator sit = sEdgesAround.begin() ; sit != sEdgesAround.end() ; ++sit) {
                            if (graphManager.getCountry(*sit) == otherCountry) {
                                mlsOtherCountryEdgesAround.addGeometry(graph.getGeometry(*sit));
                            }
                        }

                        double hDist = ign::geometry::algorithm::HausdorffDistanceOp::distance(vLs[i], mlsOtherCountryEdgesAround);
                        if (hDist >= 0 && hDist < antennaMinDist2Neighbor)
                            continue;
                    }

                    // vLs[i].startPoint().z() = previousZ;
                    // double nextZ = i != vLs.size()-1 ? _getZ(originLengthIndexedGeom, vLs[i].endPoint()) : originGeom.endPoint().z();;
                    // vLs[i].endPoint().z() = nextZ;
                    // previousZ = nextZ;

                    featOrigin.setGeometry(vLs[i]);
                    _fsEdge->createFeature(featOrigin);

                    _shapeLogger->writeFeature("ec_split_edges", featOrigin);
                }
            }
        }

        ///
        ///
        ///
        bool EdgeConnectorOp::_inCountry(std::string const& country, ign::geometry::Point const& pt ) const {
            std::map<std::string, ign::geometry::GeometryPtr>::const_iterator mit = _mCountryGeomPtr.find(country);
            if (mit == _mCountryGeomPtr.end()) {
                _logger->log(epg::log::ERROR, "Unknown country [country code] " + country);
            } else {
                return mit->second->Intersection(pt);
            }
            return false;
        }

        ///
        ///
        ///
        void EdgeConnectorOp::_removeMcoordinate(ign::geometry::LineString & ls) const {
            ign::geometry::LineString::iterator lsit;
            for ( lsit = ls.begin() ; lsit != ls.end() ; ++lsit ) {
                lsit->m() = ign::numeric::Numeric< double >::NaN();
            }
        }

        ///
        ///
        ///
        bool EdgeConnectorOp::_isCuttingPoint(
            app::calcul::detail::EdgeCleaningGraphManager & graphManager,
            vertex_descriptor v
        ) const {
            std::set<std::string> sCountry;
            std::set<std::string> sOrigin;
            bool hasTwiceSameOrigin = false;

            std::vector< oriented_edge_descriptor > vEdges;
            graphManager.getGraph().incidentEdges( v, vEdges );

            for (size_t i = 0 ; i < vEdges.size() ; ++i) {
                sCountry.insert(graphManager.getCountry(vEdges[i].descriptor));
                std::string origin = graphManager.getGraph().origins(vEdges[i].descriptor)[0];

                if(!hasTwiceSameOrigin) {
                    hasTwiceSameOrigin = sOrigin.find(origin) != sOrigin.end();
                    if(!hasTwiceSameOrigin) sOrigin.insert(origin);
                }
            }
            return hasTwiceSameOrigin && sCountry.size() > 1;
        }

        ///
        ///
        ///
        double EdgeConnectorOp::_getZ( geometry::tools::LengthIndexedLineString const& lengthIndexedLs, ign::geometry::Point const& p) const {
            double abs = lengthIndexedLs.project(p);
            ign::geometry::Point proj = lengthIndexedLs.locateAlong(abs);
            return proj.z();
        }

        ///
        ///
        ///
        void EdgeConnectorOp::_loadGraphAndPlanarize(
            app::calcul::detail::EdgeCleaningGraphManager & graphManager/*,
            std::string const& countryCode*/
            ) const 
        {
            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            // std::string const geomName = epgParams.getValue(GEOM).toString();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            graphManager.clear();

            // std::vector<std::string> vCountry;
		    // epg::tools::StringTools::Split(_borderCode, "#", vCountry);
            // std::string otherCountry = vCountry.front() == countryCode ? vCountry.back() : vCountry.front();

            // std::map<std::string, ign::geometry::GeometryPtr>::const_iterator mit = _mCountryGeomPtr.find(countryCode);
            // if (mit == _mCountryGeomPtr.end()) {
            //     _logger->log(epg::log::ERROR, "Unknown country [country code] " + countryCode);
            //     return;
            // }

            // ign::feature::FeatureFilter filter(countryCodeName +" = '"+otherCountry+"'");
            ign::feature::FeatureFilter filter;
            // epg::tools::FilterTools::addAndConditions(filter, "ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + mit->second->toString() + "'),3035))");

            int numFeatures = epg::sql::tools::numFeatures(*_fsEdge, filter);
            boost::progress_display display(numFeatures, std::cout, "[ edge_loading  % complete ]\n");

            ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(filter);
            while (itEdge->hasNext())
            {
                ++display;
                ign::feature::Feature const& fEdge = itEdge->next();
                ign::geometry::LineString const& edgeGeom = fEdge.getGeometry().asLineString();
                std::string edgeId = fEdge.getId();
                std::string country = fEdge.getAttribute(countryCodeName).toString();
                bool isCl = country.find("#") != std::string::npos;

                //DEBUG
                // if( fEdge.getId() == "9b39ec5b-4f1a-4706-b07e-06249a4f105e") {
                //     bool test = true;
                // }

                graphManager.addEdge(edgeGeom, edgeId, OriginEdgeProperties(country, isCl));

                // ign::feature::FeatureFilter filter2(countryCodeName +" = '"+countryCode+"'");
                // epg::tools::FilterTools::addAndConditions(filter2, "ST_DISTANCE(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + edgeGeom.toString() + "'),3035)) < 0.001");

                // ign::feature::FeatureIteratorPtr itEdge2 = _fsEdge->getFeatures(filter2);
                // while (itEdge2->hasNext())
                // {
                //     ign::feature::Feature const& fEdge2 = itEdge2->next();
                //     ign::geometry::LineString const& edgeGeom2 = fEdge2.getGeometry().asLineString();
                //     std::string edgeId2 = fEdge2.getId();
                //     std::string country2 = fEdge2.getAttribute(countryCodeName).toString();
                //     bool isCl2 = country2.find("#") != std::string::npos;

                //     graphManager.addEdge(edgeGeom2, edgeId2, OriginEdgeProperties(country2, isCl2));
                // }
            }

            graphManager.planarize();
        }

        ///
        ///
        ///
        void EdgeConnectorOp::_loadGraph(
            app::calcul::detail::EdgeCleaningGraphManager & graphManager, 
            bool planarize,
            ign::feature::FeatureFilter filter
            ) const 
        {
            graphManager.clear();

            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            // chargement des edges
            // patience
            int numFeatures = epg::sql::tools::numFeatures(*_fsEdge, filter);
            boost::progress_display display(numFeatures, std::cout, "[ edge_loading  % complete ]\n");

            ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(filter);
            while (itEdge->hasNext())
            {
                ++display;
                ign::feature::Feature const& fEdge = itEdge->next();
                ign::geometry::LineString const& ls = fEdge.getGeometry().asLineString();
                std::string edgeId = fEdge.getId();
                std::string country = fEdge.getAttribute(countryCodeName).toString();
                bool isCl = country.find("#") != std::string::npos;

                if (planarize) {
                    graphManager.addEdge(ls, edgeId, OriginEdgeProperties(country, isCl));
                } else {
                    graphManager.addEdgeSimple(ls, edgeId, OriginEdgeProperties(country, isCl));
                }
            }

            if (planarize) {
                graphManager.planarize();
                graphManager.createFaces();
            }
        }
    }
}

