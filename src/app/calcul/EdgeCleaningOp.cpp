// APP
#include <app/calcul/EdgeCleaningOp.h>
#include <app/params/ThemeParameters.h>
#include <app/calcul/detail/graph/concept/EdgeCleaningGraphSpecializations.h>
#include <app/tools/mergeVertices.h>


// BOOST
#include <boost/progress.hpp>

// EPG
#include <epg/Context.h>
#include <epg/params/EpgParameters.h>
#include <epg/sql/tools/numFeatures.h>
#include <epg/sql/DataBaseManager.h>
#include <epg/tools/StringTools.h>
#include <epg/tools/TimeTools.h>
#include <epg/tools/FilterTools.h>
#include <epg/graph/tools/convertPathToLineString.h>

//OME2
#include <ome2/geometry/tools/GetEndingPointsOp.h>
#include <ome2/geometry/tools/lineStringTools.h>

// SOCLE
#include <ign/geometry/graph/detail/NextEdge.h>
#include <ign/geometry/algorithm/HausdorffDistanceOp.h>
#include <ign/graph/algorithm/DijkstraShortestPaths.h>
#include <ign/math/LineT.h>
#include <ign/math/Line2T.h>


namespace app
{
    namespace calcul
    {
        ///
        ///
        ///
        EdgeCleaningOp::EdgeCleaningOp(
            std::string borderCode,
            bool verbose
        ) : 
            _countryCode(borderCode),
            _verbose(verbose),
            _tag("EdgeCleaningOpTag")
        {
            _init();
        }

        ///
        ///
        ///
        EdgeCleaningOp::~EdgeCleaningOp()
        {
            _shapeLogger->closeShape("ecl_deleted_edges");
            _shapeLogger->closeShape("ecl_country");
            _shapeLogger->closeShape("ecl_paths_out_of_country");
            _shapeLogger->closeShape("ecl_slim_surface");
            _shapeLogger->closeShape("ecl_big_face_removed");
            _shapeLogger->closeShape("ecl_slim_face_1_path");
            _shapeLogger->closeShape("ecl_slim_face_2_path_same_country");
        }

        ///
        ///
        ///
        void EdgeCleaningOp::_init()
        {
            //--
            _logger = epg::log::EpgLoggerS::getInstance();
            _logger->log(epg::log::INFO, "[START] initialization: " + epg::tools::TimeTools::getTime());

            //--
            _shapeLogger = epg::log::ShapeLoggerS::getInstance();
            _shapeLogger->addShape("ecl_deleted_edges", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("ecl_country", epg::log::ShapeLogger::POLYGON);
            _shapeLogger->addShape("ecl_paths_out_of_country", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("ecl_slim_surface", epg::log::ShapeLogger::POLYGON);
            _shapeLogger->addShape("ecl_slim_surface_medial_axis", epg::log::ShapeLogger::LINESTRING);
            _shapeLogger->addShape("ecl_big_face_removed", epg::log::ShapeLogger::POLYGON);
            _shapeLogger->addShape("ecl_slim_face_1_path", epg::log::ShapeLogger::POLYGON);
            _shapeLogger->addShape("ecl_slim_face_2_path_same_country", epg::log::ShapeLogger::POLYGON);

            //--
            epg::Context *context = epg::ContextS::getInstance();

            // epg parameters
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const boundaryTableName = epgParams.getValue(TARGET_BOUNDARY_TABLE).toString();
            std::string const idName = epgParams.getValue(ID).toString();
            std::string const geomName = epgParams.getValue(GEOM).toString();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            std::string const edgeTableName = epgParams.getValue(EDGE_TABLE).toString();
            
            // app parameters
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            std::string cpTableName = themeParameters->getValue(CP_TABLE).toString();
            if ( cpTableName == "" ) {
                std::string const cpTableSuffix = themeParameters->getValue(CP_TABLE_SUFFIX).toString();
                cpTableName = edgeTableName + cpTableSuffix;
            }
            std::string const landmaskTableName = themeParameters->getValue(LANDMASK_TABLE).toString();
            std::string const landCoverTypeName = themeParameters->getValue(LAND_COVER_TYPE).toString();
            std::string const landAreaValue = themeParameters->getValue(TYPE_LAND_AREA).toString();
            std::string clTableName = themeParameters->getValue(CL_TABLE).toString();
            if ( clTableName == "" ) {
                std::string const clTableSuffix = themeParameters->getValue(CL_TABLE_SUFFIX).toString();
                clTableName = edgeTableName + clTableSuffix;
            }
            double const landmaskBuffer = themeParameters->getValue(ECL_LANDMASK_BUFFER).toDouble();

            // on recupere un buffer autour de la frontiere
            ign::geometry::GeometryPtr boundBuffPtr(new ign::geometry::Polygon());
            ign::feature::sql::FeatureStorePostgis* fsBoundary = context->getDataBaseManager().getFeatureStore(boundaryTableName, idName, geomName);
            ign::feature::FeatureIteratorPtr itBoundary = fsBoundary->getFeatures(ign::feature::FeatureFilter(countryCodeName +" = '"+_countryCode+"'"));
            while (itBoundary->hasNext())
            {
                ign::feature::Feature const& fBoundary = itBoundary->next();
                ign::geometry::LineString const& ls = fBoundary.getGeometry().asLineString();

                ign::geometry::GeometryPtr tmpBuffPtr(ls.buffer(10000));

                boundBuffPtr.reset(boundBuffPtr->Union(*tmpBuffPtr));
            }

            //on recupere la geometry des pays
            std::vector<std::string> vCountry;
		    epg::tools::StringTools::Split(_countryCode, "#", vCountry);

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
                _mCountryGeomWithBuffPtr.insert(std::make_pair(*vit, ign::geometry::GeometryPtr(_mCountryGeomPtr[*vit]->buffer(landmaskBuffer))));

                ign::feature::Feature feat;
                feat.setGeometry(*_mCountryGeomPtr[*vit]);
                _shapeLogger->writeFeature("ecl_country", feat);
            }

            //--
            _fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);

            //--
            _fsCp = context->getDataBaseManager().getFeatureStore(cpTableName, idName, geomName);

            //--
            _logger->log(epg::log::INFO, "[END] initialization: " + epg::tools::TimeTools::getTime());
        };

        ///
        ///
        ///
        void EdgeCleaningOp::_loadGraph(
            app::calcul::detail::EdgeCleaningGraphManager & graphManager, 
            bool planarize,
            bool simplifiedPlanarization,
            ign::feature::FeatureFilter filter        
            ) const 
        {
            graphManager.clear();
            graphManager.setSimplifiedPlanarization(simplifiedPlanarization);

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
            graphManager.setSimplifiedPlanarization(false);
        }

        ///
        std::pair<double, double> EdgeCleaningOp::_getLengths( ign::geometry::Geometry const& geom, ign::geometry::Point const* startPoint ) const
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

        ///
        ///
        ///
        std::pair<double, double> EdgeCleaningOp::_addLengths(
            std::string country, 
            ign::geometry::LineString const& ls,
            double & lengthInCountry,
            double & length
        ) const {
            std::map<std::string, ign::geometry::GeometryPtr>::const_iterator mit = _mCountryGeomPtr.find(country);
            if (mit == _mCountryGeomPtr.end()) {
                _logger->log(epg::log::ERROR, "Unknown country [country code] " + country);
                return std::make_pair(0, 0);
            }

            ign::geometry::GeometryPtr resultPtr(mit->second->Intersection(ls));

            std::pair<double, double> lengths = _getLengths(*resultPtr, &ls.startPoint());

            lengthInCountry += lengths.first;
            length += ls.length();

            return std::make_pair(lengths.second, lengths.second/ls.length());
        }

        ///
        ///
        ///
        void EdgeCleaningOp::_addLengthsWithBuff(std::string country, ign::geometry::LineString const& ls , double & lengthInCountry, double & length) const 
        {
            std::map<std::string, ign::geometry::GeometryPtr>::const_iterator mit = _mCountryGeomWithBuffPtr.find(country);
            if (mit == _mCountryGeomWithBuffPtr.end()) {
                _logger->log(epg::log::ERROR, "Unknown country [country code] " + country);
                return;
            }

            ign::geometry::GeometryPtr resultPtr(mit->second->Intersection(ls));

            lengthInCountry += _getLengths(*resultPtr).first;
            length += ls.length();
        }

        ///
        ///
        ///
        double EdgeCleaningOp::_getAntennaLength(GraphType const& graph, std::list<oriented_edge_descriptor> const& lEdges) const
        {
            double length = 0;
            std::list<oriented_edge_descriptor>::const_iterator lit = lEdges.begin();
            for ( ; lit != lEdges.end() ; ++lit) {
                ign::geometry::LineString edgeGeom = graph.getGeometry(lit->descriptor);
                length += edgeGeom.length();
            }
            return length;
        }

        ///
        ///
        ///
        std::pair<double, double> EdgeCleaningOp::_getRatioAndLengthFirstPart(GraphType const& graph, std::string country, std::list<oriented_edge_descriptor> const& lEdges) const
        {
            double lengthInCountry = 0;
            double lengthInCountryFirstPart = 0;
            bool stopAddLengthInCountryFirstPart = false;
            double length = 0;
            std::list<oriented_edge_descriptor>::const_iterator lit = lEdges.begin();
            for ( ; lit != lEdges.end() ; ++lit) {
                ign::geometry::LineString edgeGeom = graph.getGeometry(*lit);
                std::pair<double, double> lengthAnRatioInCountryFirstPart = _addLengths(country, edgeGeom, lengthInCountry, length);
                if (!stopAddLengthInCountryFirstPart) {
                    lengthInCountryFirstPart += lengthAnRatioInCountryFirstPart.first;
                    if (lengthAnRatioInCountryFirstPart.second < 0.99)
                        stopAddLengthInCountryFirstPart = true;
                }
            }
			if (length == 0)
				return std::make_pair(0, 0);

            return std::make_pair(lengthInCountry / length, lengthInCountryFirstPart);
        }

        ///
        ///
        ///
        double EdgeCleaningOp::_getRatioWithBuff(GraphType const& graph, std::string country, std::list<oriented_edge_descriptor> const& lEdges) const
        {
            double lengthInCountry = 0;
            double length = 0;
            std::list<oriented_edge_descriptor>::const_iterator lit = lEdges.begin();
            for ( ; lit != lEdges.end() ; ++lit) {
                ign::geometry::LineString edgeGeom = graph.getGeometry(lit->descriptor);
                _addLengthsWithBuff(country, edgeGeom, lengthInCountry, length);
            }
			if (length == 0)
				return 0;
			
			return lengthInCountry / length;
        }

        ///
        ///
        ///
        void EdgeCleaningOp::cleanFaces() const
        {
            // app parameters
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const slimSurfaceWidth = themeParameters->getValue( ECL_SLIM_SURFACE_WIDTH ).toDouble();

            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager, true);

            std::set<edge_descriptor> sEdge2Remove;

            GraphType & graph = graphManager.getGraph();
            boost::progress_display display(graph.numFaces(), std::cout, "[ cleaning faces  % complete ]\n");
            face_iterator fit, fend;
			for( graph.faces( fit, fend ) ; fit != fend ; ++fit )
			{
                ++display;
				ign::geometry::Polygon faceGeom = graph.getGeometry( *fit );

                // DEBUG
                // _logger->log(epg::log::DEBUG, "coucou");
                _logger->log(epg::log::DEBUG, faceGeom.toString());

				if (_isSlimSurface(faceGeom, slimSurfaceWidth)) {

                    std::map<std::string, std::list<edge_descriptor>> mlEdges;

                    oriented_edge_descriptor startEdge = graph.getIncidentEdge( *fit );
                    oriented_edge_descriptor nextEdge = startEdge;
                    std::string previousCountry = graphManager.getCountry(nextEdge.descriptor);
                    bool endingPointPassed = false;
                    size_t nbPassedEndPoints = 0;
                    std::set<std::string> sHasConnection;
                    bool aborded = false;
                    do{
                        //DEBUG
                        // std::string pouet = graph.origins(nextEdge.descriptor)[0];
                        // if (pouet == "23aca9e2-5fd0-4c4f-823a-1f9c52ac54f9" || pouet == "16071668-add2-4735-bf5d-3d430987428f" || pouet == "03cb00ed-9a26-4b68-9412-4ca33b5f9014") {
                        //     bool test = true;
                        // }

                        // on ne doit pas avoir de cl dans la boucle
                        if (graphManager.isCl(nextEdge.descriptor)) {
                            _logger->log(epg::log::WARN, "Loop contains a CL [cl id] "+graph.origins(nextEdge.descriptor)[0]);
                            aborded = true;
                            break;
                        }

                        std::string country = graphManager.getCountry(nextEdge.descriptor);

                        if (endingPointPassed) {
                            if (previousCountry == country) {
                                _logger->log(epg::log::WARN, "Same country on incident edges [edge id] "+graph.origins(nextEdge.descriptor)[0]);
                                aborded = true;
                                break;
                            }
                            endingPointPassed == false;
                        } else if ( previousCountry != country ) {
                            _logger->log(epg::log::WARN, "Mixed country on path [edge id] "+graph.origins(nextEdge.descriptor)[0]);
                            aborded = true;
                            break;
                        }

                        std::map<std::string, std::list<edge_descriptor>>::iterator mlit = mlEdges.find(country);
                        if (mlit == mlEdges.end()) mlit = mlEdges.insert(std::make_pair(country, std::list<edge_descriptor>())).first;
                        mlit->second.push_back(nextEdge.descriptor);

                        if (graph.degree(graph.target(nextEdge)) > 2) {
                            if (graphManager.isTouchingCl(graph.target(nextEdge))) {
                                endingPointPassed = true;
                                ++nbPassedEndPoints;
                            } else {
                                sHasConnection.insert(country);
                            }
                        }
                        previousCountry = country;
                        nextEdge = ign::geometry::graph::detail::nextEdge( nextEdge, graph );
                    }while( nextEdge != startEdge );

                    if (aborded) continue;
                    if (nbPassedEndPoints != 2) continue;
                    if (sHasConnection.size() > 1) continue;

                    if ( mlEdges.size() > 2 ) {
                        _logger->log(epg::log::WARN, "More than 2 paths found path [edge id] "+graph.origins(startEdge.descriptor)[0]);
                        continue;
                    }

                    if ( mlEdges.size() == 1 ) {
                        _logger->log(epg::log::WARN, "Only 1 path found path [edge id] "+graph.origins(startEdge.descriptor)[0]);
                        continue;
                    }

                    if ( mlEdges.size() == 0 ) {
                        _logger->log(epg::log::WARN, "No path found path [edge id] "+graph.origins(startEdge.descriptor)[0]);
                        continue;
                    }

                    // quel chemin doit-on garder ?
                    double ratio1 = _getRatio(graph, mlEdges.begin()->first, mlEdges.begin()->second);
                    bool hasConnection1 = sHasConnection.find(mlEdges.begin()->first) != sHasConnection.end();
                    double ratio2 = _getRatio(graph, mlEdges.rbegin()->first, mlEdges.rbegin()->second);
                    bool hasConnection2 = sHasConnection.find(mlEdges.rbegin()->first) != sHasConnection.end();

                    if ( ratio1 > ratio2 && !hasConnection2 ) {
                        _removeEdges(graph, mlEdges.rbegin()->second, sEdge2Remove);
                    } else if ( !hasConnection1 ) {
                        _removeEdges(graph, mlEdges.begin()->second, sEdge2Remove);
                    }
                }
			}
            for ( std::set<edge_descriptor>::const_iterator sit = sEdge2Remove.begin() ; sit != sEdge2Remove.end() ; ++sit )
                graph.removeEdge(*sit);
        }

        ///
        ///
        ///
        // bool EdgeCleaningOp::_isSlimSurface( 
        //     ign::geometry::LineString const& closedLs, 
        //     double maxWidth
		// ) const {
        //     // trouver les 2 points les plus distants du contour
        //     // separer le contour en 2 polylignes
        //     // calculer hausdorff entre ces 2 polylignes
        //     int indexI = -1;
        //     int indexJ = -1;

        //     double maxDist = 0;
        //     for( size_t i = 0 ; i < closedLs.numPoints()-1 ; ++i )
        //     {
        //         for( size_t j = i+1 ; j < closedLs.numPoints() ; ++j )
        //         {
        //             double dist = closedLs.pointN(i).distance(closedLs.pointN(j));
        //             if ( dist > maxDist ) {
        //                 maxDist = dist;
        //                 indexI = i;
        //                 indexJ = j;
        //             }
        //         }
        //     }

        //     std::pair<ign::geometry::LineString,ign::geometry::LineString> pLs = _getSubLineStrings(indexI, indexJ, closedLs);

        //     double hausdorffDist = ign::geometry::algorithm::HausdorffDistanceOp::distance(pLs.first, pLs.second);

        //     if (hausdorffDist < 0)
        //         _logger->log(epg::log::ERROR, "Error in slim surface calculation : hausdorff distance < 0");

        //     return hausdorffDist >= 0 && hausdorffDist < maxWidth;
        // }

        ///
        ///
        ///
        bool EdgeCleaningOp::_isSlimSurface( 
            ign::geometry::Polygon const& poly, 
            double maxWidth,
            ign::geometry::Point const ** p1,
            ign::geometry::Point const ** p2
		) const {
            ign::geometry::LineString const& extRing = poly.exteriorRing();

            std::pair<bool, std::pair<size_t, size_t>> foundIndexes = ome2::geometry::tools::GetEndingPointsOp::computeIndex(extRing);

            if (!foundIndexes.first) {
                _logger->log(epg::log::ERROR, "Error in slim surface calculation : error in ending points calculation");
                _logger->log(epg::log::ERROR, poly.toString());
                return false;
            }

            std::pair<bool, std::pair<ign::geometry::LineString,ign::geometry::LineString>> foundSides = ome2::geometry::tools::getSubLineStrings(foundIndexes.second.first, foundIndexes.second.second, extRing);

            if ( !foundSides.first ) {
                _logger->log(epg::log::ERROR, "Error in slim surface calculation : error in sides calculation");
                return false;
            }

            if (p1) *p1 = &extRing.pointN(foundIndexes.second.first);
            if (p2) *p2 = &extRing.pointN(foundIndexes.second.second);

            double hausdorffDist = ign::geometry::algorithm::HausdorffDistanceOp::distance(foundSides.second.first, foundSides.second.second);

            if (hausdorffDist < 0) {
                _logger->log(epg::log::ERROR, "Error in slim surface calculation : hausdorff distance < 0");
                _logger->log(epg::log::ERROR, poly.toString());
            }

            return hausdorffDist >= 0 && hausdorffDist < maxWidth;
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_getFacePaths(
            detail::EdgeCleaningGraphManager const& graphManager, 
            face_descriptor fd, 
            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> & vpCountryEdges
        ) const {
            GraphType const& graph = graphManager.getGraph();

            oriented_edge_descriptor startEdge = graph.getIncidentEdge( fd );
            oriented_edge_descriptor currentEdge = startEdge;
            std::string currentCountry = graphManager.getCountry(currentEdge.descriptor);
            vpCountryEdges.push_back(std::make_pair(currentCountry, std::list<oriented_edge_descriptor>()));

            do{
                // on ne doit pas avoir de cl dans la boucle
                // if (graphManager.isCl(currentEdge.descriptor)) {
                //     _logger->log(epg::log::WARN, "Loop contains a CL [cl id] "+graph.origins(currentEdge.descriptor)[0]);
                //     return false;
                // }

                //DEBUG
                // ign::geometry::LineString lsTemp = graph.getGeometry(currentEdge);

                vpCountryEdges.back().second.push_back(currentEdge);

                oriented_edge_descriptor nextEdge = ign::geometry::graph::detail::nextEdge( currentEdge, graph );
                std::string nextCountry = graphManager.getCountry(nextEdge.descriptor);

                if (graph.degree(graph.target(currentEdge)) > 2 || _vertexIsCp(graph, graph.target(currentEdge)) || currentCountry != nextCountry ) {
                    if (nextEdge != startEdge)
                        vpCountryEdges.push_back(std::make_pair(nextCountry, std::list<oriented_edge_descriptor>()));
                    if ( currentCountry != nextCountry ) //log a garder ?
                        _logger->log(epg::log::WARN, "Mixed country on path [edge id] "+graph.origins( currentEdge.descriptor)[0]);
                }
                currentCountry = nextCountry;
                currentEdge = nextEdge;
            }while( currentEdge != startEdge );

            if (vpCountryEdges.size() > 1 && vpCountryEdges.front().first == vpCountryEdges.back().first) {
                vertex_descriptor v = graph.target(vpCountryEdges.back().second.back());
                if ( graph.degree(v) == 2 && !_vertexIsCp(graph, v) ) {
                    for (std::list<oriented_edge_descriptor>::const_reverse_iterator rlit = vpCountryEdges.back().second.rbegin() ; rlit != vpCountryEdges.back().second.rend() ; ++rlit) {
                        vpCountryEdges.front().second.push_front(*rlit);
                    }
                    vpCountryEdges.pop_back();
                }
            }

            if ( vpCountryEdges.size() == 0 ) {
                _logger->log(epg::log::WARN, "No path found path [edge id] "+graph.origins(startEdge.descriptor)[0]);
                return false;
            }

            return true;
        }

        ///
        ///
        ///
        double EdgeCleaningOp::_getPathLength(
            GraphType const& graph, 
            std::list<oriented_edge_descriptor> const& path
        ) const {
            double length = 0;
            for (std::list<oriented_edge_descriptor>::const_iterator lit = path.begin() ; lit != path.end() ; ++lit) {
                ign::geometry::LineString edgeGeom = graph.getGeometry(lit->descriptor);
                length += edgeGeom.length();
            }
            return length;
        }

        ///
        ///
        ///
        void EdgeCleaningOp::_removePath(
            GraphType & graph, std::list<oriented_edge_descriptor> const& path, 
            std::set<edge_descriptor>& sEdge2Remove
        ) const {
            std::list<edge_descriptor> lEdges;
            for (std::list<oriented_edge_descriptor>::const_iterator lit = path.begin() ; lit != path.end() ; ++lit)
                lEdges.push_back(lit->descriptor);

            _removeEdges(graph, lEdges, sEdge2Remove);
        }

        ///
        ///
        ///
        void EdgeCleaningOp::_removePathAndGraphEdges(
            GraphType & graph,
            std::list<oriented_edge_descriptor> const& path
        ) const {
            std::list<edge_descriptor> lEdges;
            for (std::list<oriented_edge_descriptor>::const_iterator lit = path.begin() ; lit != path.end() ; ++lit)
                lEdges.push_back(lit->descriptor);

            _removeEdgesAndGraphEdges(graph, lEdges);
        }

        ///
        ///
        ///
        void EdgeCleaningOp::cleanFaces2ByCountry() const
        {
            epg::Context* context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            std::vector<std::string> vCountry;
		    epg::tools::StringTools::Split(_countryCode, "#", vCountry);
            for (size_t i = 0 ; i < vCountry.size() ; ++i) {
                cleanFaces2(countryCodeName +" = '"+vCountry[i]+"'");
            }
        }

        ///
        ///
        ///
        void EdgeCleaningOp::cleanFacesAndAntennaByCountry(std::string const& sqlFilter) const
        {
            epg::Context* context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            bool isPlanar = true;
            bool isSimplified = false;

            std::vector<std::string> vCountry;
		    epg::tools::StringTools::Split(_countryCode, "#", vCountry);
            for (size_t i = 0 ; i < vCountry.size() ; ++i) {
                detail::EdgeCleaningGraphManager graphManager;

                ign::feature::FeatureFilter filter;
                if (sqlFilter!="") filter.setPropertyConditions(sqlFilter);
                epg::tools::FilterTools::addAndConditions(filter, countryCodeName +" LIKE '%"+vCountry[i]+"%'");

                std::set<std::string> sTreatedFeatures = _getTreatedFeatures(filter);
                std::set<std::string> sOldTreatedFeatures = sTreatedFeatures;

                _loadGraph(graphManager, isPlanar, isSimplified, filter);

                // std::set<vertex_descriptor> sTreatedDangles;

                _cleanAntennas(graphManager, /*sTreatedDangles,*/ sTreatedFeatures, isPlanar);
                bool bChangeOccured = _cleanFaces2(graphManager);

                while (bChangeOccured) {
                    bChangeOccured = _cleanAntennas(graphManager, /*sTreatedDangles,*/ sTreatedFeatures, isPlanar);
                    if (bChangeOccured)
                        bChangeOccured = _cleanFaces2(graphManager);
                }

                _tagNewTreatedFeatures(sOldTreatedFeatures, sTreatedFeatures);
            }
        }

        ///
        ///
        ///
        std::set<std::string> EdgeCleaningOp::_getTreatedFeatures(ign::feature::FeatureFilter const& filter_) const
        {
            std::set<std::string> sTreatedFeatures;

            app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG).getValue().toString();

            ign::feature::FeatureFilter filter = filter_;
            epg::tools::FilterTools::addAndConditions(filter, wTagName +" = '"+_tag+"'");

            ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(filter);
            while (itEdge->hasNext())
            {
                ign::feature::Feature const& fEdge = itEdge->next();
                std::string edgeId = fEdge.getId();
                sTreatedFeatures.insert(edgeId);
            }

            return sTreatedFeatures;
        }

        ///
        ///
        ///
        void EdgeCleaningOp::_tagNewTreatedFeatures(
            std::set<std::string> const& sOldTreatedFeatures, 
            std::set<std::string> const& sTreatedFeatures
        ) const {
            // epg parameters
            epg::Context* context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const idName = epgParams.getValue(ID).toString();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            //--
            app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG).getValue().toString();

            std::string newTreatedFeatures = "";
            for (std::set<std::string>::const_iterator sit = sTreatedFeatures.begin() ; sit != sTreatedFeatures.end() ; ++sit) {
                if( sOldTreatedFeatures.find(*sit) != sOldTreatedFeatures.end() ) continue;
                newTreatedFeatures += (newTreatedFeatures == "" ? "'" : ",'") + *sit +"'";
            }

            if (newTreatedFeatures != "")
                epg::ContextS::getInstance()->getDataBaseManager().setValueColumn(_fsEdge->getTableName(), wTagName, _tag, idName + " IN ("+newTreatedFeatures+")");
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::cleanFaces2(std::string const& sqlFilter) const
        {
            ign::feature::FeatureFilter filter;
            if (sqlFilter!="") filter.setPropertyConditions(sqlFilter);

            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager, true/*planarize*/, false/*simplified*/, filter);

            return _cleanFaces2(graphManager);
        }

        ///
        ///
        ///
        std::set<std::string> EdgeCleaningOp::_mergeFacePaths(
            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> & vpCountryEdges
        ) const {
            std::set<std::string> sHasConnection;
            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>>::iterator vpit, vpit_begin, vpit_previous;
            vpit_begin = vpCountryEdges.begin();
            vpit_previous = vpit_begin;
            ++vpit_begin;

            for ( vpit = vpit_begin ; vpit != vpCountryEdges.end() ; ) {
                if ( vpit->first == vpit_previous->first) {
                    sHasConnection.insert(vpit->first);
                    std::copy( vpit->second.begin(), vpit->second.end(), std::back_inserter( vpit_previous->second ) );
                    vpit = vpCountryEdges.erase(vpit);
                } else {
                    vpit_previous = vpit;
                    ++vpit;
                }
            }

            if (vpCountryEdges.size() > 1 && vpCountryEdges.front().first == vpCountryEdges.back().first) {
                sHasConnection.insert(vpCountryEdges.front().first);
                std::copy( vpCountryEdges.back().second.rbegin(), vpCountryEdges.back().second.rend(), std::front_inserter( vpCountryEdges.front().second ) );
                vpCountryEdges.pop_back();
            }
            return sHasConnection;
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_cleanFaces2(detail::EdgeCleaningGraphManager & graphManager) const
        {
            bool bChangeOccured = false;

            // TODO nettoyer les artefacts du a la planarization (antennes liees à des faces ?)

            while ( _cleanGraphFaces(graphManager) ) {
                bChangeOccured = true;
            }

            return bChangeOccured;
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_cleanGraphFaces(detail::EdgeCleaningGraphManager & graphManager) const {
            bool bChangeOccured = false;

            // app parameters
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const slimSurfaceWidth = themeParameters->getValue( ECL_SLIM_SURFACE_WIDTH ).toDouble();
            double const slimSurfaceMaxArea = themeParameters->getValue( ECL_SLIM_SURFACE_MAX_AREA ).toDouble();
            double const slimSurfaceMaxNbPoints = themeParameters->getValue( ECL_SLIM_SURFACE_MAX_NB_POINTS ).toDouble();

            std::set<edge_descriptor> sEdge2Remove;

            GraphType & graph = graphManager.getGraph();
            // boost::progress_display display(graph.numFaces(), std::cout, "[ cleaning graph faces  % complete ]\n");

            face_iterator fit, fend;
            for( graph.faces( fit, fend ) ; fit != fend ; ++fit )
			{
                // ++display;

				ign::geometry::Polygon faceGeom = graph.getGeometry( *fit );

                // DEBUG
                _logger->log(epg::log::DEBUG, faceGeom.toString());
                // if (faceGeom.intersects(ign::geometry::Point(4017129.38,3093853.68))) {
                //     bool test = true;
                // }
                // if (faceGeom.intersects(ign::geometry::Point(3989455.72,3143506.25))) {
                //     bool test = true;
                // }
                // double testArea = faceGeom.area();
                // size_t testNbPts = faceGeom.exteriorRing().numPoints();
                // _logger->log(epg::log::DEBUG, ign::data::Double(testArea).toString());
                // _logger->log(epg::log::DEBUG, ign::data::Integer(testNbPts).toString());


                //TODO
                ign::geometry::Point const * p1 = 0;
                ign::geometry::Point const * p2 = 0;

				if (
                    faceGeom.exteriorRing().numPoints() < slimSurfaceMaxNbPoints && //optimisation
                    faceGeom.area() < slimSurfaceMaxArea && //optimisation
                    _isSlimSurface(faceGeom, slimSurfaceWidth, &p1, &p2) // TODO recuperee v1 et v2
                ) {
                   //DEBUG
                    // _logger->log(epg::log::DEBUG, "IS SLIM!");

                    ign::feature::Feature feat;
                    feat.setGeometry(faceGeom);
                    _shapeLogger->writeFeature("ecl_slim_surface", feat);

                    std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> vpCountryEdges;
                    if (!_getFacePaths(graphManager, *fit, vpCountryEdges))
                        continue;
                    
                    if (vpCountryEdges.size() == 1) {
                        ign::feature::Feature feat;
                        feat.setGeometry(faceGeom);
                        _shapeLogger->writeFeature("ecl_slim_face_1_path", feat);

                        _removePath(graph, vpCountryEdges.front().second, sEdge2Remove);
                        bChangeOccured = true;
                        continue;
                    }

                    // traiter dans un premier temps en ne raisonnant que sur les paths
                    // TODO voir si on peut traiter les cas avec un mix 2 pays + #

                    std::set<std::string> sFaceCountries;
                    for (std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>>::const_iterator vpit = vpCountryEdges.begin() ; vpit != vpCountryEdges.end() ; ++vpit)
                        sFaceCountries.insert(vpit->first);

                    std::set<std::string> hasConnection;
                    if ( sFaceCountries.size() != 1 ) {
                        hasConnection = _mergeFacePaths(vpCountryEdges);
                    }

                    bool bUse1stMethode = vpCountryEdges.size() == 2;

                    if (bUse1stMethode) {
                        ign::geometry::LineString lsFront = epg::graph::tools::convertPathToLineString(graph, vpCountryEdges.front().second);
                        ign::geometry::LineString lsBack = epg::graph::tools::convertPathToLineString(graph, vpCountryEdges.back().second);

                        bUse1stMethode =  ign::geometry::algorithm::HausdorffDistanceOp::distance(lsFront, lsBack) <= slimSurfaceWidth;
                    }

                    // contitution des branches
                    std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>>  branch1;
                    std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>>  branch2;
                    bool hasConnection1 = false;
                    bool hasConnection2 = false;

                    if (bUse1stMethode) {
                        branch1.push_back(vpCountryEdges.front());
                        branch2.push_back(vpCountryEdges.back());

                        hasConnection1 = hasConnection.find(vpCountryEdges.front().first) != hasConnection.end();
                        hasConnection2 = hasConnection.find(vpCountryEdges.back().first) != hasConnection.end();
                    }  else {
                        // on recupere pour chacune des 2 branches entre p1 et p2 les paths qui les constituent
                        // si p1 ou p2 ne correspond pas à une extremité de chemin on ne traite pas la face
                        // si p1 et p2 correspondent aux extremites topos
                        // boucler sur les chemins, pour chaque extremité de chemin on calcul la distance à p1 et p2
                        // on selectionne les 2 extremités qui sont les plus proche de p1 et p2 (et dont la distance à ces points est < seuil)
                        double dMin1 = 1e-5;
                        double dMin2 = 1e-5;
                        int startPath1 = -1;
                        int startPath2 = -1;
                        for ( size_t i = 0; i < vpCountryEdges.size() ; ++i) {
                            double d1 = graph.getGeometry(graph.source(*vpCountryEdges[i].second.begin())).distance(*p1);
                            if ( d1 < dMin1) {
                                dMin1 = d1;
                                startPath1 = i;
                            }
                            double d2 = graph.getGeometry(graph.source(*vpCountryEdges[i].second.begin())).distance(*p2);
                            if ( d2 < dMin2) {
                                dMin2 = d2;
                                startPath2 = i;
                            }
                        }

                        if (startPath1 < 0 || startPath2 < 0) continue;

                        size_t pathMin = std::min(startPath1, startPath2);
                        size_t pathMax = std::max(startPath1, startPath2);

                        for ( size_t i = pathMin; i < vpCountryEdges.size() ; ++i) {
                            if ( i < pathMax ) 
                                branch1.push_back(vpCountryEdges[i]);
                            else
                                branch2.push_back(vpCountryEdges[i]);
                        }
                        for ( size_t i = 0; i < pathMin ; ++i) 
                            branch2.push_back(vpCountryEdges[i]);


                        hasConnection1 = _hasConnection(graph, branch1);
                        hasConnection2 = _hasConnection(graph, branch2);    
                    }

                    // pour chaque branche on calcule le ratio et on regarde si il y a des connections
                    double length1 = 0;
                    double ratio1 = 0;
                    for ( size_t i = 0; i < branch1.size() ; ++i) {
                        double ratio = _getRatio(graph, branch1[i].first, branch1[i].second);
                        double length = _getPathLength(graph, branch1[i].second);
                        ratio1 += length*ratio;
                        length1 += length;
                    }
                    ratio1 /= length1;

                    double length2 = 0;
                    double ratio2 = 0;
                    for ( size_t i = 0; i < branch2.size() ; ++i) {
                        double ratio = _getRatio(graph, branch2[i].first, branch2[i].second);
                        double length = _getPathLength(graph, branch2[i].second);
                        ratio2 += length*ratio;
                        length2 += length;
                    }
                    ratio2 /= length2;

                    // on supprime tous les chemins d'une branche si:
                    // - c'est la seule des 2 à n'avoir pas de connnexion
                    // - elle à le moins bon ration et les 2 branches n'ont pas de connexion
                    bool bChangeOccuredTmp = true;
                    if (hasConnection1 && !hasConnection2) {
                        for (size_t i = 0 ; i < branch2.size() ; ++i)
                            _removePath(graph, branch2[i].second, sEdge2Remove);
                    } else if (!hasConnection1 && hasConnection2) {
                        for (size_t i = 0 ; i < branch1.size() ; ++i)
                            _removePath(graph, branch1[i].second, sEdge2Remove);
                    } else if (!hasConnection1 && !hasConnection2) {
                        if ( ratio1 > ratio2 ) 
                            for (size_t i = 0 ; i < branch2.size() ; ++i)
                                _removePath(graph, branch2[i].second, sEdge2Remove);
                        else if ( ratio1 < ratio2 )
                            for (size_t i = 0 ; i < branch1.size() ; ++i)
                                _removePath(graph, branch1[i].second, sEdge2Remove);
                        else {
                            // 2 chemins avec le même ratio, on garde le plus court
                            if (length1 < length2) {
                                for (size_t i = 0 ; i < branch2.size() ; ++i)
                                    _removePath(graph, branch2[i].second, sEdge2Remove);
                            } else {
                                for (size_t i = 0 ; i < branch1.size() ; ++i)
                                    _removePath(graph, branch1[i].second, sEdge2Remove);
                            }
                            ign::feature::Feature feat;
                            feat.setGeometry(faceGeom);
                            _shapeLogger->writeFeature("ecl_slim_face_2_path_same_country", feat);
                        }
                    } else {
                        bChangeOccuredTmp = false;
                    }
                    if (bChangeOccuredTmp) bChangeOccured = true;
                } else {
                    // grande face
                    // si dans mauvais pays on supprime

                    std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> vpCountryEdges;
                    if (!_getFacePaths(graphManager, *fit, vpCountryEdges))
                        continue;

                    if (vpCountryEdges.size() == 1) {
                        double ratio = _getRatio(graph, vpCountryEdges.front().first, vpCountryEdges.front().second);

                        if(ratio == 0) {
                            ign::feature::Feature feat;
                            feat.setGeometry(faceGeom);
                            _shapeLogger->writeFeature("ecl_big_face_removed", feat);

                            _removePath(graph, vpCountryEdges.front().second, sEdge2Remove);
                            bChangeOccured = true;
                        }
                    }
                }
			}

            std::set<vertex_descriptor> sVertices;
            for ( std::set<edge_descriptor>::const_iterator sit = sEdge2Remove.begin() ; sit != sEdge2Remove.end() ; ++sit ) {
                if (graph.degree(graph.source(*sit)) > 2)
                    sVertices.insert(graph.source(*sit));
                if (graph.degree(graph.target(*sit)) > 2)
                    sVertices.insert(graph.target(*sit));
                graph.removeEdge(*sit);
            }

            _cleanFacesAntennas(graphManager, sVertices);

            return bChangeOccured;
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_hasConnection (GraphType const& graph, std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>>  const& branch) const {
            for ( size_t i = 0; i < branch.size() ; ++i) {
                for ( std::list<oriented_edge_descriptor>::const_iterator lit = branch[i].second.begin() ; lit != branch[i].second.end() ; ++lit ) {
                    if ( i == branch.size()-1 && std::next(lit) == branch[i].second.end() ) break;
                    if ( graph.degree(graph.target(*lit)) > 2 || _vertexIsCp(graph, graph.target(*lit)) ) {
                        return true;
                    }
                }
            }
            return false;
        };

        ///
        ///
        ///
        void EdgeCleaningOp::_cleanFacesAntennas(
            detail::EdgeCleaningGraphManager & graphManager, 
            std::set<vertex_descriptor> const& sVertices
        ) const {
            GraphType & graph = graphManager.getGraph();
            boost::progress_display display(graph.numFaces(), std::cout, "[ cleaning faces antennas % complete ]\n");

            std::set<std::string> sTreatedFeatures;

            std::set<vertex_descriptor>::const_iterator sit;
            for (sit = sVertices.begin() ; sit != sVertices.end() ; ++sit) {
                if( graph.degree( *sit ) != 1 ) continue;
                if( _vertexIsCp(graph, *sit) ) continue;

                std::pair<bool, std::list<oriented_edge_descriptor>> pAntenna = _getAntenna(graphManager, *sit, sTreatedFeatures, true/*is planar*/);

                _removePathAndGraphEdges(graph, pAntenna.second);
            }
        }

        ///
        ///
        ///
        std::pair<bool, std::list<app::calcul::detail::EdgeCleaningGraphManager::oriented_edge_descriptor>> EdgeCleaningOp::_getAntenna(
            detail::EdgeCleaningGraphManager const& graphManager,
            vertex_descriptor v,
            std::set<std::string> & sTreatedFeatures,
            bool isPlanarGraph
        ) const {
            GraphType const& graph = graphManager.getGraph();

            std::list<oriented_edge_descriptor> lAntennaEdges;
            bool isConnected2CF = false;

            std::vector< oriented_edge_descriptor > vEdges;
            graph.incidentEdges( v, vEdges );

            oriented_edge_descriptor nextEdge = vEdges.front(); // si nextEdge n'est pas une CL ?
            if (graphManager.isCl(nextEdge.descriptor)) {
                _logger->log(epg::log::WARN, "Antenna is connecting line [cl id] "+graph.origins(nextEdge.descriptor)[0]);
                return std::make_pair(isConnected2CF, lAntennaEdges);
            }

            std::string edgeFeatId = graph.origins(nextEdge.descriptor)[0];
            if (sTreatedFeatures.find(edgeFeatId) != sTreatedFeatures.end())
                return std::make_pair(isConnected2CF, lAntennaEdges);
            sTreatedFeatures.insert(edgeFeatId);

            std::string country = graphManager.getCountry(nextEdge.descriptor);

            vertex_descriptor vTarget = GraphType::nullVertex();
            
            while( true )
            {
                if (graphManager.getCountry(nextEdge.descriptor) != country) {
                    isConnected2CF = true;
                    break;
                }

                //DEBUG
                // ign::geometry::LineString ls = graph.getGeometry(nextEdge.descriptor);
                // if (ls.distance(ign::geometry::Point(3800840.443,3131644.081)) < 4) {
                //     bool test = true;
                // }

                _addAntennaEdges(graph, nextEdge, lAntennaEdges, isPlanarGraph);

                vTarget = _getTarget(graph, nextEdge, isPlanarGraph);

                if (graphManager.isCl(nextEdge.descriptor) && graph.degree(graph.source( nextEdge )) == 2 ) {
                    _logger->log(epg::log::WARN, "Antenna connected to connecting line [cl id] "+graph.origins(nextEdge.descriptor)[0]);
                    isConnected2CF = true;
                    break;
                }

                bool targetIsCp = _vertexIsCp(graph, vTarget);
                if( targetIsCp )
                    isConnected2CF = true;

                if( graph.degree( vTarget ) != 2 || targetIsCp) { // ou si nextEdge est une CL ?
                    if( graph.degree( vTarget ) == 1 /*antenne isolee*/)
                    {
                        sTreatedFeatures.insert(graph.origins(nextEdge.descriptor)[0]);
                    }
                    break;
                }

                nextEdge = _getNextEdge(graph, nextEdge.descriptor, vTarget, isPlanarGraph);
            }
            return std::make_pair(isConnected2CF, lAntennaEdges);
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_vertexIsConnected2Cl(detail::EdgeCleaningGraphManager const& graphManager, vertex_descriptor v) const
        {
            std::vector< oriented_edge_descriptor > vEdges;
            graphManager.getGraph().incidentEdges( v, vEdges );
            bool hasincidentCl = false;
            bool hasIncidentNotCl = false;
            for (std::vector< oriented_edge_descriptor >::const_iterator vit = vEdges.begin() ; vit != vEdges.end() ; ++vit) {
                if (graphManager.isCl(vit->descriptor)) {
                    hasincidentCl = true;
                } else {
                    hasIncidentNotCl = true;
                }
            }
            return hasincidentCl && hasIncidentNotCl;
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_vertexIsCp(GraphType const& graph, vertex_descriptor v) const
        {
            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const geomName = epgParams.getValue(GEOM).toString();
            
            double threshold = 1e-5;
            ign::geometry::Point const& vGeom = graph.getGeometry(v);

            // ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + ign::geometry::GeometryPtr(vGeom.buffer(threshold))->toString() + "'),3035))");
            ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_GeomFromText('" + ign::geometry::GeometryPtr(vGeom.buffer(threshold))->toString() + "'))");
            ign::feature::FeatureIteratorPtr itCp = _fsCp->getFeatures(filter);

            double dMax = threshold;
            ign::geometry::Point closestCpGeom;
            while (itCp->hasNext())
            {
                ign::feature::Feature const& fCp = itCp->next();
                ign::geometry::Point const& cpGeom = fCp.getGeometry().asPoint();

                double distance = vGeom.distance(cpGeom);
                if (distance < dMax ) {
                    dMax = distance;
                    closestCpGeom = cpGeom;
                }
            }

            if (!closestCpGeom.isNull()) {
                std::set< vertex_descriptor > sVertices;
		        graph.getVertices( closestCpGeom.getEnvelope().expandBy( threshold ), sVertices );

                double dMax2 = threshold;
                vertex_descriptor closestVertex = GraphType::nullVertex();
                for (std::set< vertex_descriptor >::const_iterator sit = sVertices.begin() ; sit != sVertices.end() ; ++sit) {
                    ign::geometry::Point const& vGeom2 = graph.getGeometry(*sit);

                    double distance2 = vGeom2.distance(closestCpGeom);
                    if (distance2 < dMax2 ) {
                        dMax2 = distance2;
                        closestVertex = *sit;
                    }
                }
                if (closestVertex == v) return true;
            }
            return false;
        }

        ///
        ///
        ///
        void EdgeCleaningOp::cleanPathsOutOfCountry() const
        {
            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const geomName = epgParams.getValue(GEOM).toString();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();


            std::map<std::string, ign::geometry::GeometryPtr>::const_iterator mit;
            for (mit = _mCountryGeomPtr.begin() ; mit != _mCountryGeomPtr.end() ; ++mit) {
                ign::feature::FeatureFilter filter("NOT ST_INTERSECTS(ST_LineInterpolatePoint(" + geomName + ", 0.5), ST_SetSRID(ST_GeomFromText('" + mit->second->toString() + "'),3035))");
                epg::tools::FilterTools::addAndConditions(filter, countryCodeName +" = '"+mit->first+"'" );
                epg::tools::FilterTools::addOrConditions(filter, countryCodeName +" LIKE '%#%'" );

                detail::EdgeCleaningGraphManager graphManager;

                _loadGraph(graphManager, false, false, filter);
                graphManager.initWeight();

                GraphType & graph = graphManager.getGraph();

                std::list<vertex_descriptor> lEndingVertices;
                // on recherche les ending vertices
                vertex_iterator vit, vend;
                for( graph.vertices( vit, vend ) ; vit != vend ; ++vit )
                {
                    if (_vertexIsConnected2Cl(graphManager, *vit) || _vertexIsCp(graph, *vit))
                        lEndingVertices.push_back(*vit);
                }

                //on retire les cl du graph pour le calcul de chemin
                edge_iterator eit, eend;
                std::set<edge_descriptor> sEdge2remove;
                for( graph.edges( eit, eend ) ; eit != eend ; ++eit )
                    if (graphManager.isCl(*eit)) sEdge2remove.insert(*eit);
                for (std::set<edge_descriptor>::const_iterator sit = sEdge2remove.begin() ; sit != sEdge2remove.end() ; ++sit )
                    graph.removeEdge(*sit);

				if (lEndingVertices.size() == 0)
					return;


                boost::progress_display display(lEndingVertices.size()-1, std::cout, "[ cleaning paths out of country  % complete ]\n");

                ign::graph::algorithm::DijkstraShortestPaths<GraphType> dijkstraOp(graph, true);
                std::set<edge_descriptor> sEdges;
                for (std::list<vertex_descriptor>::const_iterator lit1 = lEndingVertices.begin() ; lit1 != --lEndingVertices.end() ; ++lit1 ) {
                    ++display;
                    std::list<vertex_descriptor>::const_iterator lit2 = lit1;
                    ++lit2;
                    for ( ; lit2 != lEndingVertices.end() ; ++lit2 ) {
                        edges_path path;
                        dijkstraOp.getPath(*lit1, *lit2, path);

                        for (edges_path_const_iterator pit = path.begin() ; pit != path.end() ; ++pit) {
                            sEdges.insert(pit->descriptor);

                            ign::feature::Feature feat;
                            feat.setGeometry(graph.getGeometry(pit->descriptor));
                            _shapeLogger->writeFeature("ecl_paths_out_of_country", feat);
                        }
                    }
                }
                _removeEdgesAndGraphEdges(graph, std::list<edge_descriptor>(sEdges.begin(), sEdges.end()));
            }
        }

        //
        ///
        ///
        bool EdgeCleaningOp::cleanAntennas() const
        {
            ign::feature::FeatureFilter filter;

            std::set<std::string> sTreatedFeatures = _getTreatedFeatures(filter);
            std::set<std::string> sOldTreatedFeatures = sTreatedFeatures;

            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager, false);

            bool bChangeOccured = _cleanAntennas(graphManager, /*sTreatedDangles,*/ sTreatedFeatures);

            _tagNewTreatedFeatures(sOldTreatedFeatures, sTreatedFeatures);

            return bChangeOccured;
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_cleanAntennas(detail::EdgeCleaningGraphManager & graphManager, /*std::set<vertex_descriptor> & sTreatedDangles,*/ std::set<std::string> & sTreatedFeatures, bool isPlanarGraph) const
        {
            bool bChangeOccured = false;

            while ( _cleanGraphAntennas(graphManager, sTreatedFeatures, isPlanarGraph) ) {
                bChangeOccured = true;
            }

            return bChangeOccured;
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_cleanGraphAntennas(
            detail::EdgeCleaningGraphManager & graphManager, 
            /*std::set<vertex_descriptor> & sTreatedDangles,*/ 
            std::set<std::string> & sTreatedFeatures,
            bool isPlanarGraph
        ) const {
            GraphType & graph = graphManager.getGraph();

            // liste des vertex auxquels sont rattachées une ou plusieurs antennes
            std::map<vertex_descriptor, std::vector<std::list<oriented_edge_descriptor>>> mVertexAntennas;
            std::set<vertex_descriptor> sVerticesConnected2CF;

            boost::progress_display display(graph.numVertices(), std::cout, "[ cleaning graph antennas (1) % complete ]\n");
            vertex_iterator vit, vend;
            for( graph.vertices( vit, vend ) ; vit != vend ; ++vit )
            {
                ++display;

                if( graph.degree( *vit ) != 1 ) continue;

                if( _vertexIsCp(graph, *vit) ) continue;

                std::pair<bool, std::list<oriented_edge_descriptor>> pAntenna = _getAntenna(graphManager, *vit, sTreatedFeatures, isPlanarGraph);

                if (!pAntenna.second.empty()) {
                    vertex_descriptor vd = graph.target(pAntenna.second.back());
                    std::map<vertex_descriptor, std::vector<std::list<oriented_edge_descriptor>>>::const_iterator mit = mVertexAntennas.find(vd);
                    if( mit == mVertexAntennas.end() ) {
                        mVertexAntennas.insert(std::make_pair(vd, std::vector<std::list<oriented_edge_descriptor>>()));
                    }
                    mVertexAntennas[vd].push_back(pAntenna.second);

                    if (pAntenna.first)
                        sVerticesConnected2CF.insert(vd);
                }
            }

            bool bHasRemovedAntenna = false;

            boost::progress_display display2(mVertexAntennas.size(), std::cout, "[ cleaning graph antennas (2) % complete ]\n");
            std::map<vertex_descriptor, std::vector<std::list<oriented_edge_descriptor>>>::iterator mit;
            for (mit = mVertexAntennas.begin() ; mit != mVertexAntennas.end() ; mit++) {
                ++display2;
                // _logger->log(epg::log::DEBUG, "pouet1");
                // _logger->log(epg::log::DEBUG, graph.getGeometry(mit->first).toString());
                bool bConnected2CF = sVerticesConnected2CF.find(mit->first) != sVerticesConnected2CF.end() || _vertexIsConnected2Cl(graphManager, mit->first);
                // _logger->log(epg::log::DEBUG, "pouet2");

                std::vector<std::list<oriented_edge_descriptor>>::const_iterator minVit;
                bool isRemoved = false;
                do {
                    // _logger->log(epg::log::DEBUG, "pouet3");
                    std::vector<std::list<oriented_edge_descriptor>>::const_iterator vit;
                    if ( mit->second.size() > 1) {
                        // _logger->log(epg::log::DEBUG, "pouet4");
                        double minLength = std::numeric_limits<double>::max();
                        for (vit = mit->second.begin() ; vit != mit->second.end() ; ++vit) {
                            // _logger->log(epg::log::DEBUG, "pouet5");
                            // _logger->log(epg::log::DEBUG, graph.getGeometry(graph.target(vit->back())).toString());
                            // _logger->log(epg::log::DEBUG, graph.getGeometry(vit->front()).toString());
                            double length = _getAntennaLength(graph, *vit);
                            // _logger->log(epg::log::DEBUG, "pouet6");
                            if (length < minLength) {
                                minLength = length;
                                minVit = vit;
                            }
                        }
                    } else {
                        // _logger->log(epg::log::DEBUG, "pouet7");
                        minVit = mit->second.begin();
                    }
                    // _logger->log(epg::log::DEBUG, "pouet8");
                    _logger->log(epg::log::DEBUG, graph.getGeometry(minVit->begin()->descriptor).toString());
                    // _logger->log(epg::log::DEBUG, "pouet81");

                    std::string country = graphManager.getCountry(minVit->begin()->descriptor);
                    

                    // _logger->log(epg::log::DEBUG, "pouet9");
                    
                    isRemoved = _cleanAntenna(graph, country, *minVit, bConnected2CF);

                    // _logger->log(epg::log::DEBUG, "pouet10");
                    if (isRemoved) {
                         bHasRemovedAntenna = true;
                         mit->second.erase(minVit);
                    }
                    // _logger->log(epg::log::DEBUG, "pouet11");
                } while (isRemoved && graph.degree(mit->first) > 2 && mit->second.size() > 0);
                // _logger->log(epg::log::DEBUG, "pouet12");
            }
            // _logger->log(epg::log::DEBUG, "pouet13");

            return bHasRemovedAntenna;
        }

        ///
        ///
        ///
        typename app::calcul::detail::EdgeCleaningGraphManager::GraphType::oriented_edge_descriptor EdgeCleaningOp::_getNextEdge(GraphType const& graph, edge_descriptor e, vertex_descriptor vTarget, bool isPlanarGraph) const {

            if (isPlanarGraph) {
                std::string origin = graph.origins(e)[0];
                std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdges = graph.getInducedEdges(origin);

                if (foundInducedEdges.first) {
                    if (foundInducedEdges.second.size() > 1) {
                        if ( foundInducedEdges.second.begin()->descriptor == e ) 
                            e = foundInducedEdges.second.rbegin()->descriptor;
                        else
                            e = foundInducedEdges.second.begin()->descriptor;
                    }  
                }
            } 

            std::vector< oriented_edge_descriptor > vIncEdges;
            graph.incidentEdges( vTarget, vIncEdges );

            return vIncEdges.front().descriptor == e ? vIncEdges.back():vIncEdges.front();
        }

        ///
        ///
        ///
        void EdgeCleaningOp::_addAntennaEdges(GraphType const& graph, oriented_edge_descriptor oe, std::list<oriented_edge_descriptor> & lEdges, bool isPlanarGraph) const {
            if (isPlanarGraph) {
                std::string origin = graph.origins(oe.descriptor)[0];
                std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdges = graph.getInducedEdges(origin);

                if (foundInducedEdges.first) {
                    bool isSameDirection = true;
                    for ( size_t i = 0 ; i < foundInducedEdges.second.size() ; ++i ) {
                        if (foundInducedEdges.second[i].descriptor == oe.descriptor) {
                            isSameDirection = foundInducedEdges.second[i].direction == oe.direction;
                            break;
                        }
                    }
                    if (!isSameDirection)
                        foundInducedEdges.second = _getReversePath(foundInducedEdges.second);

                    for ( std::vector<oriented_edge_descriptor>::const_iterator vit = foundInducedEdges.second.begin() ; vit != foundInducedEdges.second.end() ; ++vit ) {
                        lEdges.push_back(*vit);
                    }
                }
            } else {
                lEdges.push_back(oe);
            }
        }

        ///
        ///
        ///
        typename app::calcul::detail::EdgeCleaningGraphManager::GraphType::vertex_descriptor EdgeCleaningOp::_getTarget(GraphType const& graph, oriented_edge_descriptor oe, bool isPlanarGraph) const {
            if (isPlanarGraph) {
                std::string origin = graph.origins(oe.descriptor)[0];
                std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdges = graph.getInducedEdges(origin);
                if (foundInducedEdges.first) {
                    return graph.source(oe) == graph.source(*foundInducedEdges.second.begin()) ? 
                            graph.target(*foundInducedEdges.second.rbegin()) : 
                            graph.source(*foundInducedEdges.second.begin());
                }
            } else {
                return graph.target( oe );
            }
            return app::calcul::detail::EdgeCleaningGraphManager::GraphType::nullVertex();
        }
        
        ///
        ///
        ///
        bool EdgeCleaningOp::_cleanAntenna(
            GraphType & graph,
            std::string const& country,
            std::list<oriented_edge_descriptor> const& lAntennaEdges,
            bool bAntennaIsConnected2CF
        ) const {
            // app parameters
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const antennaRatioThreshold = themeParameters->getValue( ECL_ANTENNA_RATIO_THRESHOLD ).toDouble();
            double const antennaRatioThresholdWithBuff = themeParameters->getValue( ECL_ANTENNA_RATIO_WITH_BUFFER_THRESHOLD ).toDouble();
            double const antennaMinLength = themeParameters->getValue( ECL_ANTENNA_MIN_LENGTH ).toDouble();
            double const antennaMinLengthInCountry = themeParameters->getValue( ECL_ANTENNA_MIN_LENGTH_IN_COUNTRY ).toDouble();
            
            double antennaLength = _getAntennaLength(graph, lAntennaEdges);
            if (bAntennaIsConnected2CF && antennaLength < antennaMinLength) {
                _removePathAndGraphEdges(graph, lAntennaEdges);
                return true;
            }
            
            // //DEBUG
            // _logger->log(epg::log::DEBUG, "Antenna info :");
            // std::set<std::string> sOrigins;
            // for (std::list<oriented_edge_descriptor>::const_iterator lit = lAntennaEdges.begin() ; lit != lAntennaEdges.end() ; ++lit) {
            //     for (size_t j = 0 ; j < graph.origins(*lit).size() ; ++j) {
            //         sOrigins.insert( graph.origins(*lit)[j] );
            //     }
            // }
            // _logger->log(epg::log::DEBUG, tools::StringTools::toString(sOrigins));

            std::list<oriented_edge_descriptor> reverseAntenna = _getReversePath(lAntennaEdges);

            std::pair<double, double> pRatioLengthFirstPart = _getRatioAndLengthFirstPart(graph, country, reverseAntenna);
            // //DEBUG
            // _logger->log(epg::log::DEBUG, std::to_string(pRatioLengthFirstPart.first));

            if (pRatioLengthFirstPart.second >= antennaMinLengthInCountry)
                return false;

            bool bRemovedAntenna = false;
            if (pRatioLengthFirstPart.first < antennaRatioThreshold) {
                _removePathAndGraphEdges(graph, lAntennaEdges);
                bRemovedAntenna = true;
            } 
            else {
                double ratioWithBuff = _getRatioWithBuff(graph, country, lAntennaEdges);

                if (ratioWithBuff < antennaRatioThresholdWithBuff) {
                    _removePathAndGraphEdges(graph, lAntennaEdges);
                    bRemovedAntenna = true;
                }
            }
            return bRemovedAntenna;
        }

        ///
        ///
        ///
        void EdgeCleaningOp::cleanTinyEdges() const
        {
            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager, false);
            GraphType & graph = graphManager.getGraph();

            while ( _cleanTinyEdges(graph) ) {
            }  
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_cleanTinyEdges( GraphType & graph ) const
        {
            // app parameters
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const maxLength = themeParameters->getValue( ECL_TINY_EDGE_MAX_LENGTH ).toDouble();

            std::map<edge_descriptor, edge_descriptor> mOldNewEdges;
            std::set<edge_descriptor> sEdge2Remove;
		    std::set<vertex_descriptor> sVertice2Remove;

            std::list<edge_descriptor> lEdges;
            edge_iterator eit, eend;
            for (graph.edges(eit, eend); eit != eend; ++eit)
                lEdges.push_back(*eit);

            boost::progress_display display(lEdges.size(), std::cout, "[ cleaning tiny edges  % complete ]\n");
            for (std::list<edge_descriptor>::const_iterator lit = lEdges.begin() ; lit != lEdges.end() ; ++lit)
            {
                ++display;
                if ( sEdge2Remove.find(*lit) != sEdge2Remove.end() ) continue;

                ign::geometry::LineString edgeGeom = graph.getGeometry(*lit);

                if( edgeGeom.length() > maxLength ) continue; 

                vertex_descriptor vRef = graph.degree(graph.target(*lit)) > graph.degree(graph.source(*lit)) ? graph.target(*lit) : graph.source(*lit);
                vertex_descriptor v2Merge = graph.degree(graph.target(*lit)) > graph.degree(graph.source(*lit)) ? graph.source(*lit) : graph.target(*lit);

                tools::mergeVertices( graph, v2Merge, vRef, mOldNewEdges, sEdge2Remove, sVertice2Remove );
            }

            if ( mOldNewEdges.size() == 0 ) return false;

            _persistEdges(graph, mOldNewEdges, sEdge2Remove);

            for ( std::set<edge_descriptor>::const_iterator sit = sEdge2Remove.begin() ; sit != sEdge2Remove.end() ; ++sit )
                graph.removeEdge(*sit);

            for ( std::set<vertex_descriptor>::const_iterator sit = sVertice2Remove.begin() ; sit != sVertice2Remove.end() ; ++sit )
                graph.removeVertex(*sit);
            
            return true;
        }

        ///
        ///
        ///
        void EdgeCleaningOp::_persistEdges( 
            GraphType & graph,
            std::map<edge_descriptor, edge_descriptor> const& mOldNewEdges,
            std::set<edge_descriptor> & sEdge2remove
        ) const {
            // app parameters
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const maxLength = themeParameters->getValue( ECL_TINY_EDGE_MAX_LENGTH ).toDouble();

            for ( std::map<edge_descriptor, edge_descriptor>::const_iterator mit = mOldNewEdges.begin() ; mit != mOldNewEdges.end() ; ++mit ) {
                ign::geometry::LineString edgeGeom = graph.getGeometry(mit->second);
                // _logger->log(epg::log::DEBUG, edgeGeom.toString());

                std::string edgeId = graph.origins(mit->second)[0];
                // _logger->log(epg::log::DEBUG, edgeId);

                if( graph.source(mit->second) ==  graph.target(mit->second) && edgeGeom.length() < maxLength ) {
                    _fsEdge->deleteFeature(edgeId);
                    sEdge2remove.insert(mit->second);
                } else {
                    ign::feature::Feature fEdge;
                    _fsEdge->getFeatureById(edgeId, fEdge);

                    fEdge.setGeometry(edgeGeom);

                    _fsEdge->modifyFeature(fEdge);
                } 
            }
        }

        ///
        ///
        ///
        void EdgeCleaningOp::cleanParalelleEdges() const
        {
            // app parameters
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const maxDist = themeParameters->getValue( ECL_PARALELLE_EDGE_MAX_DIST ).toDouble();

            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager, false);
            GraphType & graph = graphManager.getGraph();

            std::set<edge_descriptor> sEdge2Remove;

            std::set<edge_descriptor> sVisitedEdge;

            boost::progress_display display4(graph.numEdges(), std::cout, "[ cleaning parallele edges  % complete ]\n");
            edge_iterator eit, eend;
            for (graph.edges(eit, eend); eit != eend; ++eit)
            {
                ++display4;
                if ( sVisitedEdge.find(*eit) != sVisitedEdge.end() ) continue;
                std::vector< oriented_edge_descriptor > vParallelEdges;
                graph.edges( graph.source(*eit), graph.target(*eit), vParallelEdges );
                if(vParallelEdges.size() < 2) continue;

                double maxRatio = -1;
                edge_descriptor maxEdge;
                std::set<edge_descriptor> sVisitedParallelEdge;
                for (std::vector< oriented_edge_descriptor >::const_iterator vit = vParallelEdges.begin() ; vit != vParallelEdges.end() ; ++vit ) {
                    if ( sVisitedParallelEdge.find(vit->descriptor) != sVisitedParallelEdge.end() ) continue; /*gestion des boucles*/
                    sVisitedParallelEdge.insert(vit->descriptor);
                    std::string ownCountry = graphManager.getCountry(vit->descriptor);
                    double ratio = graphManager.isCl(vit->descriptor) ? 1.1 : _getRatio(graph, ownCountry, std::list<edge_descriptor>(1, vit->descriptor));
                    if (ratio > maxRatio) {
                        maxRatio = ratio;
                        maxEdge = vit->descriptor;
                    }
                }

                ign::geometry::LineString lsRef = graph.getGeometry(maxEdge);
                std::set<edge_descriptor> sVisitedTemp;
                sVisitedTemp.insert(maxEdge);
                for (std::vector< oriented_edge_descriptor >::const_iterator vit = vParallelEdges.begin() ; vit != vParallelEdges.end() ; ++vit ) {
                    if ( sVisitedTemp.find(vit->descriptor) != sVisitedTemp.end() ) continue; /*gestion des boucles + sert a passer le 1er element*/
                    sVisitedTemp.insert(vit->descriptor);

                    ign::geometry::LineString ls = graph.getGeometry(vit->descriptor);
                    
                    if ( ign::geometry::algorithm::HausdorffDistanceOp::distance(lsRef, ls) < maxDist ) {
                        _removeEdges(graph, std::list<edge_descriptor>(1,vit->descriptor), sEdge2Remove);
                    }
                }
                sVisitedEdge.insert(sVisitedTemp.begin(), sVisitedTemp.end());
            }

            for ( std::set<edge_descriptor>::const_iterator sit = sEdge2Remove.begin() ; sit != sEdge2Remove.end() ; ++sit )
                graph.removeEdge(*sit);
        }
    }
}