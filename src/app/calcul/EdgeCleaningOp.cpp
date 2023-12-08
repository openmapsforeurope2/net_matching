// APP
#include <app/calcul/EdgeCleaningOp.h>
#include <app/params/ThemeParameters.h>
#include <app/calcul/detail/graph/concept/EdgeCleaningGraphSpecializations.h>

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
        void EdgeCleaningOp::clean(
            std::string borderCode,
            bool verbose)
        {
            EdgeCleaningOp op(borderCode, verbose);
            op.cleanAll();
        }

        ///
        ///
        ///
        EdgeCleaningOp::EdgeCleaningOp(
            std::string borderCode,
            bool verbose
        ) : 
            _countryCode(borderCode),
            _verbose(verbose)
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
            std::string const cpTableName = themeParameters->getValue(CP_TABLE).toString();
            std::string const landmaskTableName = themeParameters->getValue(LANDMASK_TABLE).toString();
            std::string const landCoverTypeName = themeParameters->getValue(LAND_COVER_TYPE).toString();
            std::string const landAreaValue = themeParameters->getValue(TYPE_LAND_AREA).toString();
            std::string const clTableName = themeParameters->getValue(CL_TABLE).toString();
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

        ///
        double EdgeCleaningOp::_getLength( ign::geometry::Geometry const& geom ) const
        {
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
                        ign::geometry::LineString const& ls = geom.asLineString();
                        if (ls.isEmpty()) return 0;
                        double length = ls.length();
                        return length;
                    }
                    
                case ign::geometry::Geometry::GeometryTypeMultiLineString : 
                    {
                        double length = 0;
                        ign::geometry::MultiLineString const& mls = geom.asMultiLineString();
                        for( size_t i = 0 ; i < mls.numGeometries() ; ++i )
                            length += mls.lineStringN(i).length();
                        return length;
                    }
                
                case ign::geometry::Geometry::GeometryTypeGeometryCollection :
                    {
                        double length = 0;
                        ign::geometry::GeometryCollection const& collection = geom.asGeometryCollection();
                        for( size_t i = 0 ; i < collection.numGeometries() ; ++i )
                            length += _getLength( collection.geometryN(i) );
                    
                        return length;
                    }
                default :
                    IGN_THROW_EXCEPTION( "Geometry type not allowed" );
            }
        }

        ///
        ///
        ///
        void EdgeCleaningOp::_addLengths(std::string country, ign::geometry::LineString const& ls , double & lengthInCountry, double & length) const 
        {
            std::map<std::string, ign::geometry::GeometryPtr>::const_iterator mit = _mCountryGeomPtr.find(country);
            if (mit == _mCountryGeomPtr.end()) {
                _logger->log(epg::log::ERROR, "Unknown country [country code] " + country);
                return;
            }

            ign::geometry::GeometryPtr resultPtr(mit->second->Intersection(ls));

            double ratio = 0;

            lengthInCountry += _getLength(*resultPtr);
            length += ls.length();
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

            double ratio = 0;

            lengthInCountry += _getLength(*resultPtr);
            length += ls.length();
        }

        ///
        ///
        ///
        double EdgeCleaningOp::_getAntennaLength(GraphType const& graph, std::list<edge_descriptor> const& lEdges) const
        {
            double length = 0;
            std::list<edge_descriptor>::const_iterator lit = lEdges.begin();
            for ( ; lit != lEdges.end() ; ++lit) {
                ign::geometry::LineString edgeGeom = graph.getGeometry(*lit);
                length += edgeGeom.length();
            }
            return length;
        }

        ///
        ///
        ///
        double EdgeCleaningOp::_getRatio(GraphType const& graph, std::string country, std::list<oriented_edge_descriptor> const& path) const
        {
            std::list<edge_descriptor> lEdges;
            for (std::list<oriented_edge_descriptor>::const_iterator lit = path.begin() ; lit != path.end() ; ++lit)
                lEdges.push_back(lit->descriptor);

            return _getRatio(graph, country, lEdges);
        }

        ///
        ///
        ///
        double EdgeCleaningOp::_getRatio(GraphType const& graph, std::string country, std::list<edge_descriptor> const& lEdges) const
        {
            double lengthInCountry = 0;
            double length = 0;
            std::list<edge_descriptor>::const_iterator lit = lEdges.begin();
            for ( ; lit != lEdges.end() ; ++lit) {
                ign::geometry::LineString edgeGeom = graph.getGeometry(*lit);
                _addLengths(country, edgeGeom, lengthInCountry, length);
            }
			if (length == 0)
				return 0;

            return lengthInCountry / length;
        }

        ///
        ///
        ///
        double EdgeCleaningOp::_getRatioWithBuff(GraphType const& graph, std::string country, std::list<edge_descriptor> const& lEdges) const
        {
            double lengthInCountry = 0;
            double length = 0;
            std::list<edge_descriptor>::const_iterator lit = lEdges.begin();
            for ( ; lit != lEdges.end() ; ++lit) {
                ign::geometry::LineString edgeGeom = graph.getGeometry(*lit);
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

            std::list<edge_descriptor> lEdge2Remove;

            GraphType & graph = graphManager.getGraph();
            boost::progress_display display(graph.numFaces(), std::cout, "[ cleaning faces  % complete ]\n");
            face_iterator fit, fend;
			for( graph.faces( fit, fend ) ; fit != fend ; ++fit )
			{
                ++display;
				ign::geometry::Polygon faceGeom = graph.getGeometry( *fit );

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
                        _removeEdges(graph, mlEdges.rbegin()->second, lEdge2Remove);
                    } else if ( !hasConnection1 ) {
                        _removeEdges(graph, mlEdges.begin()->second, lEdge2Remove);
                    }
                }
			}
            for ( std::list<edge_descriptor>::const_iterator lit = lEdge2Remove.begin() ; lit != lEdge2Remove.end() ; ++lit )
                graph.removeEdge(*lit);
        }

        ///
        ///
        ///
        std::pair<ign::geometry::LineString,ign::geometry::LineString> EdgeCleaningOp::_getSubLineStrings(
            size_t id1, 
            size_t id2, 
            ign::geometry::LineString const& closedLs
        ) const {
            std::pair<ign::geometry::LineString, ign::geometry::LineString> pLs;
            size_t lastId = closedLs.numPoints()-1;
            if (id1 == id2 || id2 < id1 || id1 == 0 && id2 == lastId) {
                _logger->log(epg::log::ERROR, "Error in getSubLinestrings : bad indexes values");
                return pLs;
            }
                
            size_t id = id1;
            do {
                if (id == id2 || id == 0 && id2 == lastId) {
                    pLs.first.addPoint(closedLs.pointN(id));
                    pLs.second.addPoint(closedLs.pointN(id));
                } else if (id < id1)
                    pLs.second.addPoint(closedLs.pointN(id));
                else if (id < id2)
                    pLs.first.addPoint(closedLs.pointN(id));
                else
                    pLs.second.addPoint(closedLs.pointN(id));

                ++id;
                if (id == lastId) id = 0;
            } while (id != id1);

            pLs.second.addPoint(closedLs.pointN(id1));

            return pLs;
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_isSlimSurface( 
            ign::geometry::LineString const& closedLs, 
            double maxWidth
		) const {
            // trouver les 2 points les plus distants du contour
            // separer le contour en 2 polylignes
            // calculer hausdorff entre ces 2 polylignes
            int indexI = -1;
            int indexJ = -1;

            double maxDist = 0;
            for( size_t i = 0 ; i < closedLs.numPoints()-1 ; ++i )
            {
                for( size_t j = i+1 ; j < closedLs.numPoints() ; ++j )
                {
                    double dist = closedLs.pointN(i).distance(closedLs.pointN(j));
                    if ( dist > maxDist ) {
                        maxDist = dist;
                        indexI = i;
                        indexJ = j;
                    }
                }
            }

            std::pair<ign::geometry::LineString,ign::geometry::LineString> pLs = _getSubLineStrings(indexI, indexJ, closedLs);

            double hausdorffDist = ign::geometry::algorithm::HausdorffDistanceOp::distance(pLs.first, pLs.second);

            if (hausdorffDist < 0)
                _logger->log(epg::log::ERROR, "Error in slim surface calculation : hausdorff distance < 0");

            return hausdorffDist >= 0 && hausdorffDist < maxWidth;
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_isSlimSurface( 
            ign::geometry::Polygon const& poly, 
            double maxWidth
		) const {
            return _isSlimSurface( poly.exteriorRing(), maxWidth );
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

            _logger->log(epg::log::DEBUG, "cf4");

            do{
                // on ne doit pas avoir de cl dans la boucle
                if (graphManager.isCl(currentEdge.descriptor)) {
                    _logger->log(epg::log::WARN, "Loop contains a CL [cl id] "+graph.origins(currentEdge.descriptor)[0]);
                    return false;
                }

                vpCountryEdges.back().second.push_back(currentEdge);

                oriented_edge_descriptor nextEdge = ign::geometry::graph::detail::nextEdge( currentEdge, graph );
                std::string nextCountry = graphManager.getCountry(nextEdge.descriptor);

                if (graph.degree(graph.target(currentEdge)) > 2) {
                    if (nextEdge != startEdge)
                        vpCountryEdges.push_back(std::make_pair(nextCountry, std::list<oriented_edge_descriptor>()));
                } else if ( currentCountry != nextCountry ) {
                    _logger->log(epg::log::WARN, "Mixed country on path [edge id] "+graph.origins( currentEdge.descriptor)[0]);
                    return false;
                }
                currentCountry = nextCountry;
                currentEdge = nextEdge;
            }while( currentEdge != startEdge );

            if (vpCountryEdges.size() > 1 && vpCountryEdges.front().first == vpCountryEdges.back().first) {
                if ( graph.degree(graph.target(vpCountryEdges.back().second.back())) == 2 ){
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
            std::list<edge_descriptor>& lEdge2Remove
        ) const {
            std::list<edge_descriptor> lEdges;
            for (std::list<oriented_edge_descriptor>::const_iterator lit = path.begin() ; lit != path.end() ; ++lit)
                lEdges.push_back(lit->descriptor);

            _removeEdges(graph, lEdges, lEdge2Remove);
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
                ign::feature::FeatureFilter filter(countryCodeName +" = '"+vCountry[i]+"'");
                cleanFaces2(filter);
            }
        }

        ///
        ///
        ///
        void EdgeCleaningOp::cleanFacesAndAntennaByCountry() const
        {
            epg::Context* context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            bool isPlanar = true;

            std::vector<std::string> vCountry;
		    epg::tools::StringTools::Split(_countryCode, "#", vCountry);
            for (size_t i = 0 ; i < vCountry.size() ; ++i) {
                detail::EdgeCleaningGraphManager graphManager;
                ign::feature::FeatureFilter filter(countryCodeName +" LIKE '%"+vCountry[i]+"%'");
                _loadGraph(graphManager, isPlanar, filter);

                std::set<vertex_descriptor> sTreatedDangles;

                _cleanAntennas(graphManager, sTreatedDangles, isPlanar);
                bool bChangeOccured = _cleanFaces2(graphManager);

                while (bChangeOccured) {
                    bChangeOccured = _cleanAntennas(graphManager, sTreatedDangles, isPlanar);
                    if (bChangeOccured)
                        bChangeOccured = _cleanFaces2(graphManager);
                } 
            }
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::cleanFaces2(ign::feature::FeatureFilter filter) const
        {
            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager, true, filter);

            return _cleanFaces2(graphManager);
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_cleanFaces2(detail::EdgeCleaningGraphManager & graphManager) const
        {
            bool bChangeOccured = false;

            // app parameters
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const slimSurfaceWidth = themeParameters->getValue( ECL_SLIM_SURFACE_WIDTH ).toDouble();

            std::list<edge_descriptor> lEdge2Remove;

            GraphType & graph = graphManager.getGraph();
            boost::progress_display display(graph.numFaces(), std::cout, "[ cleaning faces 2  % complete ]\n");
            face_iterator fit, fend;
			for( graph.faces( fit, fend ) ; fit != fend ; ++fit )
			{
                ++display;

                _logger->log(epg::log::DEBUG, "cf1");

				ign::geometry::Polygon faceGeom = graph.getGeometry( *fit );

                _logger->log(epg::log::DEBUG, faceGeom.toString());
                _logger->log(epg::log::DEBUG, "cf2");

                //DEBUG
                // bool test = false;
                // if (faceGeom.distance(ign::geometry::Point(3823437.97,3093497.89)) < 1) {
                //     test = true;
                // }
                // // DEBUG
                // _logger->log(epg::log::DEBUG, faceGeom.toString());
                // if (faceGeom.intersects(ign::geometry::Point(3823392.232,3093616.929))) {
                //     bool test = true;
                // }

				if (_isSlimSurface(faceGeom, slimSurfaceWidth)) {

                    _logger->log(epg::log::DEBUG, "cf3");

                    ign::feature::Feature feat;
                    feat.setGeometry(faceGeom);
                    _shapeLogger->writeFeature("ecl_slim_surface", feat);

                    std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> vpCountryEdges;
                    if (!_getFacePaths(graphManager, *fit, vpCountryEdges))
                        continue;

                    _logger->log(epg::log::DEBUG, "cf5");
                    
                    if (vpCountryEdges.size() == 1) {
                        ign::feature::Feature feat;
                        feat.setGeometry(faceGeom);
                        _shapeLogger->writeFeature("ecl_slim_face_1_path", feat);

                        _removePath(graph, vpCountryEdges.front().second, lEdge2Remove);
                        bChangeOccured = true;
                        continue;
                    }

                    if (vpCountryEdges.size() != 2) continue;

                    _logger->log(epg::log::DEBUG, "cf6");

                    if (vpCountryEdges.front().first != vpCountryEdges.back().first) {
                        // quel chemin doit-on garder ?
                        double ratio1 = _getRatio(graph, vpCountryEdges.front().first, vpCountryEdges.front().second);
                        bool hasConnection1 = false; /*todo*/
                        double ratio2 = _getRatio(graph, vpCountryEdges.back().first, vpCountryEdges.back().second);
                        bool hasConnection2 = false; /*todo*/

                        _logger->log(epg::log::DEBUG, "cf7");

                        if ( ratio1 > ratio2 ) {
                            if (!hasConnection2 ) {
                                _logger->log(epg::log::DEBUG, "cf8");
                                _removePath(graph, vpCountryEdges.back().second, lEdge2Remove);
                                bChangeOccured = true;
                            }
                        } else if ( !hasConnection1 ) {
                            _logger->log(epg::log::DEBUG, "cf9");
                            _removePath(graph, vpCountryEdges.front().second, lEdge2Remove);
                            bChangeOccured = true;
                        }
                    } else {
                        // 2 chemins dans le même pays, on garde le plus court
                        double lengthFront = _getPathLength(graph, vpCountryEdges.front().second);
                        double lengthBack = _getPathLength(graph, vpCountryEdges.back().second);
                        if (lengthFront < lengthBack) {
                            _removePath(graph, vpCountryEdges.back().second, lEdge2Remove);
                            bChangeOccured = true;
                        } else {
                            _removePath(graph, vpCountryEdges.front().second, lEdge2Remove);
                            bChangeOccured = true;
                        }
                        ign::feature::Feature feat;
                        feat.setGeometry(faceGeom);
                        _shapeLogger->writeFeature("ecl_slim_face_2_path_same_country", feat);
                    }
                } else {
                    // grande face
                    // si dans mauvais pays on supprime
                    _logger->log(epg::log::DEBUG, "gf1");

                    std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> vpCountryEdges;
                    if (!_getFacePaths(graphManager, *fit, vpCountryEdges))
                        continue;

                    _logger->log(epg::log::DEBUG, "gf2");
                    
                    if (vpCountryEdges.size() == 1) {
                        double ratio = _getRatio(graph, vpCountryEdges.front().first, vpCountryEdges.front().second);

                        if(ratio == 0) {
                            ign::feature::Feature feat;
                            feat.setGeometry(faceGeom);
                            _shapeLogger->writeFeature("ecl_big_face_removed", feat);

                            _logger->log(epg::log::DEBUG, "gf3");
                            _removePath(graph, vpCountryEdges.front().second, lEdge2Remove);
                            bChangeOccured = true;
                        }
                    }
                    _logger->log(epg::log::DEBUG, "gf4");
                }
			}
            _logger->log(epg::log::DEBUG, "cf10");
            for ( std::list<edge_descriptor>::const_iterator lit = lEdge2Remove.begin() ; lit != lEdge2Remove.end() ; ++lit )
                graph.removeEdge(*lit);
            _logger->log(epg::log::DEBUG, "cf11");

            return bChangeOccured;
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

                _loadGraph(graphManager, false, filter);
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
            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager, false);

            std::set<vertex_descriptor> sTreatedDangles;
            return _cleanAntennas(graphManager, sTreatedDangles);
        }

        ///
        ///
        ///
        bool EdgeCleaningOp::_cleanAntennas(detail::EdgeCleaningGraphManager & graphManager, std::set<vertex_descriptor> & sTreatedDangles, bool isPlanarGraph) const
        {
            GraphType & graph = graphManager.getGraph();
            
            size_t oldNumTreatedDangles = sTreatedDangles.size();
            std::set< vertex_descriptor > visitedVertices;

            boost::progress_display display2(graph.numVertices(), std::cout, "[ cleaning antennas  % complete ]\n");
            vertex_iterator vit, vend;
            for( graph.vertices( vit, vend ) ; vit != vend ; ++vit )
            {
                ++display2;

                if( sTreatedDangles.find( *vit ) != sTreatedDangles.end() ) continue;
                if( visitedVertices.find( *vit ) != visitedVertices.end() ) continue;
                if( graph.degree( *vit ) != 1 ) continue;

                sTreatedDangles.insert(*vit);

                // seulement les vertex qui touchent une CL ?

                // DEBUG
                ign::geometry::Point p = graph.getGeometry(*vit);
                _logger->log(epg::log::DEBUG, p.toString());

                std::list<edge_descriptor> lAntennaEdges;
                bool isConnected2CF = false;

                std::vector< oriented_edge_descriptor > vEdges;
                graph.incidentEdges( *vit, vEdges );

                oriented_edge_descriptor nextEdge = vEdges.front(); // si nextEdge n'est pas une CL ?
                if (graphManager.isCl(nextEdge.descriptor)) {
                    _logger->log(epg::log::WARN, "Antenna is connecting line [cl id] "+graph.origins(nextEdge.descriptor)[0]);
                    continue;
                }
                std::string country = graphManager.getCountry(nextEdge.descriptor);

                //DEBUG
                // if (country=="be#fr") {
                //     bool test = true;
                // }
                // std::string id = graph.origins(nextEdge.descriptor)[0];
                // if (id == "a15662a7-be3c-434f-8124-2a3ab80691d8") {
                //     bool test = true;
                // }

                vertex_descriptor vTarget = GraphType::nullVertex();
                
                while( true )
                {
                    if (graphManager.getCountry(nextEdge.descriptor) != country) {
                        isConnected2CF = true;
                        break;
                    }

                    _addAntennaEdges(graph, nextEdge.descriptor, lAntennaEdges, isPlanarGraph);

                    vTarget = _getTarget(graph, nextEdge, isPlanarGraph);

                    if (graphManager.isCl(nextEdge.descriptor) && graph.degree(graph.source( nextEdge )) == 2 ) {
                        _logger->log(epg::log::WARN, "Antenna connected to connecting line [cl id] "+graph.origins(nextEdge.descriptor)[0]);
                        isConnected2CF = true;
                        break;
                    }

                    bool targetIsCp = _vertexIsCp(graph, vTarget);

                    if( graph.degree( vTarget ) != 2 || targetIsCp) { // ou si nextEdge est une CL ?
                        if (_vertexIsConnected2Cl(graphManager, vTarget) || targetIsCp) {
                            isConnected2CF = true;
                        }
                        if( graph.degree( vTarget ) == 1 /*antenne isolee*/)
                        {
                            visitedVertices.insert( vTarget );
                        }
                        break;
                    }

                    nextEdge = _getNextEdge(graph, nextEdge.descriptor, vTarget, isPlanarGraph);
                }
                // DEBUG
                _logger->log(epg::log::DEBUG, "ca1");

                if (!lAntennaEdges.empty()) {
                    _cleanAntenna(graph, country, lAntennaEdges, isConnected2CF, isPlanarGraph);
                }

                // DEBUG
                _logger->log(epg::log::DEBUG, "ca2");
            }
            sTreatedDangles.insert(visitedVertices.begin(), visitedVertices.end());
            return oldNumTreatedDangles != sTreatedDangles.size();
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
        void EdgeCleaningOp::_addAntennaEdges(GraphType const& graph, edge_descriptor e, std::list<edge_descriptor> & lEdges, bool isPlanarGraph) const {
            if (isPlanarGraph) {
                std::string origin = graph.origins(e)[0];
                std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdges = graph.getInducedEdges(origin);
                if (foundInducedEdges.first) {
                    for ( std::vector<oriented_edge_descriptor>::const_iterator vit = foundInducedEdges.second.begin() ; vit != foundInducedEdges.second.end() ; ++vit ) {
                        lEdges.push_back(vit->descriptor);
                    }
                }
            } else {
                lEdges.push_back(e);
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
        void EdgeCleaningOp::_cleanAntenna(
            GraphType & graph,
            std::string const& country,
            std::list<edge_descriptor> const& lAntennaEdges,
            bool bAntennaIsConnected2CF,
            bool isPlanarGraph
        ) const {
            // app parameters
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const antennaRatioThreshold = themeParameters->getValue( ECL_ANTENNA_RATIO_THRESHOLD ).toDouble();
            double const antennaRatioThresholdWithBuff = themeParameters->getValue( ECL_ANTENNA_RATIO_WITH_BUFFER_THRESHOLD ).toDouble();
            double const antennaMinLength = themeParameters->getValue( ECL_ANTENNA_MIN_LENGTH ).toDouble();
            
            double antennaLength = _getAntennaLength(graph, lAntennaEdges);
            if (bAntennaIsConnected2CF && antennaLength < antennaMinLength) {
                _removeEdgesAndGraphEdges(graph, lAntennaEdges, isPlanarGraph);
                return;
            }
            
            //DEBUG
            _logger->log(epg::log::DEBUG, "Antenna info :");
            std::list<edge_descriptor>::const_iterator lit;
            std::set<std::string> sOrigins;
            for (lit = lAntennaEdges.begin() ; lit != lAntennaEdges.end() ; ++lit) {
                for (size_t j = 0 ; j < graph.origins(*lit).size() ; ++j) {
                    sOrigins.insert( graph.origins(*lit)[j] );
                }
            }
            _logger->log(epg::log::DEBUG, tools::StringTools::toString(sOrigins));



            double ratio = _getRatio(graph, country, lAntennaEdges);
            //DEBUG
            _logger->log(epg::log::DEBUG, std::to_string(ratio));

            if (ratio < antennaRatioThreshold) {
                _removeEdgesAndGraphEdges(graph, lAntennaEdges, isPlanarGraph);
            } 
            else {
                double ratioWithBuff = _getRatioWithBuff(graph, country, lAntennaEdges);

                if (ratioWithBuff < antennaRatioThresholdWithBuff) {
                    _removeEdgesAndGraphEdges(graph, lAntennaEdges, isPlanarGraph);
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

            std::list<edge_descriptor> lEdge2Remove;

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

                ign::geometry::LineString lsRef = graph.getGeometry(*vParallelEdges.begin());
                std::set<edge_descriptor> sVisitedTemp;
                sVisitedTemp.insert(vParallelEdges.begin()->descriptor);

                std::vector< oriented_edge_descriptor >::const_iterator vit;
                for (vit = vParallelEdges.begin() ; vit != vParallelEdges.end() ; ++vit ) {
                    if ( sVisitedTemp.find(vit->descriptor) != sVisitedTemp.end() ) continue; /*gestion des boucles + sert a passer le 1er element*/
                    sVisitedTemp.insert(vit->descriptor);

                    ign::geometry::LineString ls = graph.getGeometry(vit->descriptor);
                    if ( ign::geometry::algorithm::HausdorffDistanceOp::distance(lsRef, ls) < maxDist ) {
                        _removeEdges(graph, std::list<edge_descriptor>(1,vit->descriptor), lEdge2Remove);
                    }
                }
                sVisitedEdge.insert(sVisitedTemp.begin(), sVisitedTemp.end());
            }

            for ( std::list<edge_descriptor>::const_iterator lit = lEdge2Remove.begin() ; lit != lEdge2Remove.end() ; ++lit )
                graph.removeEdge(*lit);
        }

        ///
        ///
        ///
        void EdgeCleaningOp::cleanAll() const
        {
            cleanFaces();
            cleanPathsOutOfCountry();
            cleanParalelleEdges();
            cleanFacesAndAntennaByCountry();
        };
    }
}