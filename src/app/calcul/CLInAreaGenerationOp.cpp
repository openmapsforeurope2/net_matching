#include <app/calcul/CLInAreaGenerationOp.h>

//APP
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LengthIndexedLineString.h>
#include <app/calcul/detail/graph/concept/EdgeCleaningGraphSpecializations.h>

//BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

//EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/sql/tools/numFeatures.h>
#include <epg/tools/geometry/project.h>
#include <epg/graph/tools/reverse.h>

//SOCLE
#include <ign/geometry/algorithm/HausdorffDistanceOp.h>
#include <ign/geometry/algorithm/AreaOp.h>
#include <ign/geometry/graph/detail/NextEdge.h>


namespace app
{
    namespace calcul
    {

		///
        ///
        ///
        void CLInAreaGenerationOp::compute(
            bool verbose
        ) {
            CLInAreaGenerationOp cLInAreaGenerationOp(verbose);
            cLInAreaGenerationOp._compute();
        }

		///
        ///
        ///
        CLInAreaGenerationOp::CLInAreaGenerationOp(
            bool verbose
        ) : _verbose(verbose)
        {
            _init();
        }

		///
		///
		///
		CLInAreaGenerationOp::~CLInAreaGenerationOp()
		{
			// _shapeLogger->closeShape("CLBeforeMerge");
		}

		///
        ///
        ///
        void CLInAreaGenerationOp::_init()
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

            //--
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();

			//--
            _fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);

            //--
            // _shapeLogger = epg::log::ShapeLoggerS::getInstance();
            // _shapeLogger->addShape("ec_projected_antennas", epg::log::ShapeLogger::LINESTRING);

            //--
            std::string listAttrWName = themeParameters->getValue(LIST_ATTR_W).toString();
            std::string listAttrJsonName = themeParameters->getValue(LIST_ATTR_JSON).toString();
            _attrMerger.setLists( listAttrWName, listAttrJsonName, "/");

            //--
            _logger->log(epg::log::INFO, "[END] initialization: " + epg::tools::TimeTools::getTime());
		}

        ///
        ///
        ///
        void CLInAreaGenerationOp::_createCLOnOverlappingEdges(
            GraphType const& graph,
            std::map<std::string, std::set<edge_descriptor>> & mFeatMergedEdges,
            std::multimap<std::string, detail::IncidentFeature> & mmIncidentFeatures
        ) const {
            //--
			epg::Context *context = epg::ContextS::getInstance();

			//--
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            //--
            app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG).getValue().toString();

            boost::progress_display display(graph.numEdges(), std::cout, "[ generating CL on overlapping edges % complete ]\n");

            edge_iterator eit, eend;
            for (graph.edges(eit, eend); eit != eend; ++eit)
            {
                ++display;

                std::vector<std::string> vOrigins = graph.origins(*eit);
                if( vOrigins.size() != 2 ) continue;

                //--
                ign::geometry::LineString edgeGeom = graph.getGeometry( *eit );

                //--
                ign::feature::Feature featFront;
                _fsEdge->getFeatureById(vOrigins.front(), featFront);
                ign::feature::Feature featBack;
                _fsEdge->getFeatureById(vOrigins.back(), featBack);

                //--
                std::map<std::string, std::set<edge_descriptor>>::iterator mit;
                mit = mFeatMergedEdges.find(featFront.getId());
                if( mit == mFeatMergedEdges.end() )
                    mit = mFeatMergedEdges.insert(std::make_pair(featFront.getId(), std::set<edge_descriptor>())).first;
                mit->second.insert(*eit);

                mit = mFeatMergedEdges.find(featBack.getId());
                if( mit == mFeatMergedEdges.end() )
                    mit = mFeatMergedEdges.insert(std::make_pair(featBack.getId(), std::set<edge_descriptor>())).first;
                mit->second.insert(*eit);


                // pour merger les attributs dans le bon sens
                ign::feature::Feature* fRef;
                ign::feature::Feature* f2Merge;
                fRef = featFront.getAttribute(countryCodeName).toString() < featBack.getAttribute(countryCodeName).toString() ? &featFront : &featBack;
                f2Merge = featFront.getAttribute(countryCodeName).toString() < featBack.getAttribute(countryCodeName).toString() ? &featBack : &featFront;


                //--
                _attrMerger.mergeFeatAttribute( *fRef, *f2Merge, "#" );
                fRef->setAttribute(wTagName, ign::data::String(fRef->getId()+"#"+f2Merge->getId()));
 

                // assurer la cohérence avec les données sources aux extrémités et features incidents
                std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdgesFront = graph.getInducedEdges(vOrigins.front());
                std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdgesBack = graph.getInducedEdges(vOrigins.back());
                bool isStartingPathFront = *eit == foundInducedEdgesFront.second.begin()->descriptor;
                bool isEndingPathFront = *eit == foundInducedEdgesFront.second.rbegin()->descriptor;
                bool isStartingPathBack = *eit == foundInducedEdgesBack.second.begin()->descriptor;
                bool isEndingPathBack = *eit == foundInducedEdgesBack.second.rbegin()->descriptor;
                if( (isStartingPathFront || isEndingPathFront) ) {
                    if (isStartingPathBack || isEndingPathBack) {
                        ign::geometry::LineString const& featFrontGeom = featFront.getGeometry().asLineString();
                        if (isStartingPathFront)
                            if( foundInducedEdgesFront.second.begin()->direction == ign::graph::DIRECT)
                                edgeGeom.startPoint() = featFrontGeom.startPoint();
                            else
                                edgeGeom.endPoint() = featFrontGeom.startPoint();
                        if (isEndingPathFront)
                            if( foundInducedEdgesFront.second.rbegin()->direction == ign::graph::DIRECT)
                                edgeGeom.endPoint() = featFrontGeom.endPoint();
                            else
                                edgeGeom.startPoint() = featFrontGeom.endPoint();
                    } else {
                        std::set<edge_descriptor> sEdgesBack;
                        for( std::vector<oriented_edge_descriptor>::const_iterator vit = foundInducedEdgesBack.second.begin() ; vit != foundInducedEdgesBack.second.end() ; ++vit )
                            sEdgesBack.insert(vit->descriptor);

                        if( isStartingPathFront )
                            _addIncidentFeatures(graph, graph.source(*foundInducedEdgesFront.second.begin()), sEdgesBack, mmIncidentFeatures);
                        if( isEndingPathFront )
                            _addIncidentFeatures(graph, graph.target(*foundInducedEdgesFront.second.rbegin()), sEdgesBack, mmIncidentFeatures);
                    }
                } else if (isStartingPathBack || isEndingPathBack) {
                    std::set<edge_descriptor> sEdgesFront;
                    for( std::vector<oriented_edge_descriptor>::const_iterator vit = foundInducedEdgesFront.second.begin() ; vit != foundInducedEdgesFront.second.end() ; ++vit )
                            sEdgesFront.insert(vit->descriptor);

                    if( isStartingPathBack )
                        _addIncidentFeatures(graph, graph.source(*foundInducedEdgesBack.second.begin()), sEdgesFront, mmIncidentFeatures);
                    if( isEndingPathBack )
                        _addIncidentFeatures(graph, graph.target(*foundInducedEdgesBack.second.rbegin()), sEdgesFront, mmIncidentFeatures);
                }

                fRef->setGeometry(edgeGeom);
                
				_fsEdge->createFeature(*fRef);
            }
        }

        ///
        ///
        ///
        void CLInAreaGenerationOp::_addIncidentFeatures(
            GraphType const& graph,
            vertex_descriptor v,
            std::set<edge_descriptor> const& sEdges,
            std::multimap<std::string, detail::IncidentFeature> & mmIncidentFeatures
        ) const {
            if( graph.degree(v) > 2 ) {
                ign::geometry::Point vertexGeom = graph.getGeometry(v);
                std::vector< edge_descriptor > vEdges = graph.incidentEdges(v);

                for (size_t i = 0 ; i < vEdges.size() ; ++i) {
                    if (sEdges.find(vEdges[i]) != sEdges.end()) continue;

                    std::string incidentFeatId = graph.origins(vEdges[i])[0];

                    ign::feature::Feature incidentFeat;
                    _fsEdge->getFeatureById(incidentFeatId, incidentFeat);
                    ign::geometry::LineString const& incidentFeatGeom = incidentFeat.getGeometry().asLineString();

                    bool touchStart = incidentFeatGeom.startPoint().distance(vertexGeom) < 1e-5;
                    bool touchEnd = incidentFeatGeom.endPoint().distance(vertexGeom) < 1e-5;

                    if (touchStart)
                        mmIncidentFeatures.insert(std::make_pair(incidentFeatId, detail::IncidentFeature(incidentFeatId, detail::START)));
                    if (touchEnd)
                        mmIncidentFeatures.insert(std::make_pair(incidentFeatId, detail::IncidentFeature(incidentFeatId, detail::END)));

                    if (!touchStart && !touchEnd) {
                        _logger->log(epg::log::ERROR, "Error in incident feature computation : "+vertexGeom.toString());
                    }

                }
            }
        }
        

		///
        ///
        ///
        void CLInAreaGenerationOp::_compute() const
        {
			//--
			epg::Context *context = epg::ContextS::getInstance();

			//--
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const geomName = epgParams.getValue(GEOM).toString();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

			//--
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG).getValue().toString();
            double const slimSurfaceWidth = themeParameters->getValue( CLA_SURFACE_WIDTH ).toDouble();

            //--
            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager);
            GraphType const& graph = graphManager.getGraph();

            //--
            std::multimap<std::string, detail::IncidentFeature> mmIncidentFeatures;
            std::map<std::string, std::set<edge_descriptor>> mFeatMergedEdges;

            //--
            _createCLOnOverlappingEdges(graph, mFeatMergedEdges, mmIncidentFeatures);

            //TODO :
            // créer les CL pour les edges superposés
            // ajouter les patiences manquantes
            // fusionner les cl (EDGE_LINK=w_national_identifier)

            //WARNING: on casse l'unicité des liens w_national_identifier  !!!!

            boost::progress_display display(graph.numFaces(), std::cout, "[ generating CL in area  % complete ]\n");

            

            face_iterator fit, fend;
            for( graph.faces( fit, fend ) ; fit != fend ; ++fit )
			{
                ++display;

				ign::geometry::Polygon faceGeom = graph.getGeometry( *fit );

                //DEBUG
                // if( faceGeom.intersects(ign::geometry::Point(4031796.58,2934527.57))) {
                //     bool test = true;
                // } 
                // if( faceGeom.intersects(ign::geometry::Point(3986976.57,2956178.37))) {
                //     bool test = true;
                // }

				std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> vpCountryEdges;
				if (!_getFacePaths(graphManager, *fit, vpCountryEdges))
					continue;

				if (vpCountryEdges.size() < 2) 
					continue;

				std::set<std::string> sFaceCountries;
				for (std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>>::const_iterator vpit = vpCountryEdges.begin() ; vpit != vpCountryEdges.end() ; ++vpit)
					sFaceCountries.insert(vpit->first);

				if (sFaceCountries.size() != 2) 
					continue;

				std::set<std::string> hasConnection = _mergeFacePaths(vpCountryEdges);

				// if (!hasConnection.empty())
				// 	continue;
				
				if (vpCountryEdges.size() != 2)
					continue;

                // mettre les 2 chemins dans le meme sens
                if (graph.source(*vpCountryEdges.front().second.begin()) != graph.source(*vpCountryEdges.back().second.begin()))
                    vpCountryEdges.back().second = _getReversePath(vpCountryEdges.back().second);

				ign::geometry::LineString lsFront = _convertPathToLineString(graph, vpCountryEdges.front().first, vpCountryEdges.front().second);
                ign::geometry::LineString lsBack = _convertPathToLineString(graph, vpCountryEdges.back().first, vpCountryEdges.back().second);

                if( _pathsGeomAreEqual(faceGeom, lsFront, lsBack, slimSurfaceWidth) ) {
					ign::geometry::LineString meanGeom = _computeMeanPath(lsFront, lsBack);
					geometry::tools::LengthIndexedLineString lsIndexMean(meanGeom);
                    double meanLength = meanGeom.length();

					//todo recupérer tous les edges des chemins
                    std::map<double, std::vector<detail::IncidentFeature>> mAbsIncidentFeatures;
					std::map<double, std::string> sAbsEdgeFront = _getOriginEdges(graph, vpCountryEdges.front().first, vpCountryEdges.front().second, mFeatMergedEdges, mAbsIncidentFeatures); 
					std::map<double, std::string> sAbsEdgeBack = _getOriginEdges(graph, vpCountryEdges.back().first,vpCountryEdges.back().second, mFeatMergedEdges, mAbsIncidentFeatures);

					double sStart = 0;
					do {
						double absFront = sAbsEdgeFront.begin()->first;
						double absBack = sAbsEdgeBack.begin()->first;
						ign::feature::Feature featFront;
                        _fsEdge->getFeatureById(sAbsEdgeFront.begin()->second, featFront);
						ign::feature::Feature featBack;
                        _fsEdge->getFeatureById(sAbsEdgeBack.begin()->second, featBack);

						ign::geometry::LineString newFeatGeom;
						double sTarget;
						if( absFront < absBack ) {
							sAbsEdgeFront.erase(sAbsEdgeFront.begin());
							sTarget = absFront;
						} else if ( absFront > absBack ) {
							sAbsEdgeBack.erase(sAbsEdgeBack.begin());
							sTarget = absBack;
						} else {
                            sAbsEdgeFront.erase(sAbsEdgeFront.begin());
                            sAbsEdgeBack.erase(sAbsEdgeBack.begin());
                            sTarget = absFront;
                        }
						ign::geometry::LineString lsNew = lsIndexMean.getSubLineString(sStart*meanLength, sTarget*meanLength);

                        //--
                        if (sStart != 0) {
                            std::map<double, std::vector<detail::IncidentFeature>>::iterator mit = mAbsIncidentFeatures.find(sStart);
                            if ( mit != mAbsIncidentFeatures.end() ) {
                                for( std::vector<detail::IncidentFeature>::iterator vit = mit->second.begin() ; vit != mit->second.end() ; ++vit ) {
                                    vit->ptTarget = lsNew.startPoint();
                                    mmIncidentFeatures.insert(std::make_pair(vit->originId, *vit));
                                }
                            }
                        }

                        // pour merger les attributs dans le bon sens
                        ign::feature::Feature* fRef;
                        ign::feature::Feature* f2Merge;
                        fRef = featFront.getAttribute(countryCodeName).toString() < featBack.getAttribute(countryCodeName).toString() ? &featFront : &featBack;
                        f2Merge = featFront.getAttribute(countryCodeName).toString() < featBack.getAttribute(countryCodeName).toString() ? &featBack : &featFront;

                        //--
						_attrMerger.mergeFeatAttribute( *fRef, *f2Merge, "#" );
                        fRef->setAttribute(wTagName, ign::data::String(fRef->getId()+"#"+f2Merge->getId()));

                        //PATCH
                        // lsNew.setFillZ(0);
                        if( lsNew.isEmpty() || lsNew.isNull() || !lsNew.isValid()) {
                            _logger->log(epg::log::ERROR, "Resulting inconsistent geometry [face geom]: "+faceGeom.toString());
                            continue;
                        }

						fRef->setGeometry(lsNew);

						_fsEdge->createFeature(*fRef);
						sStart = sTarget;
					} while(sAbsEdgeFront.size() > 0 && sAbsEdgeBack.size() > 0);
				}
			}

            //DEBUG
            _logger->log(epg::log::DEBUG, "yeap1");

            // on parcours tous les features qui ont des edges induits mergés : on les split si nécessaire
            boost::progress_display display2(mFeatMergedEdges.size(), std::cout, "[ splitting features  % complete ]\n");

            std::set<std::string> sIncident2delete;
            std::vector<detail::IncidentFeature> vIncident2Create;
            for( std::map<std::string, std::set<edge_descriptor>>::const_iterator mit = mFeatMergedEdges.begin() ; mit != mFeatMergedEdges.end() ; ++mit ) {
                ++display2;
                //DEBUG
                _logger->log(epg::log::DEBUG, "youp1");
                _logger->log(epg::log::DEBUG, mit->first);
                // if (mit->first == "c1e3d360-13e4-4ba5-a59f-c15a919fb883") {
                //     bool test = true;
                // }
                // if (mit->first == "03e700dc-77ed-4569-a74c-506641f8a07c") {
                //     bool test = true;
                // }

                ign::feature::Feature featOrigin;
                _fsEdge->getFeatureById(mit->first, featOrigin);
                ign::geometry::LineString geomOrigin = featOrigin.getGeometry().asLineString();
                std::string originCountry = featOrigin.getAttribute(countryCodeName).toString();
                std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdges = graph.getInducedEdges(mit->first);
                std::vector<std::list<oriented_edge_descriptor>> vPaths;
                std::list<oriented_edge_descriptor> path;
			    for( std::vector<oriented_edge_descriptor>::const_iterator vit = foundInducedEdges.second.begin() ; vit != foundInducedEdges.second.end() ; ++vit ) {
                    if (mit->second.find(vit->descriptor) == mit->second.end() ) {
                        path.push_back(*vit);
                    } else {
                        if (!path.empty()) {
                            vPaths.push_back(path);
                            path.clear();
                        }
                    }
                }
                if (!path.empty()) vPaths.push_back(path);

                //DEBUG
                _logger->log(epg::log::DEBUG, "youp2");

                //gestion incident edges
                if (mmIncidentFeatures.find(mit->first) != mmIncidentFeatures.end())
                    sIncident2delete.insert(mit->first);

                for( size_t i = 0 ; i < vPaths.size() ; ++i ) {
                    ign::geometry::LineString newGeom = _convertPathToLineString(graph, originCountry, vPaths[i]);

                    // assurer la cohérence avec les données sources aux extrémités
                    bool isStartingPath = graph.source(*vPaths[i].begin()) ==  graph.source(*foundInducedEdges.second.begin());
                    bool isEndingPath = graph.target(*vPaths[i].rbegin()) ==  graph.target(*foundInducedEdges.second.rbegin());
                    if (isStartingPath)
                        newGeom.startPoint() = geomOrigin.startPoint();
                    if (isEndingPath)
                        newGeom.endPoint() = geomOrigin.endPoint();
                    
                    //PATCH
                    // newGeom.setFillZ(0);

                    //DEBUG
                    _logger->log(epg::log::DEBUG, "youp3");

                    if( newGeom.isEmpty() || newGeom.isNull() || !newGeom.isValid()) {
                        _logger->log(epg::log::ERROR, "Resulting inconsistent geometry [orign feature id]: "+mit->first);
                        continue;
                    }

                    //DEBUG
                    _logger->log(epg::log::DEBUG, "youp4");

                    featOrigin.setGeometry(newGeom);
                    _fsEdge->createFeature(featOrigin);

                    //DEBUG
                    _logger->log(epg::log::DEBUG, "youp5");

                    //gestion incident edges
                    if ( i == 0 || i == vPaths.size()-1 ) {
                        std::pair < m_iterator, m_iterator > bounds = mmIncidentFeatures.equal_range( mit->first );
                        for ( m_iterator it = bounds.first ; it != bounds.second ; ++it )
                        {
                            if( (isStartingPath && it->second.ending == detail::START)
                                ||
                                (isEndingPath && it->second.ending == detail::END)
                            ) {
                                vIncident2Create.push_back(it->second);
                                vIncident2Create.back().originId = featOrigin.getId();
                            }
                        }
                    }

                    //DEBUG
                    _logger->log(epg::log::DEBUG, "youp6");
                }
                //DEBUG
                _logger->log(epg::log::DEBUG, "youp7");

                _fsEdge->deleteFeature(mit->first);

                //DEBUG
                _logger->log(epg::log::DEBUG, "youp8");
            }

            //DEBUG
            _logger->log(epg::log::DEBUG, "yeap2");

            //maj incident edges
            for( std::set<std::string>::const_iterator sit = sIncident2delete.begin() ; sit != sIncident2delete.end() ; ++sit ) {
                mmIncidentFeatures.erase(*sit);
            }
            //DEBUG
            _logger->log(epg::log::DEBUG, "yeap3");
            for( std::vector<detail::IncidentFeature>::const_iterator vit = vIncident2Create.begin() ; vit != vIncident2Create.end() ; ++vit ) {
                mmIncidentFeatures.insert(std::make_pair(vit->originId, *vit));
            }

            //DEBUG
            _logger->log(epg::log::DEBUG, "yeap4");

            //deplacement incident edges
            boost::progress_display display3(mmIncidentFeatures.size(), std::cout, "[ incident edges displacement  % complete ]\n");
            for( std::multimap<std::string, detail::IncidentFeature>::const_iterator mit = mmIncidentFeatures.begin() ; mit != mmIncidentFeatures.end() ; ++mit ) {
                ++display3;

                ign::feature::Feature feat;
                _fsEdge->getFeatureById(mit->first, feat);
                ign::geometry::LineString featGeom = feat.getGeometry().asLineString();

                if( mit->second.ending == detail::START )
                    featGeom.startPoint() = mit->second.ptTarget;
                else
                    featGeom.endPoint() = mit->second.ptTarget;

                if( featGeom.isEmpty() || featGeom.isNull() || !featGeom.isValid()) {
                    _logger->log(epg::log::ERROR, "Incident feature inconsistent geometry [id]: "+mit->first);
                    continue;
                }

                feat.setGeometry(featGeom);
                _fsEdge->modifyFeature(feat);
            }

            //DEBUG
            _logger->log(epg::log::DEBUG, "yeap5");

            _mergeByWTag();
		}

        ///
        ///
        ///
        void CLInAreaGenerationOp::_mergeByWTag() const {
            //--
            detail::EdgeCleaningGraphManager graphManager;
            ign::feature::FeatureFilter filter;
            _loadGraph(graphManager, filter, false);
            GraphType const& graph = graphManager.getGraph();

            //patience
            boost::progress_display display( graph.numEdges() , std::cout, "[ merging CL % complete ]\n") ;

            std::set< edge_descriptor > sMergedEdges;

            edge_iterator eit, eend ;
            for( graph.edges( eit, eend ) ; eit != eend ; ++eit )
            {
                ++display;
                if( !graphManager.isCl(*eit)) continue;
                if( graph.target(*eit) == graph.source(*eit) ) continue;
                if( sMergedEdges.find(*eit) != sMergedEdges.end() ) continue ;

                edge_descriptor edgePivot = *eit ;

                oriented_edge_descriptor tPivot[] = { 
                    oriented_edge_descriptor( *eit, ign::graph::DIRECT ), 
                    oriented_edge_descriptor( *eit, ign::graph::REVERSE ) 
                } ;

                //liste des arcs a fusionner
                edges_path path;
                path.push_back(tPivot[0]);
                
                for( size_t  i = 0 ; i < 2 ; ++i)
                {
                    oriented_edge_descriptor nextEdge = tPivot[i] ;

                    vertex_descriptor vTarget = graph.target( nextEdge );
                    vertex_descriptor vStart = graph.source( nextEdge );

                    if( graph.degree( vTarget ) != 2 ) continue ;

                    bool isLoop = false;
                            
                    while( true ) {

                        if( vTarget == vStart ) {
                            isLoop = true;
                            break;
                        }

                        std::vector< oriented_edge_descriptor > vIncidentEdges ;
                        graph.incidentEdges( vTarget, vIncidentEdges);
                        nextEdge = (vIncidentEdges.front().descriptor == nextEdge.descriptor )? vIncidentEdges.back() : vIncidentEdges.front() ;

                        if( i == 0 )
                            path.push_back( nextEdge );
                        else
                            path.push_front( epg::graph::tools::reverse( nextEdge ) );

					    vTarget = graph.target( nextEdge );

					    if( graph.degree( vTarget ) != 2 ) break;
                    }
                    if (isLoop) break;
                }

                if (path.size() > 1) {
                    ign::geometry::LineString pathGeom = _convertPathToLineString(graph, "", path);

                    ign::feature::Feature featRef;
                    _fsEdge->getFeatureById(graph.origins(path.begin()->descriptor)[0], featRef);
                    featRef.setGeometry(pathGeom);
                    _fsEdge->createFeature(featRef);

                    for (edges_path_const_iterator pit = path.begin() ; pit != path.end() ; ++pit)
                        sMergedEdges.insert(pit->descriptor);
                        
                }
            }
            for( std::set<edge_descriptor>::const_iterator sit = sMergedEdges.begin() ; sit != sMergedEdges.end() ; ++sit ) {
                _fsEdge->deleteFeature(graph.origins(*sit)[0]);
            }
        }

        ///
        ///
        ///
        ign::geometry::LineString CLInAreaGenerationOp::_convertPathToLineString(
            GraphType const& graph,
            std::string const& country,
            std::list<oriented_edge_descriptor> const& path
        ) const {
            ign::geometry::LineString pathGeom;

            ign::geometry::Point firstPoint = graph.getGeometry( graph.source( *path.begin() ) );

            std::string originId = _getOrigin(graph, country, path.begin()->descriptor);

            ign::feature::Feature featOrigin;
            _fsEdge->getFeatureById(originId, featOrigin);

            //DEBUG
            _logger->log(epg::log::DEBUG, "ouha1");
            _logger->log(epg::log::DEBUG, originId);
            _logger->log(epg::log::DEBUG, featOrigin.getId());
            _logger->log(epg::log::DEBUG, "ouha11");
            _logger->log(epg::log::DEBUG, featOrigin.getGeometry().toString());

            ign::geometry::LineString const * refGeom = &featOrigin.getGeometry().asLineString();

            _setZ(firstPoint, *refGeom),
            pathGeom.addPoint( firstPoint );
            
            edges_path_const_iterator pit;
            for( pit = path.begin() ; pit != path.end(); ++pit )
            {
                ign::geometry::LineString ls = graph.getGeometry( *pit );
                for( size_t i = 1 ; i < ls.numPoints() ; ++i ) {
                    if ( i == ls.numPoints()-1 ) {
                        std::string currentOriginId = _getOrigin(graph, country, pit->descriptor);
                        if( currentOriginId != originId ) {
                            originId = currentOriginId;
                            _fsEdge->getFeatureById(originId, featOrigin);

                            //DEBUG
                            _logger->log(epg::log::DEBUG, "ouha1");
                            _logger->log(epg::log::DEBUG, originId);
                            _logger->log(epg::log::DEBUG, featOrigin.getGeometry().toString());

                            refGeom = &featOrigin.getGeometry().asLineString();
                        }
                        _setZ(ls.pointN( i ), *refGeom);
                    }
                    pathGeom.addPoint( ign::geometry::Point( ls.pointN( i ).x(), ls.pointN( i ).y(), ls.pointN( i ).z() ) );
                }
            }

            return pathGeom;
        }

        ///
        ///
        ///
        std::string CLInAreaGenerationOp::_getOrigin(
            GraphType const& graph,
            std::string const& country_,
            edge_descriptor e
        ) const {
            std::vector<std::string> vOrigins = graph.origins(e);
            if( vOrigins.size() == 1 || country_ == "" ) return vOrigins.front();

            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            
            ign::feature::Feature featOrigin;
            _fsEdge->getFeatureById(vOrigins.front(), featOrigin);

            std::string country = featOrigin.getAttribute(countryCodeName).toString();
            if( country == country_ )
                return vOrigins.front();

            return vOrigins.back();  
        }

        ///
        ///
        ///
        void CLInAreaGenerationOp::_setZ(
            ign::geometry::Point & pt,
            ign::geometry::LineString const& refGeom
        ) const {
            ign::geometry::Point projPt;
            epg::tools::geometry::projectZ(refGeom, pt, projPt, 1e-5);
            pt.z() = projPt.z();
        }

		///
		///
		///
		std::map<double, std::string> CLInAreaGenerationOp::_getOriginEdges(
			GraphType const& graph,
            std::string const& country,
			std::list<oriented_edge_descriptor> const& path,
			std::map<std::string, std::set<edge_descriptor>> & mFeatMergedEdges,
            std::map<double, std::vector<detail::IncidentFeature>> & mIncidentFeatures
		) const {
            std::map<double, std::string> sAbsEdge;

            ign::geometry::LineString ls = _convertPathToLineString(graph, country, path);
            double length = ls.length();
            geometry::tools::LengthIndexedLineString lsIndex(ls);

            std::list<oriented_edge_descriptor>::const_iterator lit = path.begin();
            std::string previousFeat = graph.origins(lit->descriptor)[0];
			for ( ; lit != path.end() ; ++lit) {
				std::string currentFeat = graph.origins(lit->descriptor)[0];

                std::map<std::string, std::set<edge_descriptor>>::iterator mit = mFeatMergedEdges.find(currentFeat);
                if( mit == mFeatMergedEdges.end() )
                    mit = mFeatMergedEdges.insert(std::make_pair(currentFeat, std::set<edge_descriptor>())).first;
                mit->second.insert(lit->descriptor);
                    
                if ( currentFeat != previousFeat ) {

                    ign::geometry::Point ptStart = graph.getGeometry(graph.source(*lit));
                    int i = _getIndex(ls, ptStart);
                    double abs = lsIndex.getPointLocation(i)/length;

                    sAbsEdge.insert(std::make_pair(abs, previousFeat));

                    //on recupere les edges incidents
                    if( graph.degree(graph.source(*lit)) > 2 ) {
                        ign::geometry::Point vertexGeom = graph.getGeometry(graph.source(*lit));
                        std::vector< edge_descriptor > vEdges = graph.incidentEdges(graph.source(*lit));

                        for ( std::vector<edge_descriptor>::const_iterator vit = vEdges.begin() ; vit != vEdges.end() ; ++vit ) {

                            std::string incidentFeatId = graph.origins(*vit)[0];
                            if (incidentFeatId == currentFeat || incidentFeatId == previousFeat) continue;

                            std::map<double, std::vector<detail::IncidentFeature>>::iterator mit = mIncidentFeatures.find(abs);
                            if( mit == mIncidentFeatures.end() )
                                mit = mIncidentFeatures.insert(std::make_pair(abs, std::vector<detail::IncidentFeature>())).first;

                            ign::feature::Feature incidentFeat;
                            _fsEdge->getFeatureById(incidentFeatId, incidentFeat);
                            ign::geometry::LineString const& incidentFeatGeom = incidentFeat.getGeometry().asLineString();

                            bool touchStart = incidentFeatGeom.startPoint().distance(vertexGeom) < 1e-5;
                            bool touchEnd = incidentFeatGeom.endPoint().distance(vertexGeom) < 1e-5;

                            if (touchStart)
                                mit->second.push_back(detail::IncidentFeature(incidentFeatId, detail::START));
                            if (touchEnd)
                                mit->second.push_back(detail::IncidentFeature(incidentFeatId, detail::END));

                            if (!touchStart && !touchEnd) {
                                _logger->log(epg::log::ERROR, "Error in incident feature computation : "+vertexGeom.toString());
                            }
                        }
                    }
                }
                
                previousFeat = currentFeat;
			}
            sAbsEdge.insert(std::make_pair(1., previousFeat));

			return sAbsEdge;
		}

        ///
		///
		///
		int CLInAreaGenerationOp::_getIndex(
			ign::geometry::LineString const& ls,
            ign::geometry::Point const& pt
		) const {
            double minDist = std::numeric_limits<double>::max();
            int minIndex = -1;
            for (size_t i = 0; i < ls.numPoints(); ++i) {
                double dist = pt.distance(ls.pointN(i));
                if (dist < minDist) {
                    minDist = dist;
                    minIndex = i;
                }
            }
            return minIndex;
        }

		///
        ///
        ///
        void CLInAreaGenerationOp::_loadGraph(
            app::calcul::detail::EdgeCleaningGraphManager & graphManager, 
            ign::feature::FeatureFilter filter,
            bool planarize
            ) const 
        {
            graphManager.clear();

            //--
            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            //--
            app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG).getValue().toString();

            // chargement des edges
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
                    std::string wTag = fEdge.getAttribute(wTagName).toString();
                    graphManager.addEdgeSimple(ls, edgeId, OriginEdgeProperties(country, wTag, isCl));
                }
            }
            if (planarize) {
                graphManager.planarize();
                graphManager.createFaces();
            }
        }

		///
        ///
        ///
        bool CLInAreaGenerationOp::_pathsGeomAreEqual(
            ign::geometry::Polygon const& poly,
            ign::geometry::LineString & path1geom,
            ign::geometry::LineString & path2geom,
            double maxWidth,
            bool useHausdorff
        ) const {
            double hausdorffDist = ign::geometry::algorithm::HausdorffDistanceOp::distance(path1geom, path2geom);

            if (hausdorffDist < 0) {
                _logger->log(epg::log::ERROR, "Error in distance between paths calculation : hausdorff distance < 0");
                _logger->log(epg::log::ERROR, "path1 :" + path1geom.toString());
                _logger->log(epg::log::ERROR, "path2 :" + path2geom.toString());
            }

            if (useHausdorff) {
                return hausdorffDist >= 0 && hausdorffDist < maxWidth;
            }

            //--
            double lengthRation = 0.66666;
            double length1 = path1geom.length();
            double length2 = path2geom.length();
            bool bLength = length1 > lengthRation*length2 && length2 > lengthRation*length1;

            //--
            double meanWidth = 2 * ( poly.area() / poly.exteriorRing().length() );
            bool bHausdorff = hausdorffDist >= 0 && (hausdorffDist < 3*maxWidth);

            return meanWidth < maxWidth && bHausdorff && bLength;
        }

		///
        ///
        ///
        ign::geometry::LineString CLInAreaGenerationOp::_computeMeanPath(
            ign::geometry::LineString const& path1geom,
            ign::geometry::LineString const& path2geom
		) const {

			ign::geometry::LineString meanGeom;

			std::set<double> sAbsCurv;
			geometry::tools::LengthIndexedLineString lsIndex1(path1geom);
			geometry::tools::LengthIndexedLineString lsIndex2(path2geom);

			double l1 = path1geom.length();
			for (size_t i = 0; i < path1geom.numPoints(); ++i) {
				double abscurv = lsIndex1.getPointLocation(i)/l1;
				sAbsCurv.insert(abscurv);
			}
			double l2 = path2geom.length();
			for (size_t i = 0; i < path2geom.numPoints(); ++i) {
				double abscurv = lsIndex2.getPointLocation(i)/l2;
				sAbsCurv.insert(abscurv);
			}

			for (std::set<double>::iterator sit = sAbsCurv.begin(); sit != sAbsCurv.end(); ++sit) {
				ign::geometry::Point pt1 = lsIndex1.locateAlong(*sit*l1);
				ign::geometry::Point pt2 = lsIndex2.locateAlong(*sit*l2);
				ign::geometry::Point meanPt((pt1.x()+pt2.x())/2, (pt1.y()+pt2.y())/2, (pt1.z()+pt2.z())/2);
				meanGeom.addPoint(meanPt);
			}

			 return meanGeom;
		}

		///
        ///
        ///
        std::set<std::string> CLInAreaGenerationOp::_mergeFacePaths(
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
        std::list<app::calcul::detail::EdgeCleaningGraphManager::GraphType::oriented_edge_descriptor> CLInAreaGenerationOp::_getExteriorRingEdges(
            GraphType const& graph, 
            face_descriptor fd
        ) const {
            std::list<oriented_edge_descriptor> maxRing;
            double maxArea = 0;

            std::list< oriented_edge_descriptor > lEdges;

            oriented_edge_descriptor startEdge = graph.incidentEdge(fd);
            oriented_edge_descriptor nextEd = startEdge;

            do 
            {
                //arc degenere
                if ( graph[ nextEd.descriptor ].leftFace == graph[ nextEd.descriptor ].rightFace )
                {
                    nextEd = ign::geometry::graph::detail::nextEdge( nextEd, graph );
                    continue;
                }

                //on ajoute l arc
                lEdges.push_back( nextEd );

                //on cherche dans la pile si la cible de cet arc a deja ete ajoute
                //si oui, on depile en extrayant la boucle trouvee
                vertex_descriptor vTarget = graph.target( nextEd );
                if( graph.degree( vTarget ) < 3 && vTarget != graph.source( startEdge ) ) //optimisation
                {
                    nextEd = ign::geometry::graph::detail::nextEdge( nextEd, graph );
                    continue;
                }
                typename std::list< oriented_edge_descriptor >::iterator eit, eit_begin;
                for( eit = lEdges.begin() ; eit != lEdges.end() ; ++eit )
                {
                    if( graph.source( *eit ) == vTarget ) 
                    {
                        eit_begin = eit;
                        break;
                    }
                }
                //on a pas boucle, on passe a l arc suivant
                if( eit == lEdges.end() ) 
                {
                    nextEd = ign::geometry::graph::detail::nextEdge( nextEd, graph );
                    continue;
                }

                //on a boucle, on recupere la geometrie du cycle
                ign::geometry::LineString ls;
                ls.addPoint( graph[ graph.source( *eit ) ].point );
                for( ; eit != lEdges.end() ; ++eit )
                {
                    std::vector< ign::geometry::Point > const& vPoints = graph[ eit->descriptor ].intermediatePoints;
                    if( eit->direction == ign::graph::DIRECT )
                    {	
                        std::vector< ign::geometry::Point >::const_iterator it = vPoints.begin();
                        for( ; it != vPoints.end() ; ++it )
                            ls.addPoint( *it );
                    }
                    else{
                        std::vector< ign::geometry::Point >::const_reverse_iterator it = vPoints.rbegin();
                        for( ; it != vPoints.rend() ; ++it )
                            ls.addPoint( *it );
                    }

                    ls.addPoint( graph[ graph.target( *eit ) ].point );
                }

                double area = ign::geometry::algorithm::AreaOp::ComputeAlgebraicArea2d( ls );
                if ( area < maxArea ) //contour ext
                {
                    maxArea = area;
                    maxRing = lEdges;
                }
                
                //on depile
                // lEdges.resize( eit_begin - lEdges.begin() );
                // lEdges.resize( 5 );
                lEdges = std::list<oriented_edge_descriptor>(lEdges.begin(), eit_begin);
        
                nextEd = ign::geometry::graph::detail::nextEdge( nextEd, graph );
            } while ( nextEd != startEdge );
            return maxRing;
        }

		///
        ///
        ///
        bool CLInAreaGenerationOp::_getFacePaths(
            detail::EdgeCleaningGraphManager const& graphManager, 
            face_descriptor fd, 
            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> & vpCountryEdges
        ) const {
            GraphType const& graph = graphManager.getGraph();

            std::list<oriented_edge_descriptor> lEdges = _getExteriorRingEdges(graph, fd);

            //DEBUG
            if (lEdges.size() == 0) {
                return false;
            }

            std::string currentCountry = graphManager.getCountry(lEdges.begin()->descriptor);
            vpCountryEdges.push_back(std::make_pair(currentCountry, std::list<oriented_edge_descriptor>()));

            std::list<oriented_edge_descriptor>::const_iterator lit = lEdges.begin();
            std::list<oriented_edge_descriptor>::const_iterator lit_previous = lit;
            for ( ++lit ;  ; ++lit ) {
                vpCountryEdges.back().second.push_back(*lit_previous);
                if ( lit != lEdges.end() ) {
                    std::string nextCountry = graphManager.getCountry(lit->descriptor);
                    if ( graph.degree(graph.target(*lit_previous)) > 2 || currentCountry != nextCountry ) {
                        vpCountryEdges.push_back(std::make_pair(nextCountry, std::list<oriented_edge_descriptor>()));
                    }
                } else {
                    break;
                }
                lit_previous = lit;
            }

            if (vpCountryEdges.size() > 1 && vpCountryEdges.front().first == vpCountryEdges.back().first) {
                vertex_descriptor v = graph.target(vpCountryEdges.back().second.back());
                if ( graph.degree(v) == 2 ) {
                    for (std::list<oriented_edge_descriptor>::const_reverse_iterator rlit = vpCountryEdges.back().second.rbegin() ; rlit != vpCountryEdges.back().second.rend() ; ++rlit) {
                        vpCountryEdges.front().second.push_front(*rlit);
                    }
                    vpCountryEdges.pop_back();
                }
            }

            if ( vpCountryEdges.size() == 0 ) {
                _logger->log(epg::log::WARN, "No path found path [edge id] "+graph.origins(lEdges.begin()->descriptor)[0]);
                return false;
            }

            return true;
        }

	}
}