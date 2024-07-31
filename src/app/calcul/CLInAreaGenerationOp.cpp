#include <app/calcul/CLInAreaGenerationOp.h>

//APP
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LengthIndexedLineString.h>

//BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

//EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/sql/tools/numFeatures.h>
#include <epg/tools/geometry/project.h>

//SOCLE
#include <ign/geometry/algorithm/HausdorffDistanceOp.h>


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
            std::string listAttr2concatName = themeParameters->getValue(LIST_ATTR_TO_CONCAT).toString();
            std::string listAttrWName = themeParameters->getValue(LIST_ATTR_W).toString();
            std::string listAttrJsonName = themeParameters->getValue(LIST_ATTR_JSON).toString();
            _attrMerger.setLists(listAttr2concatName, listAttrWName, listAttrJsonName, "/");

            //--
            _logger->log(epg::log::INFO, "[END] initialization: " + epg::tools::TimeTools::getTime());
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
            // double const slimSurfaceWidth = themeParameters->getValue( ECL_SLIM_SURFACE_WIDTH ).toDouble();
            double const slimSurfaceWidth = 50;

            // std::vector<std::string> vCountry;
		    // epg::tools::StringTools::Split(_borderCode, "#", vCountry);

            detail::EdgeCleaningGraphManager graphManager;
            _loadGraph(graphManager);
            GraphType const& graph = graphManager.getGraph();

            boost::progress_display display(graph.numFaces(), std::cout, "[ generating CL in area  % complete ]\n");

            //--
            std::map<std::string, std::set<edge_descriptor>> mFeatMergedEdges;
            std::multimap<std::string, detail::IncidentFeature> mmIncidentFeatures;
            typedef std::multimap<std::string, detail::IncidentFeature>::const_iterator  m_iterator;

            face_iterator fit, fend;
            for( graph.faces( fit, fend ) ; fit != fend ; ++fit )
			{
                ++display;

				ign::geometry::Polygon faceGeom = graph.getGeometry( *fit );

                //DEBUG
                // if( faceGeom.intersects(ign::geometry::Point(4032840.73,3014651.08))) {
                //     bool test = true;
                // } else if( faceGeom.intersects(ign::geometry::Point(4032843.16,3014649.58))) {
                //     bool test = true;
                // } else {
                //     continue;
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

				ign::geometry::LineString lsFront = _convertPathToLineString(graph, vpCountryEdges.front().second);
                ign::geometry::LineString lsBack = _convertPathToLineString(graph, vpCountryEdges.back().second);

                if( _pathsGeomAreEqual(faceGeom, lsFront, lsBack, slimSurfaceWidth) ) {
					ign::geometry::LineString meanGeom = _computeMeanPath(lsFront, lsBack);
					geometry::tools::LengthIndexedLineString lsIndexMean(meanGeom);
                    double meanLength = meanGeom.length();

					//todo recupérer tous les edges des chemins
                    std::map<double, std::vector<detail::IncidentFeature>> mAbsIncidentFeatures;
					std::map<double, std::string> sAbsEdgeFront = _getOriginEdges(graph, vpCountryEdges.front().second, mFeatMergedEdges, mAbsIncidentFeatures); 
					std::map<double, std::string> sAbsEdgeBack = _getOriginEdges(graph, vpCountryEdges.back().second, mFeatMergedEdges, mAbsIncidentFeatures);

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
            std::set<std::string> sIncident2delete;
            std::vector<detail::IncidentFeature> vIncident2Create;
            for( std::map<std::string, std::set<edge_descriptor>>::const_iterator mit = mFeatMergedEdges.begin() ; mit != mFeatMergedEdges.end() ; ++mit ) {
                //DEBUG
                _logger->log(epg::log::DEBUG, "youp1");
                _logger->log(epg::log::DEBUG, mit->first);

                ign::feature::Feature featOrigin;
                _fsEdge->getFeatureById(mit->first, featOrigin);
                ign::geometry::LineString geomOrigin = featOrigin.getGeometry().asLineString();
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
                    ign::geometry::LineString newGeom = _convertPathToLineString(graph, vPaths[i]);

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
            for( std::multimap<std::string, detail::IncidentFeature>::const_iterator mit = mmIncidentFeatures.begin() ; mit != mmIncidentFeatures.end() ; ++mit ) {
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
		}

        ///
        ///
        ///
        ign::geometry::LineString CLInAreaGenerationOp::_convertPathToLineString(
            GraphType const& graph,
            std::list<oriented_edge_descriptor> const& path
        ) const {
            ign::geometry::LineString pathGeom;

            ign::geometry::Point firstPoint = graph.getGeometry( graph.source( *path.begin() ) );

            std::string originId = graph.origins(path.begin()->descriptor)[0];
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
                        std::string currentOriginId = graph.origins(pit->descriptor)[0];
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
			std::list<oriented_edge_descriptor> const& path,
			std::map<std::string, std::set<edge_descriptor>> & mFeatMergedEdges,
            std::map<double, std::vector<detail::IncidentFeature>> & mIncidentFeatures
		) const {
            std::map<double, std::string> sAbsEdge;

            ign::geometry::LineString ls = _convertPathToLineString(graph, path);
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
            ign::feature::FeatureFilter filter
            ) const 
        {
            graphManager.clear();

            epg::Context *context = epg::ContextS::getInstance();
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

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

                graphManager.addEdge(ls, edgeId, OriginEdgeProperties(country, isCl));
            }
			graphManager.planarize();
			graphManager.createFaces();
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
            double meanWidth = 2 * ( poly.area() / poly.exteriorRing().length() );
            bool bHausdorff = hausdorffDist >= 0 && (hausdorffDist < 3*maxWidth);

            return meanWidth < maxWidth && bHausdorff;
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
        bool CLInAreaGenerationOp::_getFacePaths(
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

                if (graph.degree(graph.target(currentEdge)) > 2 || currentCountry != nextCountry ) {
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
                if ( graph.degree(v) == 2 ) {
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

	}
}