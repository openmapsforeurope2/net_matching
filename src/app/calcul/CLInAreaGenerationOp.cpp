// APP
#include <app/calcul/CLInAreaGenerationOp.h>
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LengthIndexedLineString.h>
#include <app/calcul/detail/graph/concept/EdgeCleaningGraphSpecializations.h>

// BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

// EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/sql/tools/numFeatures.h>
#include <epg/tools/geometry/project.h>
#include <epg/graph/tools/reverse.h>

// OME2
#include <ome2/feature/sql/NotDestroyedTools.h>

// SOCLE
#include <ign/geometry/algorithm/HausdorffDistanceOp.h>
#include <ign/geometry/algorithm/AreaOp.h>
#include <ign/geometry/graph/detail/NextEdge.h>

// STL
#include <queue>


namespace app
{
    namespace calcul
    {

		///
        ///
        ///
        void CLInAreaGenerationOp::Compute(
            bool verbose
        ) {
            CLInAreaGenerationOp cLInAreaGenerationOp(verbose);
            cLInAreaGenerationOp._computeByIteration();
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
            std::string const wTagName = themeParameters->getParameter(W_TAG_NAME).getValue().toString();
            context->getDataBaseManager().resetColumn(_fsEdge->getTableName(), wTagName);

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
        std::pair<bool, std::pair<std::string, std::string>> CLInAreaGenerationOp::_findClEdgeAssociation(
            detail::EdgeCleaningGraphManager const& graphManager,
            face_descriptor f,
            edge_descriptor eCl
        ) const {
            std::string assEdgeId = "";

            //--
            GraphType const& graph = graphManager.getGraph();

            //--
            app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
            std::string const natIdName = themeParameters->getParameter(NATIONAL_IDENTIFIER_NAME).getValue().toString();

            //--
            std::string clNatId = graphManager.getNatId(eCl);

            if( graph.degree(graph.source(eCl)) == 3 ) {
                std::vector<edge_descriptor> vEdges = graph.incidentEdges(graph.source(eCl));
                for (size_t i = 0 ; i < vEdges.size() ; ++i ) {
                    if( vEdges[i] == eCl ) continue;
                    if( _belongToFace(graph, f, vEdges[i]) ) {
                        std::string edgeNatId = graphManager.getNatId(vEdges[i]);
                        if (clNatId.find(edgeNatId) != std::string::npos) {
                            assEdgeId = graph.origins(vEdges[i])[0];
                        }
                    }
                }
            }

            if( graph.degree(graph.target(eCl)) == 3 ) {
                std::vector<edge_descriptor> vEdges = graph.incidentEdges(graph.target(eCl));
                for (size_t i = 0 ; i < vEdges.size() ; ++i ) {
                    if( vEdges[i] == eCl ) continue;
                    if( _belongToFace(graph, f, vEdges[i]) ) {
                        std::string edgeNatId = graphManager.getNatId(vEdges[i]);
                        if (clNatId.find(edgeNatId) != std::string::npos) {
                            if (assEdgeId != "") {
                                return std::make_pair(false, std::make_pair("",""));
                            }
                            assEdgeId = graph.origins(vEdges[i])[0];
                        }
                    }
                }
            }
            if(assEdgeId == "")
                return std::make_pair(false, std::make_pair("",""));

            std::string clId = graph.origins(eCl)[0];
            return std::make_pair(true, std::make_pair(clId,assEdgeId));
        }

        ///
        ///
        ///
        bool CLInAreaGenerationOp::_belongToFace(
            GraphType const& graph,
            face_descriptor f,
            edge_descriptor e
        ) const {
            std::list<oriented_edge_descriptor> lEdges = _getExteriorRingEdges(graph, f);
            for( std::list<oriented_edge_descriptor>::const_iterator lit = lEdges.begin() ; lit != lEdges.end() ; ++lit ) {
                if( lit->descriptor == e ) return true;
            }
            return false;
        }

        ///
        ///
        ///
        double CLInAreaGenerationOp::_isFaceToTreat(
            detail::EdgeCleaningGraphManager const& graphManager,
            face_descriptor f,
            std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> & vpCountryEdges,
            std::vector<ign::geometry::LineString> & vPathsGeom,
            std::map<std::string, std::string> & mClToConvert,
            double slimSurfaceWidth
        ) const {
            GraphType const& graph = graphManager.getGraph();

            ign::geometry::Polygon faceGeom = graph.getGeometry( f );

            if (!_getFacePaths(graphManager, f, vpCountryEdges))
                return -1;

            if (vpCountryEdges.size() < 2) 
                return -1;

            bool foundCl = false;
            for (size_t i = 0 ; i < vpCountryEdges.size() ; ++i ) {
                if (vpCountryEdges[i].first.find("#") != std::string::npos) {
                    if( vpCountryEdges[i].second.size() == 1 ) {
                        std::pair<bool, std::pair<std::string, std::string>> foundClEdgeAss = _findClEdgeAssociation(graphManager, f, vpCountryEdges[i].second.begin()->descriptor);
                        if( foundClEdgeAss.first ) {
                            mClToConvert.insert(foundClEdgeAss.second);
                            vpCountryEdges[i] = std::make_pair(graphManager.getCountry(foundClEdgeAss.second.second), vpCountryEdges[i].second );
                        } else {
                            foundCl = true;
                            break;
                        }
                    } else {
                        foundCl = true;
                        break;
                    }
                }
            }
            if( foundCl )
                return -1;

            std::set<std::string> sFaceCountries;
            for (std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>>::const_iterator vpit = vpCountryEdges.begin() ; vpit != vpCountryEdges.end() ; ++vpit)
                sFaceCountries.insert(vpit->first);

            if (sFaceCountries.size() != 2) 
                return -1;

            _mergeFacePaths(vpCountryEdges);
            
            if (vpCountryEdges.size() != 2)
                return -1;

            // mettre les 2 chemins dans le meme sens
            if (graph.source(*vpCountryEdges.front().second.begin()) != graph.source(*vpCountryEdges.back().second.begin()))
                vpCountryEdges.back().second = _getReversePath(vpCountryEdges.back().second);

            vPathsGeom.push_back(_convertPathToLineString(graph, vpCountryEdges.front().first, vpCountryEdges.front().second));
            vPathsGeom.push_back(_convertPathToLineString(graph, vpCountryEdges.back().first, vpCountryEdges.back().second));

            return _pathsGeomAreEqual(faceGeom, vPathsGeom.front(), vPathsGeom.back(), slimSurfaceWidth);
        }

        ///
        ///
        ///
        void CLInAreaGenerationOp::_replaceClByAssociatedEdge(
            std::map<size_t, std::string> & mMeanIndexEdge,
            std::map<std::string, std::string> const& mClToConvert
        ) const {
            std::map<size_t, std::string>::iterator mit = mMeanIndexEdge.begin();
            for ( ; mit != mMeanIndexEdge.end() ; ++mit ) {
                std::map<std::string, std::string>::const_iterator mit2 = mClToConvert.find(mit->second);
                if ( mit2 == mClToConvert.end() ) continue;
                mit->second = mit2->second;
            }
        }

        ///
        ///
        ///
        bool CLInAreaGenerationOp::_createCLOnFaces(
            detail::EdgeCleaningGraphManager & graphManager,
            std::map<std::string, std::set<edge_descriptor>> & mFeatMergedEdges,
            std::multimap<std::string, detail::IncidentFeature> & mmIncidentFeatures,
            std::set<edge_descriptor> & sTreatedEdges
        ) const {
            bool hasCreatedCl = false;

            //--
			epg::Context *context = epg::ContextS::getInstance();

            //--
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
            
            //--
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG_NAME).getValue().toString();
            double const slimSurfaceWidth = themeParameters->getValue( CLA_SURFACE_WIDTH ).toDouble();

            //--
            GraphType const& graph = graphManager.getGraph();
            
            // get faces to treate ordered by width
            std::multimap< double, face_descriptor> mWidthFace;
            std::map<face_descriptor, std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>>> mvpCountryEdges;
            std::map<face_descriptor, std::vector<ign::geometry::LineString>> mvPathsGeom;
            std::map<face_descriptor, std::map<std::string, std::string>> mmClToConvert;
            face_iterator fit, fend;
            for( graph.faces( fit, fend ) ; fit != fend ; ++fit )
			{
                //DEBUG
                // ign::geometry::Polygon faceGeom = graph.getGeometry( *fit );
                // if(faceGeom.intersects(ign::geometry::Point(3370768.386,2318558.186))) {
                //     bool test = true;
                // }
                // if(faceGeom.intersects(ign::geometry::Point(3370775.604,2318549.250))) {
                //     bool test = true;
                // }

                std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> vpCountryEdges;
                std::vector<ign::geometry::LineString> vPathsGeom;
                std::map<std::string, std::string> mClToConvert;
                double width = _isFaceToTreat(graphManager, *fit, vpCountryEdges, vPathsGeom, mClToConvert, slimSurfaceWidth);

                if (width < 0) continue;

                mWidthFace.insert( std::make_pair(width, *fit));
                mvpCountryEdges.insert( std::make_pair(*fit, vpCountryEdges));
                mvPathsGeom.insert( std::make_pair(*fit, vPathsGeom));
                mmClToConvert.insert( std::make_pair(*fit, mClToConvert));
            }

            //--
            boost::progress_display display(mWidthFace.size(), std::cout, "[ generating CL in area  % complete ]\n");
            for( std::multimap< double, face_descriptor>::const_iterator mmit = mWidthFace.begin() ; mmit != mWidthFace.end() ; ++mmit )
			{
                ++display;

                //DEBUG
                // ign::geometry::Polygon faceGeom = graph.getGeometry( mmit->second );
                // if(faceGeom.intersects(ign::geometry::Point(3370768.386,2318558.186))) {
                //     bool test = true;
                // }
                // if(faceGeom.intersects(ign::geometry::Point(3370775.604,2318549.250))) {
                //     bool test = true;
                // }

                //--
                std::vector<std::pair<std::string, std::list<oriented_edge_descriptor>>> const& vpCountryEdges = mvpCountryEdges[mmit->second] ;
                std::vector<ign::geometry::LineString> const& vPathsGeom = mvPathsGeom[mmit->second] ;

                //--
                if ( _hasTreatedEdge(graph,vpCountryEdges.front().second, sTreatedEdges) || _hasTreatedEdge(graph,vpCountryEdges.back().second, sTreatedEdges) )
                    continue;

                //--
                bool isFictitiousFront = _isFictitious(graph, vpCountryEdges.front().first, vpCountryEdges.front().second);
                bool isFictitiousBack = _isFictitious(graph, vpCountryEdges.back().first, vpCountryEdges.back().second);

                std::vector<int> pathFrontgeomMeanGeomLink;
                std::vector<int> pathBackgeomMeanGeomLink;
                ign::geometry::LineString meanGeom;
                if(isFictitiousFront && !isFictitiousBack) {
                    meanGeom = _addInterMeaningFullPoint(
                        graph, 
                        vpCountryEdges.back().second,
                        vPathsGeom.front(),
                        vPathsGeom.back(),
                        pathFrontgeomMeanGeomLink,
                        pathBackgeomMeanGeomLink
                    );
                } else if (!isFictitiousFront && isFictitiousBack) {
                    meanGeom = _addInterMeaningFullPoint(
                        graph, 
                        vpCountryEdges.front().second,
                        vPathsGeom.back(),
                        vPathsGeom.front(),
                        pathBackgeomMeanGeomLink,
                        pathFrontgeomMeanGeomLink
                    );
                } else {
                    meanGeom = _computeMeanPath(
                        vPathsGeom.front(),
                        vPathsGeom.back(),
                        pathFrontgeomMeanGeomLink,
                        pathBackgeomMeanGeomLink
                    );
                }
                // geometry::tools::LengthIndexedLineString lsIndexMean(meanGeom);
                // double meanLength = meanGeom.length();

                //todo recupérer tous les edges des chemins
                std::map<size_t, std::vector<detail::IncidentFeature>> mMeanIndexIncidentFeatures;
                std::map<size_t, std::string> mMeanIndexEdgeFront = _getOriginEdges(graph, vpCountryEdges.front().first, vpCountryEdges.front().second, pathFrontgeomMeanGeomLink, mFeatMergedEdges, mMeanIndexIncidentFeatures); 
                std::map<size_t, std::string> mMeanIndexEdgeBack = _getOriginEdges(graph, vpCountryEdges.back().first,vpCountryEdges.back().second, pathBackgeomMeanGeomLink, mFeatMergedEdges, mMeanIndexIncidentFeatures);

                _replaceClByAssociatedEdge(mMeanIndexEdgeFront, mmClToConvert[mmit->second]);
                _replaceClByAssociatedEdge(mMeanIndexEdgeBack, mmClToConvert[mmit->second]);

                size_t meanIndexStart = 0;
                do {
                    size_t meanIndexFront = mMeanIndexEdgeFront.begin()->first;
                    size_t meanIndexBack = mMeanIndexEdgeBack.begin()->first;
                    ign::feature::Feature featFront;
                    _fsEdge->getFeatureById(mMeanIndexEdgeFront.begin()->second, featFront);
                    ign::feature::Feature featBack;
                    _fsEdge->getFeatureById(mMeanIndexEdgeBack.begin()->second, featBack);

                    ign::geometry::LineString newFeatGeom;
                    size_t meanIndexTarget;
                    if( meanIndexFront < meanIndexBack ) {
                        mMeanIndexEdgeFront.erase(mMeanIndexEdgeFront.begin());
                        meanIndexTarget = meanIndexFront;
                    } else if ( meanIndexFront > meanIndexBack ) {
                        mMeanIndexEdgeBack.erase(mMeanIndexEdgeBack.begin());
                        meanIndexTarget = meanIndexBack;
                    } else {
                        mMeanIndexEdgeFront.erase(mMeanIndexEdgeFront.begin());
                        mMeanIndexEdgeBack.erase(mMeanIndexEdgeBack.begin());
                        meanIndexTarget = meanIndexFront;
                    }

                    //pour eviter pbls de précision entre meanLength et dernier index de lsIndexMean
                    // double indexTarget = sTarget >= 1 ? lsIndexMean.getPointAbscisses().back() : sTarget*meanLength;

                    // ign::geometry::LineString lsNew = lsIndexMean.getSubLineString(sStart*meanLength, indexTarget);
                    ign::geometry::LineString lsNew = _getSubLineString(meanGeom, meanIndexStart, meanIndexTarget);

                    //--
                    if (meanIndexStart != 0) {
                        std::map<size_t, std::vector<detail::IncidentFeature>>::iterator mit = mMeanIndexIncidentFeatures.find(meanIndexStart);
                        if ( mit != mMeanIndexIncidentFeatures.end() ) {
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

                    if( lsNew.isEmpty() || lsNew.isNull() || !lsNew.isValid()) {
                        _logger->log(epg::log::ERROR, "Resulting inconsistent geometry [face geom]: "+graph.getGeometry( mmit->second ).toString());
                        continue;
                    }

                    hasCreatedCl = true;
                    _recordTreatedEdge(graph, vpCountryEdges.front().second, sTreatedEdges);
                    _recordTreatedEdge(graph, vpCountryEdges.back().second, sTreatedEdges);

                    fRef->setGeometry(lsNew);

                    _fsEdge->createFeature(*fRef);
                    meanIndexStart = meanIndexTarget;
                } while(mMeanIndexEdgeFront.size() > 0 && mMeanIndexEdgeBack.size() > 0);
			}

            return hasCreatedCl;
        }


        ign::geometry::LineString CLInAreaGenerationOp::_addInterMeaningFullPoint(
            GraphType const& graph,
            std::list<oriented_edge_descriptor> const& pathOther,
            ign::geometry::LineString const& lsRef,
            ign::geometry::LineString const& lsOther,
            std::vector<int> & geomRefgeomMeanLink,
            std::vector<int> & geomOthergeomMeanLink
        ) const {
            double snapDist = 1;

            geomRefgeomMeanLink.resize(lsRef.numPoints(), -1);
            geomOthergeomMeanLink.resize(lsOther.numPoints(), -1);

            double lengthOther = lsOther.length();
            geometry::tools::LengthIndexedLineString lsIndexOther(lsOther);
            
            std::queue<std::pair<double, size_t>> fifoQueueNewInterPt;
            std::list<oriented_edge_descriptor>::const_iterator lit = pathOther.begin();
			for ( ++lit ; lit != pathOther.end() ; ++lit) {
                ign::geometry::Point sourcePoint = graph.getGeometry(graph.source(*lit));
                int i = _getIndex(lsOther, sourcePoint);
                double abs = lsIndexOther.getPointLocation(i)/lengthOther;

                fifoQueueNewInterPt.push(std::make_pair(abs, i));
            }

            double lengthRef = lsRef.length();
            geometry::tools::LengthIndexedLineString lsIndexRef(lsRef);

            std::vector<double> vAbsRef = lsIndexRef.getPointAbscisses();

            ign::geometry::LineString lsMean;
            for ( size_t i = 1 ; i < vAbsRef.size() ; ++i ) {

                lsMean.addPoint(lsRef.pointN(i-1));
                geomRefgeomMeanLink[i-1] = lsMean.numPoints()-1;

                if (fifoQueueNewInterPt.empty()) 
                    continue;

                double nextAbs = fifoQueueNewInterPt.front().first*lengthRef;

                if( nextAbs < vAbsRef[i] + snapDist ) {
                    geomOthergeomMeanLink[fifoQueueNewInterPt.front().second] = lsMean.numPoints();
                    fifoQueueNewInterPt.pop();

                    if( nextAbs < vAbsRef[i] - snapDist ) {
                        ign::geometry::Point newPt = lsIndexRef.locateAlong(nextAbs);
                        lsMean.addPoint(newPt);
                    }
                }
            }
            lsMean.addPoint(lsRef.endPoint());
            geomRefgeomMeanLink.back() = lsMean.numPoints()-1;
            geomOthergeomMeanLink.back() = lsMean.numPoints()-1;

            return lsMean;
        }

        ///
        ///
        ///
        ign::geometry::LineString CLInAreaGenerationOp::_getSubLineString(
            ign::geometry::LineString const& ls,
            size_t indexStart,
            size_t indexEnd
        ) const {
            ign::geometry::LineString subLs;

            if( indexStart >= indexEnd || indexEnd > ls.numPoints()-1 ) return subLs;

            for ( size_t i = indexStart ; i <= indexEnd ; ++i )
                subLs.addPoint(ls.pointN(i));

            return subLs;
        }

        ///
        ///
        ///
        bool CLInAreaGenerationOp::_hasTreatedEdge(
			GraphType const& graph,
			std::list<oriented_edge_descriptor> const& path,
            std::set<edge_descriptor> const& sTreatedEdges
		) const {
            std::list<oriented_edge_descriptor>::const_iterator lit = path.begin();
            edge_descriptor previousEdge = lit->descriptor;
			for ( ; lit != path.end() ; ++lit) {
                edge_descriptor currentEdge = lit->descriptor;

				if ( sTreatedEdges.find(currentEdge) != sTreatedEdges.end() ) {
                    return true;
                }

                if ( currentEdge == previousEdge ) continue;

                if (graph.degree(graph.source(*lit)) > 2) {
                    std::vector< edge_descriptor > vIncidentEdges = graph.incidentEdges( graph.source(*lit) );
                    for( size_t i = 0 ; i < vIncidentEdges.size() ; ++i )  {
                        if(vIncidentEdges[i] != currentEdge && vIncidentEdges[i] != previousEdge ) {
                            if ( sTreatedEdges.find(vIncidentEdges[i]) != sTreatedEdges.end() ) {
                                return true;
                            }
                        }
                    }
                }

                previousEdge = currentEdge;
            }

            return false;
        }

        ///
        ///
        ///
        void CLInAreaGenerationOp::_recordTreatedEdge(
			GraphType const& graph,
			std::list<oriented_edge_descriptor> const& path,
            std::set<edge_descriptor> & sTreatedEdges
		) const {
            std::list<oriented_edge_descriptor>::const_iterator lit = path.begin();
            edge_descriptor previousEdge = lit->descriptor;
			for ( ; lit != path.end() ; ++lit) {
				edge_descriptor currentEdge = lit->descriptor;

                sTreatedEdges.insert(currentEdge);
                    
                if ( currentEdge == previousEdge ) continue;

                if (graph.degree(graph.source(*lit)) > 2) {
                    std::vector< edge_descriptor > vIncidentEdges = graph.incidentEdges( graph.source(*lit) );
                    for( size_t i = 0 ; i < vIncidentEdges.size() ; ++i )  {
                        if(vIncidentEdges[i] != currentEdge && vIncidentEdges[i] != previousEdge ) {
                            sTreatedEdges.insert(vIncidentEdges[i]);
                        }
                    }
                }

                previousEdge = currentEdge;
            }
        }

        ///
        ///
        ///
        bool CLInAreaGenerationOp::_createCLOnOverlappingEdges(
            GraphType const& graph,
            std::map<std::string, std::set<edge_descriptor>> & mFeatMergedEdges,
            std::multimap<std::string, detail::IncidentFeature> & mmIncidentFeatures,
            std::set<edge_descriptor> & sTreatedEdges
        ) const {
            bool hasCreatedCl = false;

            //--
			epg::Context *context = epg::ContextS::getInstance();

			//--
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            //--
            app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG_NAME).getValue().toString();

            boost::progress_display display(graph.numEdges(), std::cout, "[ generating CL on overlapping edges % complete ]\n");

            edge_iterator eit, eend;
            for (graph.edges(eit, eend); eit != eend; ++eit)
            {
                ++display;

                if (sTreatedEdges.find(*eit) != sTreatedEdges.end())
                    continue;

                std::vector<std::string> vOrigins = graph.origins(*eit);
                if( vOrigins.size() != 2 ) continue;

                //--
                ign::geometry::LineString edgeGeom = graph.getGeometry( *eit );

                //--
                ign::feature::Feature featFront;
                _fsEdge->getFeatureById(vOrigins.front(), featFront);
                ign::feature::Feature featBack;
                _fsEdge->getFeatureById(vOrigins.back(), featBack);


				//on verifie que l'on ne fusionne pas un obj d'un pays deja fusionne dans une CL
				if (featFront.getAttribute(countryCodeName).toString().find(
					featBack.getAttribute(countryCodeName).toString() ) != std::string::npos
					 || featBack.getAttribute(countryCodeName).toString().find(
						 featFront.getAttribute(countryCodeName).toString()) != std::string::npos
					) {
					_logger->log(epg::log::WARN, " fusion impossible de deux troncons d'un même pays pour la geom : " +featFront.getGeometry().asText() + " et la geom : " + featBack.getGeometry().asText());
					continue;
				}

                hasCreatedCl = true;
                sTreatedEdges.insert(*eit);


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

                _setZ(featFront.getGeometry().asLineString(), featBack.getGeometry().asLineString(), edgeGeom);

                fRef->setGeometry(edgeGeom);
                
				_fsEdge->createFeature(*fRef);
            }

            return hasCreatedCl;
        }

        ///
        ///
        ///
        void CLInAreaGenerationOp::_setZ(
            ign::geometry::LineString const& featGeom1,
            ign::geometry::LineString const& featGeom2,
            ign::geometry::LineString & edgeGeom
        ) const {
            if(ign::numeric::Numeric<double>::IsNaN(edgeGeom.startPoint().z() )){
                ign::geometry::Point projPt1;
                epg::tools::geometry::projectZ(featGeom1, edgeGeom.startPoint(), projPt1, 1e-5);

                ign::geometry::Point projPt2;
                epg::tools::geometry::projectZ(featGeom2, edgeGeom.startPoint(), projPt2, 1e-5);

                edgeGeom.startPoint().z() = (projPt1.z() + projPt2.z()) / 2;
            }
            if(ign::numeric::Numeric<double>::IsNaN(edgeGeom.endPoint().z() )){
                ign::geometry::Point projPt1;
                epg::tools::geometry::projectZ(featGeom1, edgeGeom.endPoint(), projPt1, 1e-5);

                ign::geometry::Point projPt2;
                epg::tools::geometry::projectZ(featGeom2, edgeGeom.endPoint(), projPt2, 1e-5);

                edgeGeom.endPoint().z() = (projPt1.z() + projPt2.z()) / 2;
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
        void CLInAreaGenerationOp::_computeByIteration() const
        {
            while (_compute()) {}
            _clean();
        }

        ///
        ///
        ///
        void CLInAreaGenerationOp::_clean() const
        {
            ign::feature::FeatureFilter filter;
            detail::EdgeCleaningGraphManager graphManager;

            _loadGraph(graphManager, filter);

            GraphType const& graph = graphManager.getGraph();

            vertex_iterator vit, vend;
            for( graph.vertices( vit, vend ) ; vit != vend ; ++vit )
            {
                if (graph.degree(*vit) < 3)
                    continue;

                std::vector<edge_descriptor> vCls;
                std::vector<edge_descriptor> vEdges = graph.incidentEdges(*vit);
                for( size_t i = 0 ; i < vEdges.size() ; ++i ) {
                    if( graphManager.isCl(vEdges[i]) )
                        vCls.push_back(vEdges[i]);
                }

                if( vCls.size() != 2 )
                    continue;

                std::string idNatFront = graphManager.getNatId(vCls.front());
                std::string idNatBack = graphManager.getNatId(vCls.back());

                if( idNatFront != idNatBack )
                    continue;
                
                ign::geometry::LineString lsFront = graph.getGeometry(vCls.front());
                ign::geometry::LineString lsBack = graph.getGeometry(vCls.back());

                edge_descriptor* eToMergePtr = lsFront.length() < lsBack.length() ? &vCls.front() : &vCls.back();

                vertex_descriptor vConn = graph.source(*eToMergePtr) == *vit ? graph.target(*eToMergePtr) : graph.source(*eToMergePtr);

                if( graph.degree(vConn) != 2 )
                    continue;
                
                std::vector<edge_descriptor> vEdges2 = graph.incidentEdges(vConn);
                edge_descriptor eRef = vEdges2.front() == *eToMergePtr ? vEdges2.back() : vEdges2.front();

                if ( !graphManager.isCl(eRef) )
                    continue;

                //--
                std::string idFeatRef = graph.origins(eRef)[0];
                std::string idFeatToMerge = graph.origins(*eToMergePtr)[0];

                //--
                ign::feature::Feature fRef, fToMerge;
                _fsEdge->getFeatureById(idFeatRef, fRef);
                _fsEdge->getFeatureById(idFeatToMerge, fToMerge);

                //--
                ign::geometry::LineString const& featRefGeom = fRef.getGeometry().asLineString();
                ign::geometry::LineString const& fToMergeGeom = fToMerge.getGeometry().asLineString();

                ign::geometry::LineString mergedGeom = _merge(featRefGeom, fToMergeGeom);

                //--
                fRef.setGeometry(mergedGeom);
                _fsEdge->modifyFeature(fRef);

                //--
                _fsEdge->deleteFeature(idFeatToMerge);
            }
        }

        ///
        ///
        ///
        ign::geometry::LineString CLInAreaGenerationOp::_merge(
            ign::geometry::LineString const& lsRef,
            ign::geometry::LineString const& lsToMerge
        ) const {
            ign::geometry::LineString lsMerged;
            double distEndStart = lsRef.endPoint().distance(lsToMerge.startPoint());
            double distEndEnd = lsRef.endPoint().distance(lsToMerge.endPoint());
            if( distEndStart < 1e-5 || distEndEnd < 1e-5 ) {
                lsMerged = lsRef;
                lsMerged.removePointN(lsMerged.numPoints()-1);

                if( distEndStart < 1e-5 )
                    for ( size_t i = 1 ; i < lsToMerge.numPoints() ; ++i ) 
                        lsMerged.addPoint(lsToMerge.pointN(i));
                else
                    for ( int i = lsToMerge.numPoints()-2 ; i >= 0 ; --i ) 
                        lsMerged.addPoint(lsToMerge.pointN(i));
            } else {
                double distStartStart = lsRef.startPoint().distance(lsRef.startPoint());

                lsMerged = distStartStart < 1e-5 ? lsToMerge.getReversed() : lsToMerge;
                lsMerged.removePointN(lsToMerge.numPoints()-1);

                for ( size_t i = 1 ; i < lsRef.numPoints() ; ++i ) 
                        lsMerged.addPoint(lsRef.pointN(i));

            }
            return lsMerged;
        }

		///
        ///
        ///
        bool CLInAreaGenerationOp::_compute() const
        {
			//--
			epg::Context *context = epg::ContextS::getInstance();

			//--
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

			//--
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG_NAME).getValue().toString();

            //--
            // ign::feature::FeatureFilter filter(countryCodeName+" NOT LIKE '%#%'");
            ign::feature::FeatureFilter filter;
            detail::EdgeCleaningGraphManager graphManager;

            //--
            // ign::feature::FeatureFilter filterCountCl(countryCodeName+" LIKE '%#%'");
            // size_t numCls = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterCountCl);
            // if( numCls > 0 ) {
            //     _loadGraph(graphManager, filter, false);
            //     _mergeCl(graphManager);
            // }
            
            //--
            _loadGraph(graphManager, filter);
            GraphType const& graph = graphManager.getGraph();

            //--
            std::multimap<std::string, detail::IncidentFeature> mmIncidentFeatures;
            std::map<std::string, std::set<edge_descriptor>> mFeatMergedEdges;

            std::set<edge_descriptor> sTreatedEdges;
            //--
            bool hasCreatedCl1 = _createCLOnOverlappingEdges(graph, mFeatMergedEdges, mmIncidentFeatures, sTreatedEdges);
            //--
            bool hasCreatedCl2 = _createCLOnFaces(graphManager, mFeatMergedEdges, mmIncidentFeatures, sTreatedEdges);

            //--
            if( !hasCreatedCl1 && !hasCreatedCl2 )
                return false;

            // on parcours tous les features qui ont des edges induits mergés : on les split si nécessaire
            boost::progress_display display2(mFeatMergedEdges.size(), std::cout, "[ splitting features  % complete ]\n");

            std::set<std::string> sIncident2delete;
            std::vector<detail::IncidentFeature> vIncident2Create;
            for( std::map<std::string, std::set<edge_descriptor>>::const_iterator mit = mFeatMergedEdges.begin() ; mit != mFeatMergedEdges.end() ; ++mit ) {
                ++display2;

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

                    if( newGeom.isEmpty() || newGeom.isNull() || !newGeom.isValid()) {
                        _logger->log(epg::log::ERROR, "Resulting inconsistent geometry [orign feature id]: "+mit->first);
                        continue;
                    }

                    featOrigin.setGeometry(newGeom);
                    featOrigin.setAttribute(wTagName, ign::data::String(mit->first));
                    _fsEdge->createFeature(featOrigin);

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
                }
            }

            //maj incident edges
            for( std::set<std::string>::const_iterator sit = sIncident2delete.begin() ; sit != sIncident2delete.end() ; ++sit ) {
                mmIncidentFeatures.erase(*sit);
            }
            for( std::vector<detail::IncidentFeature>::const_iterator vit = vIncident2Create.begin() ; vit != vIncident2Create.end() ; ++vit ) {
                mmIncidentFeatures.insert(std::make_pair(vit->originId, *vit));
            }

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

            //--
            for( std::map<std::string, std::set<edge_descriptor>>::const_iterator mit = mFeatMergedEdges.begin() ; mit != mFeatMergedEdges.end() ; ++mit )
                _fsEdge->deleteFeature(mit->first);

            //--
            _cleanOverLappingCl();

            //--
            _mergeByWTag();

            return true;
		}

        ///
        ///
        ///
        void CLInAreaGenerationOp::_cleanOverLappingCl() const {
            //--
			epg::Context *context = epg::ContextS::getInstance();

			//--
            epg::params::EpgParameters const& epgParams = context->getEpgParameters();
            std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

            //--
            app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG_NAME).getValue().toString();

            //--
            detail::EdgeCleaningGraphManager graphManager;
            ign::feature::FeatureFilter filter(wTagName+" LIKE '%#%'");
            _loadGraph(graphManager, filter, true);
            GraphType const& graph = graphManager.getGraph();

            //patience
            boost::progress_display display( graph.numEdges() , std::cout, "[ cleaning overlapping CL % complete ]\n") ;

            std::map<std::string, std::set<edge_descriptor>> mFeatMergedEdges;
            std::set< edge_descriptor > sMergedEdges;

            edge_iterator eit, eend ;
            for( graph.edges( eit, eend ) ; eit != eend ; ++eit )
            {
                ++display;

                edge_descriptor edgePivot = *eit ;
                std::vector<std::string> vOriginsPivot = graph.origins(edgePivot);

                if( vOriginsPivot.size() == 1 ) continue;
                if( sMergedEdges.find(edgePivot) != sMergedEdges.end() ) continue ;

                oriented_edge_descriptor tPivot[] = { 
                    oriented_edge_descriptor( edgePivot, ign::graph::DIRECT ), 
                    oriented_edge_descriptor( edgePivot, ign::graph::REVERSE ) 
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

                        if( graph.origins(nextEdge.descriptor) != vOriginsPivot ) break;

                        if( i == 0 )
                            path.push_back( nextEdge );
                        else
                            path.push_front( epg::graph::tools::reverse( nextEdge ) );

					    vTarget = graph.target( nextEdge );

					    if( graph.degree( vTarget ) != 2 ) break;
                    }
                    if (isLoop) break;
                }

                //--
                std::map<std::string, std::set<edge_descriptor>>::iterator mitFront = mFeatMergedEdges.find(vOriginsPivot.front());
                if( mitFront == mFeatMergedEdges.end() )
                    mitFront = mFeatMergedEdges.insert(std::make_pair(vOriginsPivot.front(), std::set<edge_descriptor>())).first;

                std::map<std::string, std::set<edge_descriptor>>::iterator mitBack = mFeatMergedEdges.find(vOriginsPivot.back());
                if( mitBack == mFeatMergedEdges.end() )
                    mitBack = mFeatMergedEdges.insert(std::make_pair(vOriginsPivot.back(), std::set<edge_descriptor>())).first;

                for (edges_path_const_iterator pit = path.begin() ; pit != path.end() ; ++pit) {
                    mitFront->second.insert(pit->descriptor);
                    mitBack->second.insert(pit->descriptor);
                    sMergedEdges.insert(pit->descriptor);
                }

                //--
                ign::feature::Feature featFront;
                _fsEdge->getFeatureById(vOriginsPivot.front(), featFront);
                std::string wTagNameFront = featFront.getAttribute(wTagName).toString();
                std::set<std::string> sMergedFront;
                epg::tools::StringTools::Split(wTagNameFront, "#", sMergedFront);

                ign::feature::Feature featBack;
                _fsEdge->getFeatureById(vOriginsPivot.back(), featBack);
                std::string wTagNameBack = featBack.getAttribute(wTagName).toString();
                std::set<std::string> sMergedBack;
                epg::tools::StringTools::Split(wTagNameBack, "#", sMergedBack);

                std::string commonFeatId = "";
                for( std::set<std::string>::const_iterator sit = sMergedFront.begin() ; sit != sMergedFront.end() ; ++sit ) {
                    if( sMergedBack.find(*sit) != sMergedBack.end() ) {
                        commonFeatId = commonFeatId == "" ? *sit : "";
                    }
                }

                ign::feature::Feature newFeat;
                if(commonFeatId != "")
                    _fsEdge->getFeatureById(commonFeatId, newFeat);
                else
                    newFeat = featFront;

                ign::geometry::LineString pathGeom = _convertPathToLineString(graph, "", path);
                newFeat.setGeometry(pathGeom);

                _fsEdge->createFeature(newFeat);
            }

            for( std::map<std::string, std::set<edge_descriptor>>::const_iterator mit = mFeatMergedEdges.begin() ; mit != mFeatMergedEdges.end() ; ++mit )
            {
                ign::feature::Feature featOrigin;

                //DEBUG
                _logger->log(epg::log::DEBUG, mit->first);

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

                for( size_t i = 0 ; i < vPaths.size() ; ++i )
                {
                    ign::geometry::LineString newGeom = _convertPathToLineString(graph, originCountry, vPaths[i]);

                    // assurer la cohérence avec les données sources aux extrémités
                    bool isStartingPath = graph.source(*vPaths[i].begin()) ==  graph.source(*foundInducedEdges.second.begin());
                    bool isEndingPath = graph.target(*vPaths[i].rbegin()) ==  graph.target(*foundInducedEdges.second.rbegin());
                    if (isStartingPath)
                        newGeom.startPoint() = geomOrigin.startPoint();
                    if (isEndingPath)
                        newGeom.endPoint() = geomOrigin.endPoint();

                    if( newGeom.isEmpty() || newGeom.isNull() || !newGeom.isValid()) {
                        _logger->log(epg::log::ERROR, "Resulting inconsistent geometry [orign feature id]: "+mit->first);
                        continue;
                    }

                    featOrigin.setGeometry(newGeom);
                    featOrigin.setAttribute(wTagName, ign::data::String(mit->first));
                    _fsEdge->createFeature(featOrigin);

                }

                _fsEdge->deleteFeature(mit->first);
            }
        }

        ///
        ///
        ///
        void CLInAreaGenerationOp::_mergeByWTag() const {
            //--
            app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getParameter(W_TAG_NAME).getValue().toString();

            //--
            detail::EdgeCleaningGraphManager graphManager;
            ign::feature::FeatureFilter filter(wTagName+" IS NOT NULL");
            _loadGraph(graphManager, filter, false);
            GraphType const& graph = graphManager.getGraph();

            //patience
            boost::progress_display display( graph.numEdges() , std::cout, "[ merging CL % complete ]\n") ;

            std::set< edge_descriptor > sMergedEdges;

            edge_iterator eit, eend ;
            for( graph.edges( eit, eend ) ; eit != eend ; ++eit )
            {
                ++display;
                // if( !graphManager.isCl(*eit)) continue;
                if( graph.target(*eit) == graph.source(*eit) ) continue;
                if( sMergedEdges.find(*eit) != sMergedEdges.end() ) continue ;

                edge_descriptor edgePivot = *eit ;

                oriented_edge_descriptor tPivot[] = { 
                    oriented_edge_descriptor( *eit, ign::graph::DIRECT ), 
                    oriented_edge_descriptor( *eit, ign::graph::REVERSE ) 
                } ;

                ign::feature::Feature featPivot;
                _fsEdge->getFeatureById(graph.origins(*eit)[0], featPivot);
                std::string wTagPivot = featPivot.getAttribute(wTagName).toString();

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

                        ign::feature::Feature featNext;
                        _fsEdge->getFeatureById(graph.origins(nextEdge.descriptor)[0], featNext);
                        std::string wTagNext = featNext.getAttribute(wTagName).toString();

                        if( wTagNext != wTagPivot ) break;

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
        bool CLInAreaGenerationOp::_isFictitious(
            GraphType const& graph,
            std::string const& country,
            std::list<oriented_edge_descriptor> const& path
        ) const {
            //--
            params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
            double const ficticiousRatioThreshold = themeParameters->getValue( CLA_FICTITIOUS_RATIO_THRESHOLD ).toDouble();
            double const ficticiousLengthThreshold = themeParameters->getValue( CLA_FICTITIOUS_LENGTH_THRESHOLD ).toDouble();
            
            //--
            std::pair<double, double> pRatioLength = _getFictitiousRatio(graph, country, path);

            return pRatioLength.first > ficticiousRatioThreshold || pRatioLength.second > ficticiousLengthThreshold;
        }

        ///
        ///
        ///
        std::pair<double, double> CLInAreaGenerationOp::_getFictitiousRatio(
            GraphType const& graph,
            std::string const& country,
            std::list<oriented_edge_descriptor> const& path
        ) const {
            //--
            double length = 0;
            double fictitiousLength = 0;

            //--
            params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
            std::string const fictitiousFieldName = themeParameters->getValue(EDGE_FICTITIOUS_NAME).toString();

            std::string originId = _getOrigin(graph, country, path.begin()->descriptor);

            ign::feature::Feature featOrigin;
            _fsEdge->getFeatureById(originId, featOrigin);

            std::string fictitious = featOrigin.getAttribute(fictitiousFieldName).toString();
            
            edges_path_const_iterator pit;
            for( pit = path.begin() ; pit != path.end(); ++pit )
            {
                ign::geometry::LineString edgeGeom = graph.getGeometry(pit->descriptor);
                double edgeLength = edgeGeom.length();
                length += edgeLength;

                std::string currentOriginId = _getOrigin(graph, country, pit->descriptor);
                if( currentOriginId != originId ) {
                    originId = currentOriginId;
                    _fsEdge->getFeatureById(originId, featOrigin);

                    fictitious = featOrigin.getAttribute(fictitiousFieldName).toString();
                }

                if ( fictitious == "true" ) {
                    fictitiousLength += edgeLength;
                }
                
            }
            return std::make_pair(fictitiousLength / length, fictitiousLength);
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
		std::map<size_t, std::string> CLInAreaGenerationOp::_getOriginEdges(
			GraphType const& graph,
            std::string const& country,
			std::list<oriented_edge_descriptor> const& path,
            std::vector<int> const& pathgeomMeanGeomLink,
			std::map<std::string, std::set<edge_descriptor>> & mFeatMergedEdges,
            std::map<size_t, std::vector<detail::IncidentFeature>> & mIncidentFeatures
		) const {
            std::map<size_t, std::string> sMeanIndexEdge;

            ign::geometry::LineString ls = _convertPathToLineString(graph, country, path);
            // double length = ls.length();
            // geometry::tools::LengthIndexedLineString lsIndex(ls);

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
                    size_t meanIndex = pathgeomMeanGeomLink[i];
                    // double abs = lsIndex.getPointLocation(i)/length;

                    sMeanIndexEdge.insert(std::make_pair(meanIndex, previousFeat));

                    //on recupere les edges incidents
                    if( graph.degree(graph.source(*lit)) > 2 ) {
                        ign::geometry::Point vertexGeom = graph.getGeometry(graph.source(*lit));
                        std::vector< edge_descriptor > vEdges = graph.incidentEdges(graph.source(*lit));

                        for ( std::vector<edge_descriptor>::const_iterator vit = vEdges.begin() ; vit != vEdges.end() ; ++vit ) {

                            std::string incidentFeatId = graph.origins(*vit)[0];
                            if (incidentFeatId == currentFeat || incidentFeatId == previousFeat) continue;

                            std::map<size_t, std::vector<detail::IncidentFeature>>::iterator mit = mIncidentFeatures.find(meanIndex);
                            if( mit == mIncidentFeatures.end() )
                                mit = mIncidentFeatures.insert(std::make_pair(meanIndex, std::vector<detail::IncidentFeature>())).first;

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
            sMeanIndexEdge.insert(std::make_pair(pathgeomMeanGeomLink.back(), previousFeat));

			return sMeanIndexEdge;
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
			std::string const wTagName = themeParameters->getParameter(W_TAG_NAME).getValue().toString();
            std::string const natIdName = themeParameters->getParameter(NATIONAL_IDENTIFIER_NAME).getValue().toString();

            // chargement des edges
            size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filter);
            boost::progress_display display(numFeatures, std::cout, "[ edge_loading  % complete ]\n");

            ign::feature::FeatureIteratorPtr itEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filter);
            while (itEdge->hasNext())
            {
                ++display;
                ign::feature::Feature const& fEdge = itEdge->next();
                ign::geometry::LineString const& ls = fEdge.getGeometry().asLineString();
                std::string edgeId = fEdge.getId();
                std::string country = fEdge.getAttribute(countryCodeName).toString();
                bool isCl = country.find("#") != std::string::npos;

                OriginEdgeProperties edgeProp(country, isCl);
                if (planarize) {
                    edgeProp.natId = fEdge.getAttribute(natIdName).toString();
                    graphManager.addEdge(ls, edgeId, edgeProp);
                } else {
                    edgeProp.wTag = fEdge.getAttribute(wTagName).toString();
                    graphManager.addEdgeSimple(ls, edgeId, edgeProp);
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
        double CLInAreaGenerationOp::_pathsGeomAreEqual(
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
                if( hausdorffDist >= 0 && hausdorffDist < maxWidth )
                    return 1;
                else
                    return -1;
            }

            //--
            double lengthRation = 0.66666;
            double length1 = path1geom.length();
            double length2 = path2geom.length();
            bool bLength = length1 > lengthRation*length2 && length2 > lengthRation*length1;

            //--
            double meanWidth = 2 * ( poly.area() / poly.exteriorRing().length() );
            bool bHausdorff = hausdorffDist >= 0 && (hausdorffDist < 3*maxWidth);

            if( meanWidth < maxWidth && bHausdorff && bLength)
                return meanWidth;
            else
                return -1;
        }

		///
        ///
        ///
        ign::geometry::LineString CLInAreaGenerationOp::_computeMeanPath(
            ign::geometry::LineString const& path1geom,
            ign::geometry::LineString const& path2geom,
            std::vector<int> & path1geomMeanGeomLink,
            std::vector<int> & path2geomMeanGeomLink
		) const {
            path1geomMeanGeomLink.resize(path1geom.numPoints());
            path2geomMeanGeomLink.resize(path2geom.numPoints());

			ign::geometry::LineString meanGeom;

			geometry::tools::LengthIndexedLineString lsIndex1(path1geom);
			geometry::tools::LengthIndexedLineString lsIndex2(path2geom);

            std::set<double> sAbsCurv1;
			double l1 = path1geom.length();
			for (size_t i = 0; i < path1geom.numPoints(); ++i) {
				double abscurv = lsIndex1.getPointLocation(i)/l1;
				sAbsCurv1.insert(abscurv);
			}
            std::set<double> sAbsCurv2;
			double l2 = path2geom.length();
			for (size_t i = 0; i < path2geom.numPoints(); ++i) {
				double abscurv = lsIndex2.getPointLocation(i)/l2;
				sAbsCurv2.insert(abscurv);
			}

            std::set<double> sAbsCurv = sAbsCurv1;
            sAbsCurv.insert(sAbsCurv2.begin(), sAbsCurv2.end());
			for (std::set<double>::iterator sit = sAbsCurv.begin(); sit != sAbsCurv.end(); ++sit) {
                std::set<double>::const_iterator sit1 = sAbsCurv1.find(*sit);
                if( sit1 != sAbsCurv1.end() )
                    path1geomMeanGeomLink[std::distance(sAbsCurv1.begin(), sit1)]= std::distance(sAbsCurv.begin(), sit);

                std::set<double>::const_iterator sit2 = sAbsCurv2.find(*sit);
                if( sit2 != sAbsCurv2.end() )
                    path2geomMeanGeomLink[std::distance(sAbsCurv2.begin(), sit2)]= std::distance(sAbsCurv.begin(), sit);

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
        // void CLInAreaGenerationOp::_mergeCl(
        //     detail::EdgeCleaningGraphManager & graphManager
        // ) const {
        //      GraphType & graph = graphManager.getGraph();

        //     bool mergingOccured = false;
        //     do{
        //         edge_iterator eit, eend ;
        //         for( graph.edges( eit, eend ) ; eit != eend ; ++eit )
        //         {
        //             if( !graphManager.isCl(*eit) ) continue;

        //             std::string clNatId = graphManager.getNatId(*eit);

        //             std::vector<edge_descriptor> vEdgeToMerge;
        //             vertex_descriptor vConnect;

        //             // on regarde si edge mergeable cote source
        //             std::string incidentNatIdSource = "";
        //             if( graph.degree(graph.source(*eit)) == 2 ) {
        //                 std::vector< edge_descriptor > vEdges = graph.incidentEdges(graph.source(*eit));
        //                 edge_descriptor incidentEdge = vEdges.front() == *eit ? vEdges.back() : vEdges.front();
        //                 if ( !graphManager.isCl(incidentEdge) ) {
        //                     std::string incidentNatIdSource = graphManager.getNatId(incidentEdge);
        //                     if (clNatId.find(incidentNatIdSource) != std::string::npos) {
        //                         vEdgeToMerge.push_back(incidentEdge);
        //                         vConnect = graph.source(*eit);
        //                     }
        //                 }
        //             }

        //             // on regarde si edge mergeable cote target
        //             std::string incidentNatIdTarget = "";
        //             if( graph.degree(graph.target(*eit)) == 2 ) {
        //                 std::vector< edge_descriptor > vEdges = graph.incidentEdges(graph.target(*eit));
        //                 edge_descriptor incidentEdge = vEdges.front() == *eit ? vEdges.back() : vEdges.front();
        //                 if ( !graphManager.isCl(incidentEdge) ) {
        //                     std::string incidentNatIdTarget = graphManager.getNatId(incidentEdge);
        //                     if (clNatId.find(incidentNatIdTarget) != std::string::npos) {
        //                         vEdgeToMerge.push_back(incidentEdge);
        //                         vConnect = graph.target(*eit);
        //                     }
        //                 }
        //             }

        //             if( vEdgeToMerge.size() != 1 ) continue;
        //             //todo : voir si on merge si (vEdgeToMerge.size() ==2 && incidentNatIdSource == incidentNatIdTarget)

        //             mergingOccured = true;

        //             edge_descriptor eIncident = vEdgeToMerge.front();

        //             //fusion des geometries
        //             ign::geometry::LineString clGeom = graph.getGeometry(*eit);
        //             ign::geometry::LineString incidentGeom = graph.getGeometry(eIncident);
                    
        //             ign::graph::EdgeDirection incidentDirection = graph.source(eIncident) == vConnect ? ign::graph::DIRECT : ign::graph::REVERSE;
        //             ign::graph::EdgeDirection clDirection = graph.source(*eit) == vConnect ? ign::graph::DIRECT : ign::graph::REVERSE;

        //             if( incidentDirection == clDirection )
        //                 clGeom.reverse();

        //             //suppression cl
        //             graphManager.removeEdge(*eit);

        //             //remplacement edge incident
        //             vertex_descriptor vNewConnect = clDirection == ign::graph::DIRECT ? graph.target(*eit) : graph.source(*eit);

        //             vertex_descriptor vSource = incidentDirection == ign::graph::DIRECT ? vNewConnect : graph.source(eIncident);
        //             vertex_descriptor vTarget = incidentDirection == ign::graph::DIRECT ? graph.target(eIncident) : vNewConnect;

        //             std::vector<ign::geometry::Point> intermediatePoints = incidentDirection == ign::graph::DIRECT ? _getIntermediatePoints(clGeom, incidentGeom) : _getIntermediatePoints(incidentGeom, clGeom);

        //             std::vector<std::string> vOrigins = graph.origins(eIncident);
        //             oriented_edge_descriptor oeNewIncident = graph.addEdge(vSource, vTarget, intermediatePoints);

        //             for (std::vector<std::string>::const_iterator vit = vOrigins.begin() ; vit != vOrigins.end() ; ++vit) {
        //                 std::pair<bool, std::vector<oriented_edge_descriptor>> foundInducedEdges = graph.getInducedEdges(*vit);
        //                 if (foundInducedEdges.first) {
        //                     for (typename std::vector<oriented_edge_descriptor>::iterator vit2 = foundInducedEdges.second.begin() ; vit2 != foundInducedEdges.second.end() ; ++vit2) {
        //                         if (vit2->descriptor == eIncident) {
        //                             vit2->descriptor = oeNewIncident.descriptor;
        //                         }
        //                     }
        //                     graph.setInducedEdges(*vit, foundInducedEdges.second);
        //                     graph.addOrigin( oeNewIncident.descriptor , *vit );
        //                 }
        //             }
        //             graph.setOrigins( eIncident , std::vector<std::string>() );
        //             graph.removeEdge( eIncident, true );
        //         }
        //     } while ( mergingOccured );
        // }

        ///
        ///
        ///
        // std::vector<ign::geometry::Point> CLInAreaGenerationOp::_getIntermediatePoints(
        //     ign::geometry::LineString const& lsFront, 
        //     ign::geometry::LineString const& lsBack
        // ) const {
        //     std::vector<ign::geometry::Point> vPoints;
        //     for( size_t i = 1 ; i < lsFront.numPoints() ; ++i ) {
        //         vPoints.push_back(lsFront.pointN(i));
        //     }
        //     for( size_t i = 1 ; i < lsBack.numPoints()-1 ; ++i ) {
        //         vPoints.push_back(lsBack.pointN(i));
        //     }
        //     return vPoints;
        // }

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
                        currentCountry = nextCountry;
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