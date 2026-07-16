// APP
#include <app/calcul/JunctionMatchingOp.h>
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LineStringSplitter.h>

// BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

// EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/tools/StringTools.h>
#include <epg/tools/geometry/angle.h>
#include <epg/sql/tools/numFeatures.h>
#include <epg/tools/geometry/project.h>

// OME2
#include <ome2/feature/sql/NotDestroyedTools.h>


///
///
///
app::calcul::JunctionMatchingOp::JunctionMatchingOp(
	std::string const& borderCode,
	bool verbose 
): 
	_borderCode(borderCode),
	_verbose(verbose)
{
	_init();
}

///
///
///
app::calcul::JunctionMatchingOp::~JunctionMatchingOp()
{
	_shapeLogger->closeShape("displacements");
	_shapeLogger->closeShape("projections");
	_shapeLogger->closeShape("projected_points");
	_shapeLogger->closeShape("split_edges");
}

///
///
///
void app::calcul::JunctionMatchingOp::Compute(
	std::string const& borderCode,
	bool verbose
){
	JunctionMatchingOp op(borderCode, verbose);
	op._matchJunctions();
}

///
///
///
void app::calcul::JunctionMatchingOp::_init()
{
	_logger = epg::log::EpgLoggerS::getInstance();
	_logger->log(epg::log::TITLE, "[ BEGIN INITIALIZATION ] : " + epg::tools::TimeTools::getTime());

	//--
	epg::Context* context = epg::ContextS::getInstance();

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const edgeTableName = context->getEpgParameters().getValue(EDGE_TABLE).toString();

	//--
	_shapeLogger = epg::log::ShapeLoggerS::getInstance();
    _shapeLogger->addShape("displacements", epg::log::ShapeLogger::LINESTRING);
    _shapeLogger->addShape("projections", epg::log::ShapeLogger::LINESTRING);
    _shapeLogger->addShape("split_edges", epg::log::ShapeLogger::LINESTRING);
    _shapeLogger->addShape("projected_points", epg::log::ShapeLogger::POINT);

	//--
	epg::tools::StringTools::Split(_borderCode, "#", _vCountry);

	//--
	_fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);


	_logger->log(epg::log::TITLE, "[ END INITIALIZATION ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::JunctionMatchingOp::_matchJunctions() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN MATCHING JUNCTIONS " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

	//--
	app::calcul::detail::EdgeCleaningGraphManager graphManager1, graphManager2;
	_loadGraph(_vCountry[0], graphManager1);
	_loadGraph(_vCountry[1], graphManager2);

	GraphType const& graph1 = graphManager1.getGraph();
	GraphType const& graph2 = graphManager2.getGraph();

	//--
	std::map< vertex_descriptor, vertex_descriptor> mJ1J2;
	_computeOrientedMatching(mJ1J2, graph1, graph2);

	std::map < vertex_descriptor, vertex_descriptor> mJ2J1;
	_computeOrientedMatching(mJ2J1, graph2, graph1);

	//compare mJ1J2 et mJ2J1
	//verification que les "meilleur noeuds" dans country1 et dans country 2 sont réciproques
	boost::progress_display display(mJ1J2.size(), std::cout, "[ matching best junction candidates % complete ]\n");
	
	std::map<std::string, ign::feature::Feature> mEdges2Modify;

	for (std::map< vertex_descriptor, vertex_descriptor>::const_iterator mitJ1J2 = mJ1J2.begin(); mitJ1J2 != mJ1J2.end(); ++mitJ1J2)
	{
		++display;

		std::map< vertex_descriptor, vertex_descriptor>::const_iterator mitJ2J1 = mJ2J1.find(mitJ1J2->second);
		//si le candidat J2 n'a pas de meilleur résultat associé dans 1
		if (mitJ2J1 == mJ2J1.end())
			continue;

		//si le meilleur candidat n'est pas réciproque (J2 associé à un autre J1' que le J1 qui s'associe à lui)
		vertex_descriptor vJ1BestCandidateFromJ2 = mitJ2J1->second;
		if (mitJ2J1->second != mitJ1J2->first)
			continue;

		//si ce sont les mêmes meilleurs candidats réciproques, on modifie la geom des carrefour
		//modification des edges lié à ce point dans J1 et J2
		ign::geometry::Point ptJ1 = graph1.getGeometry(mitJ2J1->second);
		ign::geometry::Point ptJ2 = graph2.getGeometry(mitJ1J2->second);

		bool isFictitious1 = _isFictitious(mitJ2J1->second, graphManager1);
		bool isFictitious2 = _isFictitious(mitJ1J2->second, graphManager2);

		if( isFictitious1 && !isFictitious2 ) {
			_moveVertex(graph2, mitJ1J2->second, ptJ1, mEdges2Modify);
		} else if ( !isFictitious1 && isFictitious2 ) {
			_moveVertex(graph1, mitJ2J1->second, ptJ2, mEdges2Modify);
		} else {
			ign::geometry::MultiPoint mpt;
			mpt.addGeometry(ptJ1);
			mpt.addGeometry(ptJ2);

			ign::geometry::Point ptJNew = mpt.getCentroid();
			ptJNew.setFillZ((ptJ1.z() + ptJ2.z())*0.5);

			_moveVertex(graph1, mitJ2J1->second, ptJNew, mEdges2Modify);
			_moveVertex(graph2, mitJ1J2->second, ptJNew, mEdges2Modify);
		}
	}
		
	for (std::map<std::string, ign::feature::Feature>::iterator mit = mEdges2Modify.begin(); mit != mEdges2Modify.end(); ++mit) 
		_fsEdge->modifyFeature(mit->second);
	

	_logger->log(epg::log::TITLE, "[ END CL MERGING FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
bool app::calcul::JunctionMatchingOp::_isFictitious(
	vertex_descriptor v,
	app::calcul::detail::EdgeCleaningGraphManager const& graphManager
) const {
    GraphType const& graph = graphManager.getGraph();

	std::vector< edge_descriptor > vIncidentEdges = graph.incidentEdges( v );

	typename std::vector< edge_descriptor >::const_iterator eit;
	for( eit = vIncidentEdges.begin() ; eit != vIncidentEdges.end() ; ++eit )
	{
		if ( graphManager.getWTag(*eit) == "true" ) return true;
	}
	return false;
}

///
///
///
void app::calcul::JunctionMatchingOp::_loadGraph(
	std::string const& country,
	app::calcul::detail::EdgeCleaningGraphManager & graphManager
) const {
	//--
	graphManager.clear();

	//--
	epg::Context* context = epg::ContextS::getInstance();
	epg::params::EpgParameters const& epgParams = context->getEpgParameters();
	std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

	//--
	app::params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
	std::string const fictitiousFieldName = themeParameters->getValue(EDGE_FICTITIOUS_NAME).toString();
	
	//--
	ign::feature::FeatureFilter filter(countryCodeName + " LIKE '%" + country + "%'");

	//--
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filter);
	boost::progress_display display(numFeatures, std::cout, "[ loading graph " + country + " % complete ]\n");

	ign::feature::FeatureIteratorPtr itEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filter);
	while (itEdge->hasNext())
	{
		++display;
		ign::feature::Feature fEdge = itEdge->next();
		ign::geometry::LineString const& ls = fEdge.getGeometry().asLineString();
		std::string edgeId = fEdge.getId();
		std::string fictitious = fEdge.getAttribute(fictitiousFieldName).toString();

		graphManager.addEdgeSimple(ls, edgeId, OriginEdgeProperties(country, fictitious, false));
	}
}

///
///
///
void app::calcul::JunctionMatchingOp::_computeOrientedMatching(
	std::map< vertex_descriptor,vertex_descriptor> & mJunct1Junct2,
	GraphType const& graph1,
	GraphType const& graph2
) const {
	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distMaxJunctions = themeParameters->getValue(JM_MAX_DIST).toDouble();

	//--
	GraphType::vertex_iterator vit1, vit1End;
	graph1.vertices(vit1, vit1End);

	//--
	boost::progress_display display(graph1.numVertices(), std::cout, "[ matching junctions by country % complete ]\n");
	while (vit1 != vit1End) {
		++display;

		//on s'assure que le noeud est un carrefour de degre au moins 3
		size_t degreeV1 = graph1.degree(*vit1);
		if (degreeV1 < 3) {
			++vit1;
			continue;
		}

		ign::geometry::Point ptV1 = graph1.getGeometry(*vit1);

		//recuperation des noeuds proches du country2 
		ign::geometry::Envelope bboxV1(ptV1);
		bboxV1.expandBy(distMaxJunctions);
		std::set<vertex_descriptor> sGraph2Candidates;
		graph2.getVertices(bboxV1, sGraph2Candidates);

		double distMin = distMaxJunctions;
		for (std::set<vertex_descriptor>::const_iterator sit2 = sGraph2Candidates.begin(); sit2 != sGraph2Candidates.end(); ++sit2) {
			//on verifie le degree des noeuds
			size_t degreeV2 = graph2.degree(*sit2);
			if (degreeV1 != degreeV2)
				continue;
			ign::geometry::Point ptV2 = graph2.getGeometry(*sit2);
			//on verifie la distance entre les noeuds
			double distV1V2 = ptV2.distance(ptV1);

			//ajout d'une note modulant la dist selon l'orientation des edges?
			if (distV1V2 < distMin) {
				distMin = distV1V2;
				mJunct1Junct2[*vit1] = *sit2;
			}
		}
		++vit1;
	}
}

///
///
///
bool app::calcul::JunctionMatchingOp::_IsSimilarIncidentsEdgesOnJunctions(
	std::set<double> const& sAnglEdgesJ1,
	std::set<double> const& sAnglEdgesJ2
) const {
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const angleMaxOrientEdgJunctions = themeParameters->getValue(JM_MAX_ANGLE).toDouble()* M_PI / 180;

	std::set<double>::const_iterator sit1 = sAnglEdgesJ1.begin();
	std::set<double>::const_iterator sit2 = sAnglEdgesJ2.begin();
	while (sit1 != sAnglEdgesJ1.end()) {
		double diffAngle = fabs(*sit1 - *sit2);
		++sit1;
		++sit2;

		if (diffAngle > angleMaxOrientEdgJunctions) 
			return false;	
	}
	return true;
}

///
///
///
void app::calcul::JunctionMatchingOp::_moveVertex(
	GraphType const& graph,
	vertex_descriptor v,
	ign::geometry::Point const& newVertexGeom,
	std::map<std::string, ign::feature::Feature> & mEdges2Modify
) const {

	std::vector< oriented_edge_descriptor > vIncidentEdges;
	graph.incidentEdges(v, vIncidentEdges);

	for (std::vector< oriented_edge_descriptor >::const_iterator oeit = vIncidentEdges.begin(); oeit != vIncidentEdges.end(); ++oeit) {
		//recuperation du edge feature associe et modication de la geom de l'edge dans la table
		std::string idEdge = graph.origins(oeit->descriptor)[0];
		ign::feature::Feature fEdge;
		if (mEdges2Modify.find(idEdge) != mEdges2Modify.end())
			fEdge = mEdges2Modify.find(idEdge)->second;
		else
			_fsEdge->getFeatureById(idEdge, fEdge);
		
		//modification de la nouvelle geometrie de l'edge
		ign::geometry::LineString lsEdge2modify = fEdge.getGeometry().asLineString();

		//projection de l'edge, et recupération de l'abs curviligne
		//suppression des points de la ls entre la proj du nouveau point et le nouveau point (start ou end selon l'orientation) 
		app::geometry::tools::LineStringSplitter lsSplitter2modify(lsEdge2modify);
		ign::geometry::Point proj = epg::tools::geometry::project(lsEdge2modify, newVertexGeom);
		lsSplitter2modify.addCuttingGeometry(proj);
		std::vector< ign::geometry::LineString > vLs2modify = lsSplitter2modify.getSubLineStringsZ();

		//DEBUG
		{
			ign::feature::Feature feat;
			feat.setGeometry(ign::geometry::LineString(newVertexGeom, proj));
			_shapeLogger->writeFeature("projections", feat);
			feat.setGeometry(proj);
			_shapeLogger->writeFeature("projected_points", feat);
		}
		
		if (vLs2modify.size() > 1) {
			if (oeit->direction == ign::graph::DIRECT)
				lsEdge2modify = vLs2modify[1];
			else
				lsEdge2modify = vLs2modify[0];
		}

		//DEBUG
		{
			ign::feature::Feature feat;
			feat.setGeometry(lsEdge2modify);
			_shapeLogger->writeFeature("split_edges", feat);
		}

		if (oeit->direction == ign::graph::DIRECT){
			//DEBUG
			ign::feature::Feature feat;
			feat.setGeometry(ign::geometry::LineString(newVertexGeom, lsEdge2modify.startPoint()));
			_shapeLogger->writeFeature("displacements", feat);

			lsEdge2modify.setPointN(newVertexGeom, 0);
		}	
		else {
			//DEBUG
			ign::feature::Feature feat;
			feat.setGeometry(ign::geometry::LineString(newVertexGeom, lsEdge2modify.endPoint()));
			_shapeLogger->writeFeature("displacements", feat);

			lsEdge2modify.setPointN(newVertexGeom, lsEdge2modify.numPoints() - 1);
		}

		fEdge.setGeometry(lsEdge2modify);
		mEdges2Modify[idEdge] = fEdge;
	}
}