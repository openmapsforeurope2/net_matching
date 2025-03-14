
//APP
#include <app/calcul/JunctionMatchingOp.h>
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LineStringSplitter.h>

//BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

//SOCLE
#include <ign/geometry/graph/builder/SimpleGraphBuilder.h>

//EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/tools/StringTools.h>
#include <epg/utils/replaceTableName.h>
#include <epg/tools/geometry/angle.h>
#include <epg/sql/tools/numFeatures.h>

///
///
///
app::calcul::JunctionMatchingOp::JunctionMatchingOp(
	std::string const& countryCodeDouble,
	bool verbose
){
	_init(countryCodeDouble, verbose);
}

///
///
///
app::calcul::JunctionMatchingOp::~JunctionMatchingOp()
{

}

///
///
///
void app::calcul::JunctionMatchingOp::MatchJunctions(
	std::string const& countryCodeDouble,
	bool verbose
){
	JunctionMatchingOp op(countryCodeDouble, verbose);
	op._matchJunctions();
}

///
///
///
void app::calcul::JunctionMatchingOp::_init(
	std::string const& countryCodeDouble,
	bool verbose
) {
	_logger = epg::log::EpgLoggerS::getInstance();
	_logger->log(epg::log::TITLE, "[ BEGIN INITIALIZATION ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const edgeTableName = context->getEpgParameters().getValue(EDGE_TABLE).toString();

	_countryCodeDouble = countryCodeDouble;
	epg::tools::StringTools::Split(_countryCodeDouble, "#", _vCountriesCodeName);
	_verbose = verbose;


	///recuperation des features
	std::string const boundaryTableName = epg::utils::replaceTableName(context->getEpgParameters().getValue(TARGET_BOUNDARY_TABLE).toString());
	_fsBoundary = context->getDataBaseManager().getFeatureStore(boundaryTableName, idName, geomName);

	std::string const landmaskTableName = epg::utils::replaceTableName(themeParameters->getValue(LANDMASK_TABLE).toString());
	_fsLandmask = context->getDataBaseManager().getFeatureStore(landmaskTableName, idName, geomName);

	_fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);


	_logger->log(epg::log::TITLE, "[ END INITIALIZATION ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::JunctionMatchingOp::_matchJunctions() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN MATCHING JUNCTIONS " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	//--
	app::calcul::detail::EdgeCleaningGraphManager graphManager1, graphManager2;
	_loadGraph(_vCountriesCodeName[0], graphManager1);
	_loadGraph(_vCountriesCodeName[1], graphManager2);

	GraphType const& graphEdgCountry1 = graphManager1.getGraph();
	GraphType const& graphEdgCountry2 = graphManager2.getGraph();

	//--
	std::map< vertex_descriptor, vertex_descriptor> mMatchedJ1WithBestJ2;
	_getMatchedJunctBest(mMatchedJ1WithBestJ2, graphEdgCountry1, graphEdgCountry2);

	std::map < vertex_descriptor, vertex_descriptor> mMatchedJ2BestWithJ1;
	_getMatchedJunctBest(mMatchedJ2BestWithJ1, graphEdgCountry2, graphEdgCountry1);

	//compare mMatchedBestJ1 et mMatchedBestJ2
	//verification que les "meilleur noeuds" dans country1 et dans country 2 sont réciproques
	boost::progress_display display(mMatchedJ1WithBestJ2.size(), std::cout, "[ MATCH BEST JUNCTIONS CANDIDATES]\n");
	
	std::map<std::string, ign::feature::Feature> mEdgesModifiedGeom;

	for (std::map< vertex_descriptor, vertex_descriptor>::const_iterator mitJ1 = mMatchedJ1WithBestJ2.begin();
		mitJ1 != mMatchedJ1WithBestJ2.end(); ++mitJ1) {
		++display;

		vertex_descriptor vJ1MatchedWithJ2 = mitJ1->first;
		vertex_descriptor vJ2BestCandidateFromJ1 = mitJ1->second;

		std::map< vertex_descriptor, vertex_descriptor>::const_iterator mit2 = mMatchedJ2BestWithJ1.find(vJ2BestCandidateFromJ1);
		//si le candidat J2 n'a pas de meilleur résultat associé dans 1
		if (mit2 == mMatchedJ2BestWithJ1.end())
			continue;

		//si le meilleur candidat n'est pas réciproque (J2 associé à un autre J1' que le J1 qui s'associe à lui)
		vertex_descriptor vJ1BestCandidateFromJ2 = mit2->second;
		if (vJ1BestCandidateFromJ2 != vJ1MatchedWithJ2)
			continue;

		//si ce sont les mêmes meilleurs candidats réciproques, on modifie la geom des carrefour
		//modification des edges lié à ce point dans J1 et J2
		ign::geometry::Point ptJ1 = graphEdgCountry1.getGeometry(vJ1BestCandidateFromJ2);
		ign::geometry::Point ptJ2 = graphEdgCountry2.getGeometry(vJ2BestCandidateFromJ1);

		bool isFictitious1 = _isFictitious(vJ1BestCandidateFromJ2, graphManager1);
		bool isFictitious2 = _isFictitious(vJ2BestCandidateFromJ1, graphManager2);

		if( isFictitious1 && !isFictitious2 ) {
			_setNewGeomJunction(graphEdgCountry2, vJ2BestCandidateFromJ1, ptJ1, mEdgesModifiedGeom);
		} else if ( isFictitious1 && !isFictitious2 ) {
			_setNewGeomJunction(graphEdgCountry1, vJ1BestCandidateFromJ2, ptJ2, mEdgesModifiedGeom);
		} else {
			ign::geometry::MultiPoint mpJunctions2match;
			mpJunctions2match.addGeometry(ptJ1);
			mpJunctions2match.addGeometry(ptJ2);

			ign::geometry::Point ptJNew = mpJunctions2match.getCentroid();
			ptJNew.setFillZ((ptJ1.z() + ptJ2.z())*0.5);

			_setNewGeomJunction(graphEdgCountry1, vJ1BestCandidateFromJ2, ptJNew, mEdgesModifiedGeom);
			_setNewGeomJunction(graphEdgCountry2, vJ2BestCandidateFromJ1, ptJNew, mEdgesModifiedGeom);
		}
	}
		
	for (std::map<std::string, ign::feature::Feature>::iterator mit = mEdgesModifiedGeom.begin(); mit != mEdgesModifiedGeom.end(); ++mit) 
		_fsEdge->modifyFeature(mit->second);
	

	_logger->log(epg::log::TITLE, "[ END CL MERGING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
bool app::calcul::JunctionMatchingOp::_isFictitious(
	vertex_descriptor vJunction,
	app::calcul::detail::EdgeCleaningGraphManager const& graphManager
) const {
    GraphType const& graph = graphManager.getGraph();

	std::vector< edge_descriptor > vIncidentEdges = graph.incidentEdges( vJunction );

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
	graphManager.clear();

	//--
	epg::Context* context = epg::ContextS::getInstance();
	epg::params::EpgParameters const& epgParams = context->getEpgParameters();
	std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

	//--
	app::params::ThemeParameters *themeParameters = params::ThemeParametersS::getInstance();
	std::string const fictitiousFieldName = themeParameters->getValue(EDGE_FICTITIOUS).toString();
	
	//--
	ign::feature::FeatureFilter filter(countryCodeName + " = '" + country + "'");

	//--
	int numFeatures = epg::sql::tools::numFeatures(*_fsEdge, filter);
	boost::progress_display display(numFeatures, std::cout, "[ LOAD GRAPH EDGE " + country + " ]\n");

	ign::feature::FeatureIteratorPtr itEdge = _fsEdge->getFeatures(filter);
	while (itEdge->hasNext())
	{
		++display;
		ign::feature::Feature const& fEdge = itEdge->next();
		ign::geometry::LineString const& ls = fEdge.getGeometry().asLineString();
		std::string edgeId = fEdge.getId();
		std::string fictitious = fEdge.getAttribute(fictitiousFieldName).toString();

		graphManager.addEdgeSimple(ls, edgeId, OriginEdgeProperties(country, fictitious, false));
	}
}

///
///
///
void app::calcul::JunctionMatchingOp::_getMatchedJunctBest(
	std::map< vertex_descriptor,vertex_descriptor> & mMatchedJuncRefWithBestJuncMatched,
	GraphType const& graphRef,
	GraphType const& graph2match
) const {
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distMaxJunctions = themeParameters->getValue(DIST_MAX_JUNCTIONS).toDouble();


	GraphType::vertex_iterator vitRef, vitEndRef;
	graphRef.vertices(vitRef, vitEndRef);

	boost::progress_display display(graphRef.numVertices(), std::cout, "[ MATCH JUNCTIONS BY COUNTRY ]\n");
	while (vitRef != vitEndRef) {
		++display;

		//on s'assure que le noeud est un carrefour de degre au moins 3
		size_t degreeJRef = graphRef.degree(*vitRef);
		if (degreeJRef < 3) {
			++vitRef;
			continue;
		}

		ign::geometry::Point ptJRef = graphRef.getGeometry(*vitRef);

		//recuperation des noeuds proches du country2 
		ign::geometry::Envelope envlpDistMaxJunctRef(ptJRef);
		envlpDistMaxJunctRef.expandBy(distMaxJunctions);
		std::set<vertex_descriptor> sVCandidate2matchArroundJRef;
		graph2match.getVertices(envlpDistMaxJunctRef, sVCandidate2matchArroundJRef);

		double distMin = distMaxJunctions;
		for (std::set<vertex_descriptor>::const_iterator sitV2match = sVCandidate2matchArroundJRef.begin(); sitV2match != sVCandidate2matchArroundJRef.end(); ++sitV2match) {
			//on verifie le degree des noeuds
			size_t degreeCandidateJ2match = graph2match.degree(*sitV2match);
			if (degreeJRef != degreeCandidateJ2match)
				continue;
			ign::geometry::Point ptCandidateJ2match = graph2match.getGeometry(*sitV2match);
			//on verifie la distance entre les noeuds
			double distJ1J2Candidate = ptCandidateJ2match.distance(ptJRef);
			if (distJ1J2Candidate > distMaxJunctions)
				continue;

			//ajout d'une note modulant la dist selon l'orientation des edges?
			if (distJ1J2Candidate < distMin) {
				distMin = distJ1J2Candidate;
				mMatchedJuncRefWithBestJuncMatched[*vitRef] = *sitV2match;
			}
		}
		++vitRef;
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
	double const angleMaxOrientEdgJunctions = themeParameters->getValue(ANGLE_MAX_ORIENTATION_EDGES).toDouble()* M_PI / 180;

	std::set<double>::const_iterator sit1 = sAnglEdgesJ1.begin();
	std::set<double>::const_iterator sit2 = sAnglEdgesJ2.begin();
	while (sit1 != sAnglEdgesJ1.end()) {
		double diffAngl = fabs(*sit1 - *sit2);
		++sit1;
		++sit2;

		if (diffAngl > angleMaxOrientEdgJunctions) 
			return false;	
	}
	return true;
}

///
///
///
void app::calcul::JunctionMatchingOp::_setNewGeomJunction(
	GraphType const& graph,
	vertex_descriptor vJunction,
	ign::geometry::Point const& ptNewGeomJunction,
	std::map<std::string, ign::feature::Feature> & mEdgesModifiedGeom
) const {
	std::vector< oriented_edge_descriptor > vEdgesIncidentsJunction;
	graph.incidentEdges(vJunction, vEdgesIncidentsJunction);

	for (std::vector< oriented_edge_descriptor >::const_iterator oeit = vEdgesIncidentsJunction.begin(); oeit != vEdgesIncidentsJunction.end(); ++oeit) {
		//recuperation du edge feature associe et modication de la geom de l'edge dans la table
		std::string idEdge2modify = graph.origins(oeit->descriptor)[0];
		ign::feature::Feature featEdge2modify;
		if (mEdgesModifiedGeom.find(idEdge2modify) != mEdgesModifiedGeom.end())
			featEdge2modify = mEdgesModifiedGeom.find(idEdge2modify)->second;
		else
			_fsEdge->getFeatureById(idEdge2modify, featEdge2modify);
		
		//modification de la nouvelle geometrie de l'edge
		ign::geometry::LineString lsEdge2modify = featEdge2modify.getGeometry().asLineString();

		//projection de l'edge, et recup�ration de l'abs curviligne
		//suppression des points de la ls entre la proj du nouveau point et le nouveau point (start ou end selon l'orientation) 
		app::geometry::tools::LineStringSplitter lsSplitter2modify(lsEdge2modify);
		lsSplitter2modify.addCuttingGeometry(ptNewGeomJunction);
		std::vector< ign::geometry::LineString > vLs2modify = lsSplitter2modify.getSubLineStringsZ();
		if (vLs2modify.size() > 1) {
			if (oeit->direction == ign::graph::DIRECT)
				lsEdge2modify = vLs2modify[1];
			else
				lsEdge2modify = vLs2modify[0];
		}

		if (oeit->direction == ign::graph::DIRECT)
			lsEdge2modify.setPointN(ptNewGeomJunction, 0);
		else
			lsEdge2modify.setPointN(ptNewGeomJunction, lsEdge2modify.numPoints() - 1);


		featEdge2modify.setGeometry(lsEdge2modify);
		mEdgesModifiedGeom[idEdge2modify] = featEdge2modify;
	}
}