
#include <app/calcul/JunctionMatchingOp.h>

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

//APP
#include <app/params/ThemeParameters.h>

///
///
///
app::calcul::JunctionMatchingOp::JunctionMatchingOp(std::string countryCodeDouble, bool verbose)
{
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
void app::calcul::JunctionMatchingOp::MatchJunctions(std::string countryCodeDouble, bool verbose)
{
	JunctionMatchingOp op(countryCodeDouble, verbose);
	op._matchJunctions();
}


void app::calcul::JunctionMatchingOp::_init(std::string countryCodeDouble, bool verbose)
{
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
void app::calcul::JunctionMatchingOp::_matchJunctions()
{
	_logger->log(epg::log::TITLE, "[ BEGIN MATCHING JUNCTIONS " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	//--

	GraphType graphEdgCountry1, graphEdgCountry2;
	_loadGraphEdges(_vCountriesCodeName[0],graphEdgCountry1);
	_loadGraphEdges(_vCountriesCodeName[1],graphEdgCountry2);

	std::map< vertex_descriptor, vertex_descriptor> mMatchedJ1WithBestJ2;
	_getMatchedJunctBest(mMatchedJ1WithBestJ2, graphEdgCountry1, graphEdgCountry2);

	std::map < vertex_descriptor, vertex_descriptor> mMatchedJ2BestWithJ1;
	_getMatchedJunctBest(mMatchedJ2BestWithJ1, graphEdgCountry2, graphEdgCountry1);

	//compare mMatchedBestJ1 et mMatchedBestJ2
	//verification que les "meilleur noeuds" dans country1 et dans country 2 sont réciproques
	boost::progress_display display(mMatchedJ1WithBestJ2.size(), std::cout, "[ MATCH BEST JUNCTIONS CANDIDATES]\n");
	
	std::map<std::string, ign::feature::Feature> mEdgesModifiedGeom;

	for (std::map< vertex_descriptor, vertex_descriptor>::iterator mitJ1 = mMatchedJ1WithBestJ2.begin();
		mitJ1 != mMatchedJ1WithBestJ2.end(); ++mitJ1) {
		vertex_descriptor vJ1MatchedWithJ2 = mitJ1->first;
		vertex_descriptor vJ2BestCandidateFromJ1 = mitJ1->second;

		//si le candidat J2 n'a pas de meilleur résultat associé dans 1
		if (mMatchedJ2BestWithJ1.find(vJ2BestCandidateFromJ1) == mMatchedJ2BestWithJ1.end())
			continue;

		vertex_descriptor vJ1BestCandidateFromJ2 = mMatchedJ2BestWithJ1.find(vJ2BestCandidateFromJ1)->second;
		//si le meilleur candidat n'est pas réciproque (J2 associé à un autre J1' que le J1 qui s'associe à lui)
		if (vJ1BestCandidateFromJ2 != vJ1MatchedWithJ2)
			continue;

		//si ce sont les mêmes meilleurs candidats réciproques, on modifie la geom des carrefour
		//modification des edges lié à ce point dans J1 et J2
		ign::geometry::Point ptJ1 = graphEdgCountry1.getGeometry(vJ1BestCandidateFromJ2);
		ign::geometry::Point ptJ2 = graphEdgCountry2.getGeometry(vJ2BestCandidateFromJ1);

		ign::geometry::MultiPoint mpJunctions2match;
		mpJunctions2match.addGeometry(ptJ1);
		mpJunctions2match.addGeometry(ptJ2);

		ign::geometry::Point ptJNew = mpJunctions2match.getCentroid();
		ptJNew.setFillZ((ptJ1.z() + ptJ2.z())*0.5);

		_setNewGeomJunction(graphEdgCountry1, vJ1BestCandidateFromJ2, ptJNew, mEdgesModifiedGeom);
		_setNewGeomJunction(graphEdgCountry2, vJ2BestCandidateFromJ1, ptJNew, mEdgesModifiedGeom);

	}
		
	for (std::map<std::string, ign::feature::Feature>::iterator mit = mEdgesModifiedGeom.begin(); mit != mEdgesModifiedGeom.end(); ++mit) 
		_fsEdge->modifyFeature(mit->second);
	

	_logger->log(epg::log::TITLE, "[ END CL MERGING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}



void app::calcul::JunctionMatchingOp::_loadGraphEdges(std::string countryCodeSimple, GraphType& graphEdges)
{
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	ign::feature::FeatureFilter filterEdgeCountryCode(countryCodeName + " = '" + countryCodeSimple + "'");
	ign::feature::FeatureIteratorPtr it = _fsEdge->getFeatures(filterEdgeCountryCode);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsEdge, filterEdgeCountryCode);
	boost::progress_display display(numFeatures, std::cout, "[ LOAD GRAPH EDGE " + countryCodeSimple + " ]\n");
	ign::geometry::graph::builder::SimpleGraphBuilder<GraphType> graphBuilder(graphEdges, 0.01);
	while (it->hasNext()) {
		++display;
		ign::feature::Feature fedge = it->next();
		graphBuilder.addEdge(fedge.getGeometry().asLineString(), fedge.getId());
	}
}



void app::calcul::JunctionMatchingOp::_getMatchedJunctBest(
	std::map< vertex_descriptor,vertex_descriptor>& mMatchedJuncRefWithBestJuncMatched,
	GraphType& graphRef,
	GraphType& graph2match)
{
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
		//DEBUG
		/*::geometry::Point ptDEbug;
		ptDEbug.setX(3940556);
		ptDEbug.setY(3004362);
		if (ptDEbug.distance(ptJRef) > 6) {
			++vitRef;
			continue;
		}*/
		//DEBUG
				
		//recupération des angles incidents
		/*std::vector< oriented_edge_descriptor > vEdgIncidenJref;
		graphRef.incidentEdges(*vitRef, vEdgIncidenJref);
		std::set<double> sAnglEdgesJRef;
		for (std::vector< oriented_edge_descriptor >::iterator eitRef = vEdgIncidenJref.begin(); eitRef != vEdgIncidenJref.end(); ++eitRef) {
			double angleEdgeIncident = _getAngleEdgeIncident(graphRef, *eitRef);
			sAnglEdgesJRef.insert(angleEdgeIncident);
		}*/

		//recuperation des noeuds proches du country2 
		ign::geometry::Envelope envlpDistMaxJunctRef(ptJRef);
		envlpDistMaxJunctRef.expandBy(distMaxJunctions);
		std::set<vertex_descriptor> sVCandidate2matchArroundJRef;
		graph2match.getVertices(envlpDistMaxJunctRef, sVCandidate2matchArroundJRef);

		/*if (sVCandidate2matchArroundJRef.size() == 1) {
			vertex_descriptor vCandidate2match = *sVCandidate2matchArroundJRef.begin();
			size_t degreeCandidateJ2match = graph2match.degree(vCandidate2match);
			if ((degreeJRef == degreeCandidateJ2match)) {
				mMatchedJuncRefWithBestJuncMatched[*vitRef] = vCandidate2match;
			}
			//si ce n'est pas le meme degré le candidat n'est pas gardé
			++vitRef;
			continue;
		}*/

		double distMin = distMaxJunctions;
		for (std::set<vertex_descriptor>::iterator sitV2match = sVCandidate2matchArroundJRef.begin(); sitV2match != sVCandidate2matchArroundJRef.end(); ++sitV2match) {
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

			//on donne une note selon l'orientation des edges
			/*std::vector< oriented_edge_descriptor > vEdgIncidenJ2match;
			graph2match.incidentEdges(*sitV2match, vEdgIncidenJ2match);
			std::set<double> sAnglEdgesJ2match;
			for (std::vector< oriented_edge_descriptor >::iterator eit2match = vEdgIncidenJ2match.begin(); eit2match != vEdgIncidenJ2match.end(); ++eit2match) {
				double angleEdgeIncident = _getAngleEdgeIncident(graph2match, *eit2match);
				sAnglEdgesJ2match.insert(angleEdgeIncident);
			}*/

			/*if ( !_IsSimilarIncidentsEdgesOnJunctions(sAnglEdgesJ1, sAnglEdgesJ2) )
				continue;*/

			//en fonction de la dist et de la note lie a l'orientation des edges ont ajoute le meilleur candidat
		}
		++vitRef;
	}
}

bool app::calcul::JunctionMatchingOp::_IsSimilarIncidentsEdgesOnJunctions(
	std::set<double> sAnglEdgesJ1,
	std::set<double> sAnglEdgesJ2)
{
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const angleMaxOrientEdgJunctions = themeParameters->getValue(ANGLE_MAX_ORIENTATION_EDGES).toDouble()* M_PI / 180;

	std::set<double>::iterator sit1 = sAnglEdgesJ1.begin();
	std::set<double>::iterator sit2 = sAnglEdgesJ2.begin();
	while (sit1 != sAnglEdgesJ1.end()) {
		double diffAngl = fabs(*sit1 - *sit2);
		++sit1;
		++sit2;

		if (diffAngl > angleMaxOrientEdgJunctions) 
			return false;	
	}
	return true;
	/*for (std::vector< oriented_edge_descriptor >::iterator eit1 = vEdgIncJ1.begin(); eit1 != vEdgIncJ1.end(); ++eit1) {
		//debug
		std::string idE1 = graphEdgCountry1.origins(eit1->descriptor)[0];

		ign::geometry::LineString ls1 = graphEdgCountry1.getGeometry(*eit1);
		if (eit1->direction == ign::graph::REVERSE)
			ls1.reverse();
		ign::math::Vec2d vec1(ls1.endPoint().x() - ls1.startPoint().x(), ls1.endPoint().y() - ls1.startPoint().y());
		double azimuth1 = atan2(vec1.y(), vec1.x());
		if (azimuth1 < 0) azimuth1 += 2 * M_PI;
		double angleEdge1 = fmod(fabs(azimuth1), 2 * M_PI);
		

		for (std::vector< oriented_edge_descriptor >::iterator eit2 = vEdgIncJ2.begin(); eit2 != vEdgIncJ2.end(); ++eit2) {
			std::string idE2 = graphEdgCountry2.origins(eit2->descriptor)[0];
			ign::geometry::LineString ls2 = graphEdgCountry2.getGeometry(*eit2);
			if (eit2->direction == ign::graph::REVERSE)
				ls2.reverse();

			ign::math::Vec2d vec2(ls2.endPoint().x() - ls2.startPoint().x(), ls2.endPoint().y() - ls2.startPoint().y());
			double angleEdge, azimuth1, azimuth2;
			azimuth1 = atan2(vec1.y(), vec1.x());
			azimuth2 = atan2(vec2.y(), vec2.x());
			
			if (azimuth2 < 0) azimuth2 += 2 * M_PI;
			angleEdge = fmod(fabs(azimuth2 - azimuth1), 2 * M_PI);
			double angleEdgeEpg = epg::tools::geometry::angle(vec1, vec2);
		}
	}*/


}
/*
double app::calcul::JunctionMatchingOp::_getAngleEdgeIncident(
	GraphType& graphEdgCountry,
	oriented_edge_descriptor& oeit)
{
	ign::geometry::LineString ls = graphEdgCountry.getGeometry(oeit.descriptor);
	if (oeit.direction == ign::graph::REVERSE)
		ls.reverse();
	ign::math::Vec2d vec(ls.endPoint().x() - ls.startPoint().x(), ls.endPoint().y() - ls.startPoint().y());
	double azimuth = atan2(vec.y(), vec.x());
	if (azimuth < 0) azimuth += 2 * M_PI;
	double angleEdgeIncident = fmod(fabs(azimuth), 2 * M_PI);
	return angleEdgeIncident;
}*/


void app::calcul::JunctionMatchingOp::_setNewGeomJunction(
	GraphType& graph,
	vertex_descriptor& vJunction,
	ign::geometry::Point ptNewGeomJunction,
	std::map<std::string, ign::feature::Feature>& mEdgesModifiedGeom)
{
	std::vector< oriented_edge_descriptor > vEdgesIncidentsJunction;
	graph.incidentEdges(vJunction, vEdgesIncidentsJunction);

	for (std::vector< oriented_edge_descriptor >::iterator oeit = vEdgesIncidentsJunction.begin(); oeit != vEdgesIncidentsJunction.end(); ++oeit) {
		//recuperation du edge feature associe et modication de la geom de l'edge dans la table
		std::string idEdge2modify = graph.origins(oeit->descriptor)[0];
		ign::feature::Feature featEdge2modify;
		if (mEdgesModifiedGeom.find(idEdge2modify) != mEdgesModifiedGeom.end())
			featEdge2modify = mEdgesModifiedGeom.find(idEdge2modify)->second;
		else
			_fsEdge->getFeatureById(idEdge2modify, featEdge2modify);
		
		//modification de la nouvelle geometrie de l'edge
		ign::geometry::LineString lsEdge2modify = featEdge2modify.getGeometry().asLineString();
		//ign::geometry::LineString lsEdge2modify = graph.getGeometry(oeit->descriptor);
		if (oeit->direction == ign::graph::DIRECT)
			lsEdge2modify.setPointN(ptNewGeomJunction, 0);
		else
			lsEdge2modify.setPointN(ptNewGeomJunction, lsEdge2modify.numPoints() - 1);

		featEdge2modify.setGeometry(lsEdge2modify);
		mEdgesModifiedGeom[idEdge2modify] = featEdge2modify;
	}
}