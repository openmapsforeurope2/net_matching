
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
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();

	GraphType graphEdgCountry1, graphEdgCountry2;
	_loadGraphEdges(_vCountriesCodeName[0],graphEdgCountry1);
	_loadGraphEdges(_vCountriesCodeName[1],graphEdgCountry2);

	GraphType::vertex_iterator vit1, vitEnd1;
	graphEdgCountry1.vertices(vit1, vitEnd1);

	boost::progress_display display1(graphEdgCountry1.numVertices(), std::cout, "[ MATCH JUNCTIONS ]\n");
	std::map<std::string, std::string> mMatchedJunctionsVertexId;

	while (vit1 != vitEnd1) {
		++display1;
		if (graphEdgCountry1.degree(*vit1) < 3) {
			++vit1;
			continue;
		}

		//recuperation des noeuds du pays 2 proche


		++vit1;
	}

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