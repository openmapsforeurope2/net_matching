// APP
#include <app/calcul/CFeatGenerationOp.h>
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LengthIndexedLineString.h>
#include <app/geometry/tools/LineStringSplitter.h>

// BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

// SOCLE
#include <ign/geometry/algorithm/BufferOpGeos.h>
#include <ign/math/Line2T.h>
#include <ign/math/LineT.h>
#include <ign/geometry/graph/builder/SimpleGraphBuilder.h>
#include <ign/geometry/graph/tools/SnapRoundPlanarizer.h>
#include <ign/geometry/algorithm/OptimizedHausdorffDistanceOp.h>
//#include <ign/geometry/tools/LengthIndexedLineString.h>
#include <ign/geometry/algorithm/OptimizedHausdorffDistanceOp.h>

// EPG
#include <epg/Context.h>
#include <epg/io/query/CountryQueriesReader.h>
#include <epg/io/query/utils/queryUtils.h>
#include <epg/tools/TimeTools.h>
#include <epg/tools/FilterTools.h>
#include <epg/tools/StringTools.h>
#include <epg/utils/CopyTableUtils.h>
#include <epg/utils/replaceTableName.h>
#include <epg/tools/geometry/interpolate.h>
#include <epg/tools/geometry/project.h>
#include <epg/tools/geometry/getSubLineString.h>
#include <epg/tools/geometry/angle.h>
#include <epg/tools/geometry/SegmentIndexedGeometry.h>
#include <epg/tools/geometry/LineIntersector.h>
#include<epg/graph/tools/merge.h>
#include<epg/graph/tools/createPath.h>
#include <epg/calcul/merging/MergingByAttributes.h>
#include <epg/graph/EpgFeatureGraph.h>

// OME2
#include <ome2/calcul/detail/ClMerger.h>

// BOOST
#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>

// QT
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>


///
///
///
app::calcul::CFeatGenerationOp::CFeatGenerationOp(std::string countryCodeDouble, bool verbose)
{
	_init(countryCodeDouble, verbose);
}

///
///
///
app::calcul::CFeatGenerationOp::~CFeatGenerationOp()
{
	_shapeLogger->closeShape("CLBeforeMerge");
	_shapeLogger->closeShape("ClMergedBeforeUpdate");
	_shapeLogger->closeShape("ClDeletedNoCandidatefound");
	_shapeLogger->closeShape("ClDoublon");
	_shapeLogger->closeShape("ClDeleteByAngleDistEdges");
	_shapeLogger->closeShape("ClDebug");
	_shapeLogger->closeShape("edgeClCutByCp");
	_shapeLogger->closeShape("lsBorderCutByAngle");

	delete _mlsBorderSmoothed;	
}

///
///
///
void app::calcul::CFeatGenerationOp::GenerateConnectingLinesByCountry(
	std::string countryCodeDouble,
	bool verbose
) {
	CFeatGenerationOp op(countryCodeDouble, verbose);
    op._generateConnectingLinesByCountry();
}

///
///
///
void app::calcul::CFeatGenerationOp::_generateConnectingLinesByCountry() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN CL GENERATION BY COUNTRY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distBuffer = themeParameters->getValue(CL_BUFFER_DIST).toDouble();
	double const thresholdNoCL = themeParameters->getValue(CL_THRESHOLD_NO_CL).toDouble();
	double const ratioInBuff = themeParameters->getValue(CL_RATIO_IN_BUFFER).toDouble();
	double const snapOnVertexBorder = themeParameters->getValue(CL_SNAP_ON_VERTEX_BORDER_DIST).toDouble();
	double angleMaxBorder = themeParameters->getValue(CL_BORDER_MAX_ANGLE).toDouble();
	double angleMaxEdges = themeParameters->getValue(CL_EDGE_MAX_ANGLE).toDouble();
	double const distMergeCL = themeParameters->getValue(CL_MERGE_CL_DIST).toDouble();
	double angleMaxToCutBorder = themeParameters->getValue(ANGLE_MAX_2_CUT_BORDER).toDouble();
	
	angleMaxBorder = angleMaxBorder * M_PI / 180;
	angleMaxEdges = angleMaxEdges * M_PI / 180;
	angleMaxToCutBorder = angleMaxToCutBorder * M_PI / 180;

	//--
	ign::feature::FeatureIteratorPtr itBoundary = _fsBoundary->getFeatures(ign::feature::FeatureFilter(countryCodeName + " = '" + _countryCodeDouble + "'"));
	while (itBoundary->hasNext())
	{
		ign::feature::Feature fBoundary = itBoundary->next();
		_logger->log(epg::log::INFO, "id boundary :"+ fBoundary.getId());

		std::string boundaryType = fBoundary.getAttribute("boundary_type").toString();
		// On ne traite que les frontières de type international_boundary ou coastline_sea_limit
		// On aurait pu filtrer en entrée mais le filtre semble très long, peut-être à cause de l'enum qui oblige à utiliser boundary_type::text like '%coastline_sea_limit%'
		if (boundaryType != "international_boundary" && boundaryType.find("coastline_sea_limit") == -1)
			continue;

		ign::geometry::LineString lsBoundary = fBoundary.getGeometry().asLineString();
		std::vector<ign::geometry::LineString> vLsBorderCutByAngle;
		_getBorderCutByAngle(lsBoundary, vLsBorderCutByAngle, angleMaxToCutBorder);
		for (size_t i = 0; i < vLsBorderCutByAngle.size(); ++i) {
			ign::geometry::LineString lsBoundaryCutByAngle = vLsBorderCutByAngle[i];

			ign::geometry::algorithm::BufferOpGeos buffOp;
			ign::geometry::GeometryPtr buffBorder(buffOp.buffer(lsBoundaryCutByAngle, distBuffer, 0, ign::geometry::algorithm::BufferOpGeos::CAP_FLAT));

			ign::feature::Feature fBuff;
			fBuff.setGeometry(buffBorder->clone());

			_getCLfromBorder(lsBoundaryCutByAngle, buffBorder, distBuffer, thresholdNoCL, angleMaxBorder, ratioInBuff, snapOnVertexBorder);
		}
	}

	_logger->log(epg::log::TITLE, "[ END CL GENERATION BY COUNTRY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::MergeConnectingLinesOnBorder(
	std::string countryCodeDouble,
	bool verbose
) {
	CFeatGenerationOp op(countryCodeDouble, verbose);
    op._mergeConnectingLinesOnBorder();
}

///
///
///
void app::calcul::CFeatGenerationOp::_mergeConnectingLinesOnBorder() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN CL MERGING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distMaxEdges = themeParameters->getValue(CL_EDGE_MAX_DIST).toDouble();
	double const snapProjCl2edge = themeParameters->getValue(CL_SNAP_PROJ_CL_2_EDGE_DIST).toDouble();

	_mergeIntersectingClWithGraph(distMaxEdges, snapProjCl2edge);

	_logger->log(epg::log::TITLE, "[ END CL MERGING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::SnapConnectingLines(std::string countryCodeDouble, bool verbose)
{
	CFeatGenerationOp op(countryCodeDouble, verbose);
    op._snapConnectingLines();
}

///
///
///
void app::calcul::CFeatGenerationOp::_snapConnectingLines() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN CL SNAPPING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distMaxClClosest = themeParameters->getValue(CL_CL_CLOSEST_MAX_DIST).toDouble();

	_snapCl2Cl( distMaxClClosest);

	_logger->log(epg::log::TITLE, "[ END CL SNAPPING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::DeleteConnectingLines(std::string countryCodeDouble, bool verbose)
{
	CFeatGenerationOp op(countryCodeDouble, verbose);
    op._deleteConnectingLines();
}

///
///
///
void app::calcul::CFeatGenerationOp::_deleteConnectingLines() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN CL DELETING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double angleMaxEdges = themeParameters->getValue(CL_EDGE_MAX_ANGLE).toDouble();
	double const distMaxEdges = themeParameters->getValue(CL_EDGE_MAX_DIST).toDouble();
	double const snapProjCl2edge = themeParameters->getValue(CL_SNAP_PROJ_CL_2_EDGE_DIST).toDouble();

	angleMaxEdges = angleMaxEdges * M_PI / 180;

	_deleteClByAngleAndDistEdges(angleMaxEdges, distMaxEdges, snapProjCl2edge);

	_deleteCLUnderThreshold();

	_logger->log(epg::log::TITLE, "[ END CL DELETING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::UpdateGeomConnectingLines(
	std::string countryCodeDouble,
	bool verbose
) {
	CFeatGenerationOp op(countryCodeDouble, verbose);
    op._updateGeomConnectingLines();
}

///
///
///
void app::calcul::CFeatGenerationOp::_updateGeomConnectingLines() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN LOAD GRAPH CL FOR CONTINUITY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	GraphType graphCL;
	_loadGraphCL(graphCL);

	_logger->log(epg::log::TITLE, "[ END LOAD GRAPH CL FOR CONTINUITY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());


	_logger->log(epg::log::TITLE, "[ BEGIN CL UPDATE GEOMETRY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const snapProjCl2edge = themeParameters->getValue(CL_SNAP_PROJ_CL_2_EDGE_DIST).toDouble();
	epg::Context* context = epg::ContextS::getInstance();
	std::string idName = context->getEpgParameters().getValue(ID).toString();
	std::string geomName = context->getEpgParameters().getValue(GEOM).toString();

	_updateGeomCL( snapProjCl2edge);

	_logger->log(epg::log::TITLE, "[ END CL UPDATE GEOMETRY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	//-- copie intermediaire pour debug
	epg::utils::CopyTableUtils::copyTable(
		themeParameters->getValue(CL_TABLE).toString(),
		idName,
		geomName,
		ign::geometry::Geometry::GeometryTypeLineString,
		themeParameters->getValue(CL_TABLE).toString()+"_no_cont",
		"", false, true
	);

	_logger->log(epg::log::TITLE, "[ BEGIN CL UPDATE CONTINUITY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	_setContinuityCl(graphCL);

	_logger->log(epg::log::TITLE, "[ END CL UPDATE CONTINUITY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}



///
///
///
void app::calcul::CFeatGenerationOp::ComputeCL(
	std::string countryCodeDouble,
	bool verbose
) {
	CFeatGenerationOp op(countryCodeDouble, verbose);
    op._computeCL();
}

///
///
///
void app::calcul::CFeatGenerationOp::_computeCL() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN CL GENERATION FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	_generateConnectingLinesByCountry();
	_mergeConnectingLinesOnBorder();
    _snapConnectingLines();
	_deleteConnectingLines();
	_updateGeomConnectingLines();
	_logger->log(epg::log::TITLE, "[ END CL GENERATION FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

}

///
///
///
void app::calcul::CFeatGenerationOp::ComputeCP(
	std::string countryCodeDouble,
	bool verbose
) {
	CFeatGenerationOp op(countryCodeDouble, verbose);
    op._computeCP();
}

///
///
///
void app::calcul::CFeatGenerationOp::_computeCP() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN CP GENERATION FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	epg::Context* context = epg::ContextS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distBuffer = themeParameters->getValue(CP_BUFFER_DIST).toDouble();
	double const distCLIntersected = themeParameters->getValue(CP_INTERSECTED_CL_DIST).toDouble();
	double const distUnderShoot = themeParameters->getValue(CP_UNDERSHOOT_DIST).toDouble();
	double const distCp2snapCl = themeParameters->getValue(CP_CP_2_CL_SNAP_DIST).toDouble();
	double const snapDistOnVertexFromCl = themeParameters->getValue(CP_VERTEX_CL_SNAP_DIST).toDouble();

	ign::feature::FeatureIteratorPtr itBoundary = _fsBoundary->getFeatures(ign::feature::FeatureFilter(countryCodeName + " = '" + _countryCodeDouble + "'"));
	while (itBoundary->hasNext()) {
		ign::feature::Feature fBoundary = itBoundary->next();
		_logger->log(epg::log::INFO, "id boundary :" + fBoundary.getId());
		std::string boundaryType = fBoundary.getAttribute("boundary_type").toString();
		if (boundaryType != "international_boundary" && boundaryType.find("coastline_sea_limit") == -1)
			continue;
		ign::geometry::LineString lsBoundary = fBoundary.getGeometry().asLineString();
		ign::geometry::algorithm::BufferOpGeos buffOp;
		ign::geometry::GeometryPtr buffBorder(buffOp.buffer(lsBoundary, distUnderShoot, 0, ign::geometry::algorithm::BufferOpGeos::CAP_FLAT));

		_getCPfromIntersectBorder(lsBoundary, distCLIntersected);

		_addToUndershootNearBorder(lsBoundary, buffBorder, distUnderShoot);
	}

	_snapCPNearBy(0);

	std::map<std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint>> mClSplittedByCp;
	_snapCpOnClNearBy(distCp2snapCl, snapDistOnVertexFromCl, mClSplittedByCp);
	_cutClByCp(mClSplittedByCp);

	_logger->log(epg::log::TITLE, "[ END CP GENERATION FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}
	
///
///
///
void app::calcul::CFeatGenerationOp::_init(
	std::string countryCodeDouble,
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
	std::string const cpTableName = themeParameters->getValue(CP_TABLE).toString();
	std::string const clTableName = themeParameters->getValue(CL_TABLE).toString();
	
	_shapeLogger = epg::log::ShapeLoggerS::getInstance();
	_shapeLogger->addShape("CLBeforeMerge", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClMergedBeforeUpdate", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDeletedNoCandidatefound", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDoublon", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDeleteByAngleDistEdges", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDebug", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("edgeClCutByCp", epg::log::ShapeLogger::LINESTRING); 
	_shapeLogger->addShape("lsBorderCutByAngle", epg::log::ShapeLogger::LINESTRING); 

	_countryCodeDouble = countryCodeDouble;
	epg::tools::StringTools::Split(_countryCodeDouble, "#", _vCountriesCodeName);
	_verbose = verbose;

	///recuperation de la liste des attributs à concatener, de w et json dans la fusion des attributs
	std::string listAttrWName = themeParameters->getValue(LIST_ATTR_W).toString();
	std::string listAttrJsonName = themeParameters->getValue(LIST_ATTR_JSON).toString();
	_attrMergerOnBorder.setLists( listAttrWName, listAttrJsonName, "/");
	

	std::string listValueFormwayBigDist2merge = themeParameters->getValue(CP_VALUE_FORMWAY_BIGDIST2MERGE).toString();
	std::vector<std::string> vValueFormwayBigDist2merge;
	epg::tools::StringTools::Split(listValueFormwayBigDist2merge, "/", vValueFormwayBigDist2merge);
	for (size_t i = 0; i < vValueFormwayBigDist2merge.size(); ++i) {
		_sFormwayValues4BigDist2Merge.insert(vValueFormwayBigDist2merge[i]);
	}

	///recuperation des features
	std::string const boundaryTableName = epg::utils::replaceTableName(context->getEpgParameters().getValue(TARGET_BOUNDARY_TABLE).toString());
	_fsBoundary = context->getDataBaseManager().getFeatureStore(boundaryTableName, idName, geomName);

	std::string const landmaskTableName = epg::utils::replaceTableName(themeParameters->getValue(LANDMASK_TABLE).toString());
	_fsLandmask= context->getDataBaseManager().getFeatureStore(landmaskTableName, idName, geomName);

	_fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);

	_reqFilterEdges2generateCF = themeParameters->getValue(SQL_FILTER_EDGES_2_GENERATE_CF).toString();

	//--
	_fsCP = context->getDataBaseManager().getFeatureStore(cpTableName, idName, geomName);
	_fsCL = context->getDataBaseManager().getFeatureStore(clTableName, idName, geomName);

	// id generator
	_idGeneratorCP = epg::sql::tools::IdGeneratorInterfacePtr(epg::sql::tools::IdGeneratorFactory::getNew(*_fsCP, "CONNECTINGPOINT"));
	_idGeneratorCL = epg::sql::tools::IdGeneratorInterfacePtr(epg::sql::tools::IdGeneratorFactory::getNew(*_fsCL, "CONNECTINGLINE"));



	std::string  boundarySmoothedTableName;
	if (themeParameters->getValue(BOUNDARY_SMOOTHED_TABLE).toString() != "")
		boundarySmoothedTableName = epg::utils::replaceTableName(themeParameters->getValue(BOUNDARY_SMOOTHED_TABLE).toString());
	else
		boundarySmoothedTableName = boundaryTableName;
	ign::feature::FeatureFilter filter(countryCodeName + " = '" + _countryCodeDouble + "'");
	_mlsBorderSmoothed = new epg::tools::MultiLineStringTool(filter, *context->getDataBaseManager().getFeatureStore(boundarySmoothedTableName, idName, geomName));
	

	_logger->log(epg::log::TITLE, "[ END INITIALIZATION ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::_getBorderCutByAngle(
	ign::geometry::LineString & lsBorder,
	std::vector<ign::geometry::LineString> & vLsBorderCutByAngle,
	double angleMaxToCutBorder
) const {
	ign::geometry::LineString lsCurr(lsBorder.pointN(0), lsBorder.pointN(1));

	for (size_t i = 1; i < lsBorder.numPoints(); ++i) {

		if (i == lsBorder.numPoints() - 1) {
			vLsBorderCutByAngle.push_back(lsCurr);
			break;
		}
		ign::math::Vec2d vecLsCurr(lsCurr.endPoint().x() - lsCurr.startPoint().x(), lsCurr.endPoint().y() - lsCurr.startPoint().y());
		
		ign::geometry::LineString lsNext(lsBorder.pointN(i), lsBorder.pointN(i + 1));
		ign::math::Vec2d vecLsNext(lsNext.endPoint().x() - lsNext.startPoint().x(), lsNext.endPoint().y() - lsNext.startPoint().y());

		double angleBorderWithNextPortion = epg::tools::geometry::angle(vecLsCurr, vecLsNext);

		if (angleBorderWithNextPortion > angleMaxToCutBorder) {
			vLsBorderCutByAngle.push_back(lsCurr);
			lsCurr = lsNext;
		}
		else
			lsCurr.addPoint(lsBorder.pointN(i + 1));
	}

	for (size_t i = 0; i < vLsBorderCutByAngle.size(); ++i) {
		ign::feature::Feature fShaplog;
		fShaplog.setGeometry(vLsBorderCutByAngle[i]);
		_shapeLogger->writeFeature("lsBorderCutByAngle", fShaplog);
	}
}

///
///
///
void app::calcul::CFeatGenerationOp::_getCLfromBorder(
	ign::geometry::LineString & lsBorder,
	ign::geometry::GeometryPtr & buffBorder,
	double distBuffer,
	double thresholdNoCL,
	double angleMax,
	double ratioInBuff,
	double snapOnVertexBorder
) const {
	epg::Context* context = epg::ContextS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	std::vector<ign::feature::FeatureAttributeType> listAttrEdge = _fsEdge->getFeatureType().attributes();

	ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + buffBorder->toString() + "'),3035))");
	epg::tools::FilterTools::addAndConditions(filter, countryCodeName +" NOT LIKE '%#%'");
	//ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_GeomFromText('" + buffBorder->toString() + "'))");
	if (_reqFilterEdges2generateCF != "")
		epg::tools::FilterTools::addAndConditions(filter, _reqFilterEdges2generateCF);

	ign::feature::FeatureIteratorPtr eit = _fsEdge->getFeatures(filter);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsEdge, filter);
	if (numFeatures == 0)
		return;
	boost::progress_display display(numFeatures, std::cout, "[ CREATE CONNECTING LINES ]\n");

	while (eit->hasNext())
	{
		++display;
		//
		ign::feature::Feature fEdge = eit->next();
		ign::geometry::LineString lsEdge = fEdge.getGeometry().asLineString();
	
		//ign::geometry::Geometry* geomIntersect = lsEdge.Intersection(*buff);
		std::vector<ign::geometry::LineString> vLsProjectedOnBorder;
		app::geometry::tools::LineStringSplitter lsSplitter(lsEdge);
		lsSplitter.addCuttingGeometry(*buffBorder);
		std::vector<ign::geometry::LineString> subEdgesBorder = lsSplitter.getSubLineStringsZ();

		//pas d'intersection par le buffer
		if (subEdgesBorder.size() == 1) {			
			double angleEdgBorder = _getAngleEdgeWithBorder(lsEdge,lsBorder);		
			//si l'edge est "proche" on considere qu'il est entierement dans le buffer et longe la frontiere
			if (lsEdge.distance(lsBorder) < distBuffer && (angleEdgBorder < angleMax || angleEdgBorder > (M_PI - angleMax) ) ) {
				ign::geometry::LineString lsCL;
				_getGeomCL(lsCL, *_mlsBorderSmoothed, lsEdge, distBuffer,snapOnVertexBorder);
				if (lsCL.numPoints() >= 2) {
					vLsProjectedOnBorder.push_back(lsCL);
				}
			}
		}
		else 
		{
			int numfirstSubInBuff = -1;
			int numlastSubInBuff = -1;
			int lengthInBuff = 0;
			int lengthNearByBuff = 0;

			for (size_t i = 0; i < subEdgesBorder.size(); ++i) {
				ign::geometry::LineString lsSubEdgeCurr = subEdgesBorder[i];

 				double angleSubEdgBorder = _getAngleEdgeWithBorder(lsSubEdgeCurr, lsBorder);

				int numSeg = static_cast<int>(std::floor(lsSubEdgeCurr.numSegments() / 2.));
				ign::geometry::Point interiorPointSEC = epg::tools::geometry::interpolate(lsSubEdgeCurr, numSeg, 0.5);
				bool isSubSegInBuff = false;

				if (buffBorder->contains(interiorPointSEC) && (angleSubEdgBorder < angleMax || angleSubEdgBorder > (M_PI - angleMax) ) ) {
					isSubSegInBuff = true;
					numlastSubInBuff = i;

					lengthInBuff += lsSubEdgeCurr.length();
					if (numfirstSubInBuff < 0)
						numfirstSubInBuff = i;
				}

				if (isSubSegInBuff || lsSubEdgeCurr.length() <= thresholdNoCL)
					lengthNearByBuff += lsSubEdgeCurr.length();

				if ((lsSubEdgeCurr.length() > thresholdNoCL && !isSubSegInBuff) || i == subEdgesBorder.size() - 1) {
					if (lengthInBuff > ratioInBuff * lengthNearByBuff) {
						//recup ptStart, ptFin et proj des pt sur la border
						//recup de la border entre ces points pour recup de la geom CL
						ign::geometry::LineString lsCL;
						_getGeomCL(lsCL, *_mlsBorderSmoothed, subEdgesBorder[numfirstSubInBuff], distBuffer, snapOnVertexBorder);
						if (lsCL.numPoints() >= 2 ) {
							vLsProjectedOnBorder.push_back(lsCL);
						}
					
					}

					//reset
					numfirstSubInBuff = -1;
					numlastSubInBuff = -1;
					lengthInBuff = 0;
					lengthNearByBuff = 0;
				}
			}
		}

		for (size_t i = 0; i < vLsProjectedOnBorder.size(); ++i) {
			//create CL
			//generation de l'id CL
			//creation du feat en copie du featEdge puis modif de la geom et de l'id
			ign::feature::Feature featCL = _fsCL->newFeature();
			vLsProjectedOnBorder[i].setFillZ(0);
			featCL.setGeometry(vLsProjectedOnBorder[i]);
			std::string idCL = _idGeneratorCL->next();
			featCL.setId(idCL);
			featCL.setAttribute(linkedFeatIdName, ign::data::String(fEdge.getId()));
			for (std::vector<ign::feature::FeatureAttributeType>::iterator vit = listAttrEdge.begin();
				vit != listAttrEdge.end(); ++vit) {
				std::string attrName = vit->getName();
				if (attrName == geomName || attrName == idName || !_fsCL->getFeatureType().hasAttribute(attrName) )
					continue;
				featCL.setAttribute(attrName,fEdge.getAttribute(attrName));
			}
			_fsCL->createFeature(featCL, idCL);
			//shape CL
			{
				ign::feature::Feature fShaplog = featCL;
				ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
				lsSphaplog.clearZ();
				fShaplog.setGeometry(lsSphaplog);
				_shapeLogger->writeFeature("CLBeforeMerge", fShaplog);
			}
		}

	}	
}

///
///
///
void app::calcul::CFeatGenerationOp::_addToUndershootNearBorder(
	ign::geometry::LineString & lsBorder,
	ign::geometry::GeometryPtr & buffBorder,
	double distUnderShoot
) const {
	epg::Context* context = epg::ContextS::getInstance();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	std::vector<ign::feature::FeatureAttributeType> listAttrEdge = _fsEdge->getFeatureType().attributes();

	//ign::feature::FeatureFilter filterBuffBorder("ST_INTERSECTS(" + geomName + ", ST_GeomFromText('" + buffBorder->toString() + "'))");
	ign::feature::FeatureFilter filterBuffBorder("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + buffBorder->toString() + "'),3035))");
	if (_reqFilterEdges2generateCF != "")
		epg::tools::FilterTools::addAndConditions(filterBuffBorder, _reqFilterEdges2generateCF);
	//on ne prend que les edes ayant un cc simple pour ne pas créer de CP là où il y a des CLs
	epg::tools::FilterTools::addAndConditions(filterBuffBorder, "(" + countryCodeName + " = '" + _vCountriesCodeName[0] + "' or " + countryCodeName + " = '" + _vCountriesCodeName[1] + "')");

	ign::feature::FeatureIteratorPtr eit = _fsEdge->getFeatures(filterBuffBorder);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsEdge, filterBuffBorder);
	boost::progress_display display(numFeatures, std::cout, "[ GET CONNECTING POINTS BY UNDERSHOOT LINES ]\n");

	epg::tools::geometry::SegmentIndexedGeometry segIndexLsBorder(&lsBorder);

	//recuperation des troncons qui intersect le buff de 5m
	while (eit->hasNext())
	{
		++display;
		ign::feature::Feature fEdge = eit->next();
		ign::geometry::LineString lsEdge = fEdge.getGeometry().asLineString();
		//si intersection border -> on fait rien
		if (lsBorder.intersects(lsEdge))
			//if (!lsBorder.Intersection(lsEdge)->isNull())
			//est ce que intersects plus long que intersection?
			continue;

		double distBorder2StartPt = lsBorder.distance(lsEdge.startPoint());
		double distBorder2EndPt = lsBorder.distance(lsEdge.endPoint());
		ign::geometry::Point ptClosestBorder;
		ign::math::Vec2d vecToBorder;
		if (distBorder2StartPt < distBorder2EndPt) {
			ptClosestBorder = lsEdge.startPoint();
		}
		else {
			ptClosestBorder = lsEdge.endPoint();
			lsEdge.reverse();
		}
		vecToBorder.x() = lsEdge.pointN(1).x() - ptClosestBorder.x();
		vecToBorder.y() = lsEdge.pointN(1).y() - ptClosestBorder.y();
		
		//on verifie que le point est un dangle, sinon on fait rien
		ign::feature::FeatureFilter filterArroundPt;
		if (_reqFilterEdges2generateCF != "")
			filterArroundPt.setPropertyConditions(_reqFilterEdges2generateCF);
		filterArroundPt.setExtent(ptClosestBorder.getEnvelope().expandBy(1));
		ign::feature::FeatureIteratorPtr eitArroundPt = _fsEdge->getFeatures(filterArroundPt);
		bool isPtADangle = true;
		while (eitArroundPt->hasNext()) {
			ign::feature::Feature featArroundPt = eitArroundPt->next();
			if (featArroundPt.getId() == fEdge.getId())
				continue;
			double dist = featArroundPt.getGeometry().distance(ptClosestBorder);
			if (dist > 0)
				continue;
			isPtADangle = false;
			break;
		}
		if (!isPtADangle)
			continue;

		ign::geometry::Point projPt;
		std::vector<ign::geometry::LineString> vBorderSegments;
		segIndexLsBorder.getSegments( ptClosestBorder.getEnvelope().expandBy(distUnderShoot*1.5), vBorderSegments );
		if (vBorderSegments.empty())
			continue;
		ign::geometry::MultiLineString mlsBorderSegments(vBorderSegments);
		////////////////////todo recup le bon segment de la frontiere
		std::vector< ign::geometry::Point > vPtIntersect = epg::tools::geometry::LineIntersector::compute(ptClosestBorder, lsEdge.pointN(1), mlsBorderSegments);
		double distMin = 100000;
		for (std::vector< ign::geometry::Point >::iterator vit = vPtIntersect.begin(); vit != vPtIntersect.end(); ++vit) {
			double dist = ptClosestBorder.distance(*vit);
			if (dist < distMin) {
				projPt = *vit;
				distMin = dist;
			}
		}

		// On privilégie une projection classique si le projeté est (beaucoup) plus près
		ign::geometry::Point projPt2 = epg::tools::geometry::project(mlsBorderSegments, ptClosestBorder, 0);
		double distance2 = ptClosestBorder.distance(projPt2);
		if (distance2 < distMin/3) {
			distMin = distance2;
			projPt = projPt2;
		}

		if (distMin > distUnderShoot )
				continue;

		projPt.setZ(0);
		ign::feature::Feature fCF = _fsCP->newFeature();
		fCF.setGeometry(projPt);
		std::string idCP = _idGeneratorCP->next();
		fCF.setId(idCP);
		fCF.setAttribute(linkedFeatIdName, ign::data::String(fEdge.getId()));
		for (std::vector<ign::feature::FeatureAttributeType>::iterator vit = listAttrEdge.begin();
			vit != listAttrEdge.end(); ++vit) {
			std::string attrName = vit->getName();
			if (attrName == geomName || attrName == idName || !_fsCP->getFeatureType().hasAttribute(attrName))
				continue;
			fCF.setAttribute(attrName, fEdge.getAttribute(attrName));
		}
		//fCF.setAttribute("w_calcul", ign::data::String("1"));
		_fsCP->createFeature(fCF, idCP);
		{
			ign::feature::Feature fShaplog = fCF;
			ign::geometry::Point ptProjNoZ = projPt;
			ptProjNoZ.clearZ();
			fShaplog.setGeometry(ptProjNoZ);
			_shapeLogger->writeFeature("CPUndershoot", fShaplog);
		}
	}
}

///
///
///
void app::calcul::CFeatGenerationOp::_getCPfromIntersectBorder(
	ign::geometry::LineString & lsBorder,
	double distCLIntersected
) const {
	epg::Context* context = epg::ContextS::getInstance();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	std::vector<ign::feature::FeatureAttributeType> listAttrEdge = _fsEdge->getFeatureType().attributes();

	//on prend les edges qui intersectent les frontières et ou ceux qui intersectent les Cls
	ign::feature::FeatureFilter filterFeaturesToMatch("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + lsBorder.toString() + "'),3035))"
		" OR ST_INTERSECTS(" + geomName + ", (SELECT ST_Union(array(SELECT "+geomName+" FROM "+ _fsEdge->getTableName()+" WHERE "+ countryCodeName + " = '" + _countryCodeDouble + "'))))");

	if (_reqFilterEdges2generateCF != "")
		epg::tools::FilterTools::addAndConditions(filterFeaturesToMatch, _reqFilterEdges2generateCF);
	//on ne prend que les edes ayant un cc simple pour ne pas créer de CP là où il y a des CLs
	epg::tools::FilterTools::addAndConditions(filterFeaturesToMatch,"("+ countryCodeName+" = '" +_vCountriesCodeName[0]+"' or "+countryCodeName + " = '" + _vCountriesCodeName[1] +"')");
	ign::feature::FeatureIteratorPtr itFeaturesToMatch = _fsEdge->getFeatures(filterFeaturesToMatch);

	int numFeatures = context->getDataBaseManager().numFeatures(*_fsEdge, filterFeaturesToMatch);
	boost::progress_display display(numFeatures, std::cout, "[ CREATE CONNECTING POINTS ]\n");

	while (itFeaturesToMatch->hasNext())
	{
		++display;

		ign::feature::Feature fToMatch = itFeaturesToMatch->next();
		ign::geometry::LineString lsFToMatch = fToMatch.getGeometry().asLineString();
		ign::geometry::GeometryPtr geomPtr(lsFToMatch.Intersection(lsBorder));

		ign::feature::Feature fClArround;
		bool hasClConnected = _isEdgeConnected2cl(lsFToMatch, lsFToMatch.getEnvelope().expandBy(distCLIntersected), fClArround, distCLIntersected);
		
		if (hasClConnected) {
			ign::geometry::GeometryPtr geomIntersectCl(fClArround.getGeometry().Intersection(lsFToMatch));

			//si il existe une intersection entre l'edge et une CL, et qu'elle est sous un seuil, on prend cette intersection à la place de celle avec la frontière
			if (!geomIntersectCl->isNull() && !geomIntersectCl->isEmpty()) {
				double distance = geomIntersectCl->distance(lsBorder);

				if(distance >= 0 && distance < distCLIntersected ) {
					if (geomPtr->isNull() || geomPtr->distance(*geomIntersectCl) < distCLIntersected) //utile ce test ou le precedent suffit?
						geomPtr.reset(geomIntersectCl.release());
					}
			}
		}

		if (geomPtr->isNull() || geomPtr->isEmpty())
			continue;

		ign::feature::Feature fCF = _fsCP->newFeature();
		fCF.setAttribute(linkedFeatIdName, ign::data::String(fToMatch.getId()));
		for (std::vector<ign::feature::FeatureAttributeType>::iterator vit = listAttrEdge.begin();
			vit != listAttrEdge.end(); ++vit) {
			std::string attrName = vit->getName();
			if (attrName == geomName || attrName == idName || !_fsCP->getFeatureType().hasAttribute(attrName))
				continue;
			fCF.setAttribute(attrName, fToMatch.getAttribute(attrName));
		}

		if (geomPtr->isPoint())
		{
			fCF.setGeometry(geomPtr->asPoint());
			std::string idCP = _idGeneratorCP->next();
			_fsCP->createFeature(fCF, idCP);
		}

		else if (geomPtr->isGeometryCollection())
		{
			ign::geometry::GeometryCollection geomCollect = geomPtr->asGeometryCollection();
			for (size_t i = 0; i < geomCollect.numGeometries(); ++i)
			{
				if (geomCollect.geometryN(i).isPoint())
				{
					ign::geometry::Point ptIntersect = geomCollect.geometryN(i).asPoint();		
					fCF.setGeometry(ptIntersect);
					std::string idCP = _idGeneratorCP->next();
					_fsCP->createFeature(fCF, idCP);
				}
			}
		}
	}

}

///
///
///
double app::calcul::CFeatGenerationOp::_getAngleEdgeWithBorder(
	ign::geometry::LineString & lsEdge,
	ign::geometry::LineString & lsBorder
) const {
	double angleEdgWBorder;
	ign::math::Vec2d vecEdge(lsEdge.endPoint().x() - lsEdge.startPoint().x(), lsEdge.endPoint().y() - lsEdge.startPoint().y());

	ign::geometry::Point ptStartProjOnBorder= epg::tools::geometry::project(lsBorder, lsEdge.startPoint(), 0);
	ign::geometry::Point ptEndProjOnBorder = epg::tools::geometry::project(lsBorder, lsEdge.endPoint(), 0);
	
	ign::math::Vec2d vecBorder(ptEndProjOnBorder.x() - ptStartProjOnBorder.x(), ptEndProjOnBorder.y() - ptStartProjOnBorder.y());

	angleEdgWBorder = epg::tools::geometry::angle(vecBorder, vecEdge);

	return angleEdgWBorder;
	//	//double angleSubEdgBorder = epg::tools::geometry::angle(vecBorder, vecSubEdge);
}

///
///
///
void app::calcul::CFeatGenerationOp::_getGeomCL(
	ign::geometry::LineString & lsCL,
	epg::tools::MultiLineStringTool & mslBorder,
	ign::geometry::LineString & lsStart2EndToPrject,
	double distMaxBorder,
	double snapOnVertexBorder
) const {

	ign::geometry::Point ptStartToProject = lsStart2EndToPrject.startPoint();
	ign::geometry::Point ptEndToProject = lsStart2EndToPrject.endPoint();
	// std::pair< bool, ign::geometry::LineString > pathFound = mslBorder.getPathAlong(ptStartToProject, ptEndToProject, lsStart2EndToPrject, 2* distMaxBorder, distMaxBorder+1);
	std::pair< bool, ign::geometry::LineString > pathFound = mslBorder.getPathAlong(ptStartToProject, ptEndToProject, lsStart2EndToPrject, 1000, 1000, snapOnVertexBorder);

	if (pathFound.first)
		lsCL = pathFound.second;

}

///
///
///
void app::calcul::CFeatGenerationOp::_snapCPNearBy(
	double snapOnVertexBorder
) const {

	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distMergeCP = themeParameters->getValue(CP_MERGE_DIST_CP).toDouble();
	double const distMergeTractorCP = themeParameters->getValue(CP_MERGE_DIST_TRACTOR_CP).toDouble();
	double maxDistMerge = std::max(distMergeCP, distMergeTractorCP);

	ign::feature::FeatureFilter filterCP;
	
	for (size_t i = 0; i < _vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCP, countryCodeName + " = '" + _vCountriesCodeName[i] + "'");
	}
	ign::feature::FeatureIteratorPtr itCP = _fsCP->getFeatures(filterCP);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCP, filterCP);
	// boost::progress_display display(numFeatures, std::cout, "[ FUSION CONNECTING POINTS WITH #]\n");

	std::set<std::string> sCP2Snap;
	std::string separator = "#";

	while (itCP->hasNext())
	{
		// ++display;

		ign::feature::Feature fCPCurr = itCP->next();
		std::string idCP = fCPCurr.getId();

		if (sCP2Snap.find(idCP) != sCP2Snap.end())
			continue;

		std::map<std::string, ign::feature::Feature> mCPNear;
		bool hasNearestCP = _getNearestCP(fCPCurr, maxDistMerge, mCPNear);

		std::list<std::string> lCp2Delete;

		if (!hasNearestCP) {
			_fsCP->deleteFeature(fCPCurr.getId());
			continue;
		}

		std::set<std::string> s1, s2;
		std::map<std::string, ign::feature::Feature> mLinkedEdgeFeature;
		for(std::map<std::string, ign::feature::Feature>::iterator mit = mCPNear.begin(); mit != mCPNear.end(); ++mit) {
			if(mit->second.getAttribute(countryCodeName).toString() == _vCountriesCodeName[0]) s1.insert(mit->first);
			else s2.insert(mit->first);
			sCP2Snap.insert(mit->first);

			mLinkedEdgeFeature.insert(std::make_pair(mit->first, ign::feature::Feature()));

			_fsEdge->getFeatureById(mit->second.getAttribute(linkedFeatIdName).toString(), mLinkedEdgeFeature[mit->first]);
		}

		// map pour optimisation
		std::map<std::string, std::map<std::string, std::pair<bool, double>>> mmAreMergeable;
		for (std::set<std::string>::const_iterator sit1 = s1.begin() ; sit1 != s1.end() ; ++sit1) {
			std::map<std::string, std::pair<bool, double>> mAreMergeable;
			for (std::set<std::string>::const_iterator sit2 = s2.begin() ; sit2 != s2.end() ; ++sit2) {
				double distance = mCPNear[*sit1].getGeometry().asPoint().distance(mCPNear[*sit2].getGeometry().asPoint());
				bool areMergeable = _areMergeable(mLinkedEdgeFeature[*sit1], mLinkedEdgeFeature[*sit2], distance);
				mAreMergeable.insert(std::make_pair(*sit2, std::make_pair(areMergeable, distance)));
			}
			mmAreMergeable.insert(std::make_pair(*sit1, mAreMergeable));
		}

		std::map<std::string, std::string> m1;
		for (std::set<std::string>::const_iterator sit1 = s1.begin() ; sit1 != s1.end() ; ++sit1) {
			std::string samicopain = "";
			double distanceMax = std::numeric_limits<double>::max();
			for (std::set<std::string>::const_iterator sit2 = s2.begin() ; sit2 != s2.end() ; ++sit2) {
				double distance = mmAreMergeable[*sit1][*sit2].second;
				if (distance < distanceMax && mmAreMergeable[*sit1][*sit2].first) {
					samicopain = *sit2;
					distanceMax = distance;
				}
			}
			if (!samicopain.empty()) {
				m1.insert(std::make_pair(*sit1, samicopain));
			} else {
				lCp2Delete.push_back(*sit1);
			}
		}

		std::map<std::string, std::string> m2;
		for (std::set<std::string>::const_iterator sit2 = s2.begin() ; sit2 != s2.end() ; ++sit2) {
			std::string samicopain = "";
			double distanceMax = std::numeric_limits<double>::max();
			for (std::set<std::string>::const_iterator sit1 = s1.begin() ; sit1 != s1.end() ; ++sit1) {
				double distance = mmAreMergeable[*sit1][*sit2].second;;
				if (distance < distanceMax && mmAreMergeable[*sit1][*sit2].first) {
					samicopain = *sit1;
					distanceMax = distance;
				}
			}
			if (!samicopain.empty()) {
				m2.insert(std::make_pair(*sit2, samicopain));
			} else {
				lCp2Delete.push_back(*sit2);
			}
		}

		typedef boost::bimap<boost::bimaps::set_of<std::string>, boost::bimaps::multiset_of<size_t>> bimap_t;
		typedef bimap_t::value_type value_type;

		bimap_t mapCpGroup;
		size_t group = 0;
		size_t nMapCpGroup = mapCpGroup.size();
		for (std::map<std::string, std::string>::const_iterator mit1 = m1.begin() , next_mit1 = mit1 ; mit1 != m1.end() ; mit1 = next_mit1) {
			++next_mit1;
			for (std::map<std::string, std::string>::const_iterator mit2 = m2.begin() ; mit2 != m2.end() ; ++mit2) {
				if( mit1->second == mit2->first && mit2->second == mit1->first) {
					mapCpGroup.insert(value_type(mit1->first, ++group));
					mapCpGroup.insert(value_type(mit2->first, group));
					m1.erase(mit1);
					m2.erase(mit2);
					break;
				}
			}
		}

		for (std::map<std::string, std::string>::const_iterator mit1 = m1.begin() , next_mit1 = mit1 ; mit1 != m1.end() ; mit1 = next_mit1) {
			++next_mit1;
			std::map<std::string, std::string>::const_iterator mit2 = m2.find(mit1->second);
			if( mit2 == m2.end()) continue; /*cp de m2 déjà affecté à un groupe*/

			mapCpGroup.insert(value_type(mit1->first,++group));
			mapCpGroup.insert(value_type(mit2->first,group));
			m1.erase(mit1);
			m2.erase(mit2);
		}

		for (std::map<std::string, std::string>::const_iterator mit2 = m2.begin() , next_mit2 = mit2 ; mit2 != m2.end() ; mit2 = next_mit2) {
			++next_mit2;
			std::map<std::string, std::string>::const_iterator mit1 = m1.find(mit2->second);
			if( mit1 == m1.end()) continue; /*cp de m1 déjà affecté à un groupe*/

			mapCpGroup.insert(value_type(mit2->first,++group));
			mapCpGroup.insert(value_type(mit1->first,group));
			m2.erase(mit2);
			m1.erase(mit1);
		}

		// il reste les cl reliés à des groupes
		for (std::map<std::string, std::string>::const_iterator mit2 = m2.begin() ; mit2 != m2.end() ; ++mit2) {
			auto l_mit = mapCpGroup.left.find(mit2->second);
			IGN_ASSERT(l_mit != mapCpGroup.left.end());
			mapCpGroup.insert(value_type(mit2->first, l_mit->second));
		}
		for (std::map<std::string, std::string>::const_iterator mit1 = m1.begin() ; mit1 != m1.end() ; ++mit1) {
			auto l_mit = mapCpGroup.left.find(mit1->second);
			IGN_ASSERT(l_mit != mapCpGroup.left.end());
			mapCpGroup.insert(value_type(mit1->first, l_mit->second));
		}

		for (size_t i = 1 ; i <= group ; ++i) {
			ign::geometry::MultiPoint multiPtCP;
			auto range = mapCpGroup.right.equal_range(i);
			for (auto r_it = range.first; r_it != range.second; ++r_it) {
				multiPtCP.addGeometry(mCPNear[r_it->second].getGeometry().asPoint());
			}

			//geom
			ign::geometry::Point ptCentroidCP = multiPtCP.asMultiPoint().getCentroid();
			ign::feature::FeatureFilter filterBorderNearCP;// (countryCodeName + " = 'be#fr'");
			filterBorderNearCP.setExtent(ptCentroidCP.getEnvelope().expandBy(maxDistMerge));
			ign::geometry::LineString lsBorderClosest;
			double distMinBorder = 2 * maxDistMerge;
			ign::feature::FeatureIteratorPtr fitBorder = _fsBoundary->getFeatures(filterBorderNearCP);
			while (fitBorder->hasNext()) {
				ign::geometry::LineString lsBorder = fitBorder->next().getGeometry().asLineString();
				double dist = lsBorder.distance(ptCentroidCP);
				if (dist < distMinBorder) {
					distMinBorder = dist;
					lsBorderClosest = lsBorder;
				}
			}
			ign::geometry::Point ptCentroidOnBorderCP = epg::tools::geometry::project(lsBorderClosest, ptCentroidCP, snapOnVertexBorder);
			ptCentroidOnBorderCP.setZ(0);

			//modif features
			for (auto r_it = range.first; r_it != range.second; ++r_it) {
				mCPNear[r_it->second].setGeometry(ptCentroidOnBorderCP);
				_fsCP->modifyFeature(mCPNear[r_it->second]);
			}
		}

		for (std::list<std::string>::const_iterator lit = lCp2Delete.begin() ; lit != lCp2Delete.end() ; ++lit) {
			_fsCP->deleteFeature(*lit);
		}
	}
}

///
///
///
bool app::calcul::CFeatGenerationOp::_areMergeable(
	ign::feature::Feature const& feat1,
	ign::feature::Feature const& feat2,
	double distance
) const {
	bool areCollinear = _areCollinear(feat1.getGeometry().asLineString(), feat2.getGeometry().asLineString());
	bool areDistanceTypeCompatible = _areDistanceTypeCompatible(feat1, feat2, distance);
	
	return !areCollinear && areDistanceTypeCompatible;
}

///
///
///
bool app::calcul::CFeatGenerationOp::_areDistanceTypeCompatible(
	ign::feature::Feature const& feat1,
	ign::feature::Feature const& feat2,
	double distance
) const {
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string const typeName = themeParameters->getValue(FORM_OF_WAY).toString();
	double const distMergeCP = themeParameters->getValue(CP_MERGE_DIST_CP).toDouble();
	double const distMergeTractorCP = themeParameters->getValue(CP_MERGE_DIST_TRACTOR_CP).toDouble();

	std::string const& type1 = feat1.getAttribute(typeName).toString();
	std::string const& type2 = feat2.getAttribute(typeName).toString();

	bool isWalkwayOrTractor1 = _sFormwayValues4BigDist2Merge.find(type1) != _sFormwayValues4BigDist2Merge.end();
	bool isWalkwayOrTractor2 = _sFormwayValues4BigDist2Merge.find(type2) != _sFormwayValues4BigDist2Merge.end();

	return isWalkwayOrTractor1 && isWalkwayOrTractor2 ? distance < distMergeTractorCP : distance < distMergeCP;
}

///
///
///
bool app::calcul::CFeatGenerationOp::_areCollinear(
	ign::geometry::LineString const& ls1,
	ign::geometry::LineString const& ls2
) const {
	return false;

	ign::geometry::algorithm::OptimizedHausdorffDistanceOp op(ls1, ls2, -1, 10 /*TODO a rendre parametrable*/);
	double dAB = op.getDemiHausdorff(ign::geometry::algorithm::OptimizedHausdorffDistanceOp::DhdFromAtoB);
	if (dAB >= 0) return true;
	double dBA = op.getDemiHausdorff(ign::geometry::algorithm::OptimizedHausdorffDistanceOp::DhdFromBtoA);
	if (dBA < 0) return false;
	return true;
}

///
///
///
bool app::calcul::CFeatGenerationOp::_isEdgeConnected2cl(
	ign::geometry::Geometry & geomObjNearCl,
	ign::geometry::Envelope & envArroundGeom,
	ign::feature::Feature & fCl2SnapOn,
	double distMinCl
) const {
	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	ign::feature::FeatureFilter filterArroundCp(countryCodeName + " = '" + _countryCodeDouble + "'");
	filterArroundCp.setExtent(envArroundGeom);
	ign::feature::FeatureIteratorPtr itClArround = _fsEdge->getFeatures(filterArroundCp);

	bool hasEdgConnected2Cl = false;
	while (itClArround->hasNext())
	{
		ign::feature::Feature fClArround = itClArround->next();
		ign::geometry::LineString lsClArround = fClArround.getGeometry().asLineString();
		double dist = lsClArround.distance(geomObjNearCl);
		if (dist < distMinCl) {
			hasEdgConnected2Cl = true;
			distMinCl = dist;
			fCl2SnapOn = fClArround;
		}
	}
	return hasEdgConnected2Cl;
}

///
///
///
void app::calcul::CFeatGenerationOp::_snapCpOnClNearBy(
	double distCp2snapCl,
	double snapDistOnVertexFromCl,
	std::map< std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint> > & mClSplittedByCp
) const {
	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	ign::feature::FeatureFilter filterCp;
	for (size_t i = 0; i < _vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCp, countryCodeName + " = '" + _vCountriesCodeName[i] + "'");
	}
	ign::feature::FeatureIteratorPtr itCp = _fsCP->getFeatures(filterCp);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCP, filterCp);
	boost::progress_display display(numFeatures, std::cout, "[ SNAP CP ON CL]\n");
	std::map<std::string, ign::feature::Feature> mSnappedCpOnCl;
	while (itCp->hasNext())
	{
		++display;
		ign::feature::Feature fCp= itCp->next();
		ign::geometry::Point ptCp = fCp.getGeometry().asPoint();

		ign::feature::Feature fCl2SnapOn;
		double  distMinCl = distCp2snapCl * 2;
		bool hasClNearBy = _isEdgeConnected2cl(ptCp,ptCp.getEnvelope().expandBy(distCp2snapCl), fCl2SnapOn, distMinCl);

		/*ign::feature::FeatureFilter filterArroundCp(countryCodeName + " = '" + _countryCodeDouble+ "'");
		filterArroundCp.setExtent(ptCp.getEnvelope().expandBy(distCp2snapCl));
		ign::feature::FeatureIteratorPtr itClArround = _fsEdge->getFeatures(filterArroundCp);

		ign::feature::Feature fCl2SnapOn;
		double distMinCl = distCp2snapCl * 2;

		while (itClArround->hasNext())
		{
			ign::feature::Feature fClArround = itClArround->next();
			ign::geometry::LineString lsClArround = fClArround.getGeometry().asLineString();
			double dist = ptCp.distance(lsClArround);
			if (dist < distMinCl) {
				distMinCl = dist;
				fCl2SnapOn = fClArround;
			}
		}
		//pas de Cl proche pour snapper
		if (fCl2SnapOn.getId().empty())
			continue;*/
			
		//pas de Cl proche pour snapper
		if (!hasClNearBy)
			continue;

		//on projette la Cp sur la Cl identifiée
		ign::geometry::LineString lsCl2SnapOn = fCl2SnapOn.getGeometry().asLineString();
		ign::geometry::Point ptSnap = epg::tools::geometry::project(lsCl2SnapOn, ptCp, snapDistOnVertexFromCl);
		ptSnap.setZ(0);

		//on modifie la geom du Cp
		fCp.setGeometry(ptSnap);
		mSnappedCpOnCl[fCp.getId()] = fCp;

		//on stocke la geom du Cp comme coupure de la Cl
		std::string idCl2SnapOn = fCl2SnapOn.getId();
		if (mClSplittedByCp.find(idCl2SnapOn) == mClSplittedByCp.end()) { //pas de decoupe pour cette Cl
			ign::geometry::MultiPoint mpt;
			mpt.addGeometry(ptSnap);
			std::pair<ign::feature::Feature, ign::geometry::MultiPoint> pairClSnappedOn = std::make_pair(fCl2SnapOn, mpt);
			mClSplittedByCp[idCl2SnapOn] = pairClSnappedOn;
		}
		else { //il y a déjà des points de coupure sur cette CL
			std::pair<ign::feature::Feature, ign::geometry::MultiPoint> pairClSnappedOn = mClSplittedByCp.find(idCl2SnapOn)->second;
			pairClSnappedOn.second.addGeometry(ptSnap);
			mClSplittedByCp[idCl2SnapOn] = pairClSnappedOn;
		}
	}

	for (std::map<std::string, ign::feature::Feature>::iterator mit = mSnappedCpOnCl.begin(); mit != mSnappedCpOnCl.end(); ++mit) 
		_fsCP->modifyFeature(mit->second);	
}

///
///
///
void app::calcul::CFeatGenerationOp::_cutClByCp(
	std::map<std::string,std::pair<ign::feature::Feature, ign::geometry::MultiPoint>> & mClSplittedByCp
) const {
	std::set<std::string> sCl2delete;
	boost::progress_display display(mClSplittedByCp.size(), std::cout, "[ CUT CL BY CP]\n");
	for (std::map<std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint>>::iterator mit = mClSplittedByCp.begin();
		mit != mClSplittedByCp.end(); ++mit) {
		++display;
		ign::feature::Feature fEdgDoubl2Cut = mit->second.first;
		ign::geometry::LineString lsCl = fEdgDoubl2Cut.getGeometry().asLineString();
		app::geometry::tools::LineStringSplitter lsSplitterClSnappedOn(lsCl, 0.1); //on snappe à 10cm si il y a plusieurs coupures
		//boucle sur mpt
		ign::geometry::MultiPoint mpt = mit->second.second;
		for (size_t i = 0; i < mpt.numGeometries(); ++i)
			lsSplitterClSnappedOn.addCuttingGeometry(mpt.geometryN(i).asPoint());

		std::vector<ign::geometry::LineString> subCl2cut = lsSplitterClSnappedOn.getSubLineStringsZ();

		if (subCl2cut.size() == 1)
			continue;
			
		sCl2delete.insert(fEdgDoubl2Cut.getId());
		for (size_t i = 0; i < subCl2cut.size(); ++i) {
			ign::feature::Feature fNewEdgeDouble = fEdgDoubl2Cut;
			fNewEdgeDouble.setGeometry(subCl2cut[i]);
			fNewEdgeDouble.setId("");
			_fsEdge->createFeature(fNewEdgeDouble);
			_shapeLogger->writeFeature("edgeClCutByCp", fEdgDoubl2Cut);
			_logger->log(epg::log::INFO, "edge create by _cutClByCp: " + fNewEdgeDouble.getId());
		}
	}

	for (std::set<std::string>::iterator sit = sCl2delete.begin(); sit != sCl2delete.end(); ++sit)
		_fsEdge->deleteFeature(*sit);
}

///
///
///
bool app::calcul::CFeatGenerationOp::_getNearestCP(
	ign::feature::Feature const& fCP,
	double distMergeCP,
	std::map<std::string,ign::feature::Feature> & mCPNear
) const {
	epg::Context* context = epg::ContextS::getInstance();
	mCPNear[fCP.getId()] = fCP;
	
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	ign::feature::FeatureFilter filterArroundCP;
	for (std::map<std::string, ign::feature::Feature>::iterator mit = mCPNear.begin(); mit != mCPNear.end(); ++mit) {
		epg::tools::FilterTools::addAndConditions(filterArroundCP, idName + " <> '" + mit->first + "'");	//(idName + " <> '" + fCP.getId() + "'");
	}
	filterArroundCP.setExtent(fCP.getGeometry().getEnvelope().expandBy(distMergeCP));
	ign::feature::FeatureIteratorPtr itArroundCP = _fsCP->getFeatures(filterArroundCP);
	if (!itArroundCP->hasNext())
		return false;
	while (itArroundCP->hasNext())
	{
		ign::feature::Feature fCPArround = itArroundCP->next();
		_getNearestCP(fCPArround, distMergeCP, mCPNear);
		mCPNear[fCPArround.getId()] = fCPArround;
	}
	return true;
}

///
///
///
void  app::calcul::CFeatGenerationOp::_snapCl2Cl(
	double distMaxClClosest
) const {
	_logger->log(epg::log::TITLE, "[ BEGIN SNAP CONNECTING LINES 2 CONNECTING LINES ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	ign::feature::FeatureFilter filterCL(countryCodeName + " = '" + _countryCodeDouble + "'");

	ign::feature::FeatureIteratorPtr itCL = _fsCL->getFeatures(filterCL);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCL);
	boost::progress_display displayLoad(numFeatures, std::cout, "[ SNAP CL 2 CL]\n");

	GraphType graphCl;
	_loadGraphCL(graphCl);

	std::map<std::string,ign::feature::Feature> mFClModified;
	while (itCL->hasNext()) {
		++displayLoad;
		ign::feature::Feature fCL = itCL->next();
		ign::geometry::LineString lsCl = fCL.getGeometry().asLineString();

		edge_descriptor edCl = graphCl.getInducedEdges(fCL.getId()).second[0].descriptor;

		if(graphCl.degree(graphCl.source(edCl)) == 1) {
			ign::feature::Feature fCl2snapStart;
			bool isClosestStartCl2snapStart;
			bool hasClCandidate2snapStart = _hasClExtremityClose( distMaxClClosest, fCL, lsCl.startPoint(), fCl2snapStart, isClosestStartCl2snapStart);

			if (hasClCandidate2snapStart) {
				if (mFClModified.find(fCl2snapStart.getId()) != mFClModified.end()) //modif deja faite
					continue;
				ign::geometry::LineString lsCl2snapStart = fCl2snapStart.getGeometry().asLineString();
				ign::geometry::MultiPoint mp;
				mp.addGeometry(lsCl.startPoint());
				if(isClosestStartCl2snapStart)
					mp.addGeometry(lsCl2snapStart.startPoint());
				else 
					mp.addGeometry(lsCl2snapStart.endPoint());

				ign::geometry::Point pt2Modif = mp.asMultiPoint().getCentroid();
				pt2Modif.setZ(0);

				lsCl.setPointN(pt2Modif, 0);
				if (isClosestStartCl2snapStart)
					lsCl2snapStart.setPointN(pt2Modif, 0);
				else
					lsCl2snapStart.setPointN(pt2Modif, lsCl2snapStart.numPoints() - 1);

				fCL.setGeometry(lsCl);
				fCl2snapStart.setGeometry(lsCl2snapStart);
				mFClModified[fCL.getId()] = fCL;
				mFClModified[fCl2snapStart.getId()]= fCl2snapStart;
			}

		}
		if (graphCl.degree(graphCl.target(edCl)) == 1) {
			ign::feature::Feature fCl2snapEnd;
			bool isClosestStartCl2snapEnd;
			bool hasClCandidate2snapEnd = _hasClExtremityClose(distMaxClClosest, fCL, lsCl.endPoint(), fCl2snapEnd, isClosestStartCl2snapEnd);
			if (hasClCandidate2snapEnd) {
				if (mFClModified.find(fCl2snapEnd.getId()) != mFClModified.end()) //modif deja faite
					continue;
				ign::geometry::LineString lsCl2snapEnd = fCl2snapEnd.getGeometry().asLineString();
				ign::geometry::MultiPoint mp;
				mp.addGeometry(lsCl.endPoint());
				if (isClosestStartCl2snapEnd)
					mp.addGeometry(lsCl2snapEnd.startPoint());
				else
					mp.addGeometry(lsCl2snapEnd.endPoint());

				ign::geometry::Point pt2Modif = mp.asMultiPoint().getCentroid();
				pt2Modif.setZ(0);

				lsCl.setPointN(pt2Modif, lsCl.numPoints() - 1);
				if (isClosestStartCl2snapEnd)
					lsCl2snapEnd.setPointN(pt2Modif, 0);
				else
					lsCl2snapEnd.setPointN(pt2Modif, lsCl2snapEnd.numPoints() - 1);

				fCL.setGeometry(lsCl);
				fCl2snapEnd.setGeometry(lsCl2snapEnd);
				mFClModified[fCL.getId()] = fCL;
				mFClModified[fCl2snapEnd.getId()] = fCl2snapEnd;
			}
		}
	}

	for (std::map<std::string, ign::feature::Feature>::iterator mit = mFClModified.begin(); mit != mFClModified.end(); ++mit) {
		_fsCL->modifyFeature(mit->second);
	}

}

///
///
///
bool  app::calcul::CFeatGenerationOp::_hasClExtremityClose(
	double distMaxClClosest,
	ign::feature::Feature fClCurr,
	ign::geometry::Point ptClCurr,
	ign::feature::Feature & fCl2snap,
	bool isClosestStartCl2snap 
) const {
	//recup des CL proches si dist = 0 rien faire, si dist inf distMaxClClosest snapper
	epg::Context* context = epg::ContextS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	ign::feature::FeatureFilter filterArroundCl (idName + " <> '" + fClCurr.getId() + "'");
	epg::tools::FilterTools::addAndConditions(filterArroundCl, countryCodeName + " = '" + _countryCodeDouble + "'");
	filterArroundCl.setExtent(ptClCurr.getEnvelope().expandBy(distMaxClClosest));
	ign::feature::FeatureIteratorPtr itClArround = _fsCL->getFeatures(filterArroundCl);
	while (itClArround->hasNext()) {
		ign::feature::Feature fClArround = itClArround->next();
		ign::geometry::LineString lsClArround = fClArround.getGeometry().asLineString();
		double distClCurr = lsClArround.distance(fClCurr.getGeometry());

		if (distClCurr == 0 || distClCurr > distMaxClClosest)
			continue;

		double distStart = lsClArround.startPoint().distance(ptClCurr);
		double distEnd= lsClArround.endPoint().distance(ptClCurr);

		if(distEnd < distMaxClClosest && distStart < distMaxClClosest )
			_logger->log(epg::log::WARN, "CL has two extremities close : " + fClCurr.getId());

		if (distEnd < distStart && distEnd < distMaxClClosest ) {
			fCl2snap = fClArround;
			isClosestStartCl2snap = false;
			return true;
		}

		if (distStart < distEnd && distStart < distMaxClClosest) {
			fCl2snap = fClArround;
			isClosestStartCl2snap = true;
			return true;
		}
	}
	return false;
}

///
///
///
void app::calcul::CFeatGenerationOp::_mergeIntersectingClWithGraph(
	double distMaxEdges,
	double snapProjCl2edge
) const {
	_logger->log(epg::log::TITLE, "[ BEGIN FUSION CONNECTING LINES ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	std::string separator = "#";

	GraphType graphCl;
	ign::geometry::graph::tools::SnapRoundPlanarizer< GraphType >  planarizerCl(graphCl, 10000);//scale =100 -> precision de 0.01
	ign::feature::FeatureFilter filterCL;
	for (size_t i = 0; i < _vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCL, countryCodeName + " = '" + _vCountriesCodeName[i] + "'");
	}
	ign::feature::FeatureIteratorPtr itCL = _fsCL->getFeatures(filterCL);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCL);
	boost::progress_display displayLoad(numFeatures, std::cout, "[ LOAD GRAPH PLANARIZE CL ]\n");
	while (itCL->hasNext()) {
		++displayLoad;
		ign::feature::Feature fCL = itCL->next();
		planarizerCl.addEdge(fCL.getGeometry().asLineString(), fCL.getId());
	}
	planarizerCl.planarize();

	//fusion des edges si le même adjacent avec les mêmes edges d'origines
	_mergingEdgesByOrigin(graphCl);

	GraphType::edge_iterator eit, eitEnd;
	graphCl.edges(eit, eitEnd);
	boost::progress_display display(graphCl.numEdges(), std::cout, "[ FUSION CONNECTING LINES ]\n");
	while (eit != eitEnd) {
		++display;

		std::vector<std::string> vClOrigins = graphCl.origins(*eit);
		ign::geometry::LineString lsCl = graphCl.getGeometry(*eit);

		if (graphCl.origins(*eit).size() == 1) {
			++eit;
			continue;
		}

		std::map<std::string, ign::feature::Feature> mIdClOriginsCountry1, mIdClOriginsCountry2;

		//recuperation des edges liés aux CLs
		for (std::vector<std::string>::iterator vit = vClOrigins.begin(); vit != vClOrigins.end(); ++vit) {
			ign::feature::Feature fCl;
			_fsCL->getFeatureById(*vit, fCl);
			//_logger->log(epg::log::WARN, "CL DOES NOT EXIST  : " + *vit);
			if (fCl.getId().empty()) {
				//sClDoublonToDelete.insert(*vit);
				continue;
			}
			std::string countryCodeCl = fCl.getAttribute(countryCodeName).toString();
			if (countryCodeCl == _vCountriesCodeName[0])
				mIdClOriginsCountry1[fCl.getId()] = fCl;
			else if (countryCodeCl == _vCountriesCodeName[1])
				mIdClOriginsCountry2[fCl.getId()] = fCl;
			else //ne devrait pas arriver
				continue;
		}

		if (mIdClOriginsCountry1.size() == 0 || mIdClOriginsCountry2.size() == 0) {
			++eit;
			continue;//pas de fusion, CL de un seul pays
		}

		//recuperation des portions d'edges associées et selection des CLs à fusionner
		std::set<std::string> sEdgesMerged;
		int nbEdgesMerged= -1;

		while(nbEdgesMerged != sEdgesMerged.size() ){
			nbEdgesMerged = sEdgesMerged.size();

			double distMin = 100000;
			std::pair<std::string, std::string> cl2merge;

			for (std::map<std::string, ign::feature::Feature>::iterator mit1 = mIdClOriginsCountry1.begin(); mit1 != mIdClOriginsCountry1.end(); ++mit1) {

				std::string idEdgeLinked1 = mit1->second.getAttribute(linkedFeatIdName).toString();
				if (sEdgesMerged.find(mit1->first) != sEdgesMerged.end())//deja utilise pour merge
					continue;
				ign::feature::Feature fEdge1;
				_fsEdge->getFeatureById(idEdgeLinked1, fEdge1);
				ign::geometry::LineString lsEdge1 = fEdge1.getGeometry().asLineString();
				ign::geometry::LineString lsClEdge1;
				_getGeomProjClOnEdge(lsCl, lsEdge1, lsClEdge1, snapProjCl2edge);
				if (lsClEdge1.isEmpty())
					continue;

				for (std::map<std::string, ign::feature::Feature>::iterator mit2 = mIdClOriginsCountry2.begin(); mit2 != mIdClOriginsCountry2.end(); ++mit2) {

					std::string idEdgeLinked2 = mit2->second.getAttribute(linkedFeatIdName).toString();
					if (sEdgesMerged.find(mit2->first) != sEdgesMerged.end())//deja utilise pour merge
						continue;
					ign::feature::Feature fEdge2;
					_fsEdge->getFeatureById(idEdgeLinked2, fEdge2);
					ign::geometry::LineString lsEdge2 = fEdge2.getGeometry().asLineString();
					ign::geometry::LineString lsClEdge2;
					_getGeomProjClOnEdge(lsCl, lsEdge2, lsClEdge2, snapProjCl2edge);

					if (lsClEdge2.isEmpty())
						continue;

					double hausdorffDistEdges = ign::geometry::algorithm::OptimizedHausdorffDistanceOp::distance(lsClEdge1, lsClEdge2);

					if (hausdorffDistEdges > distMaxEdges) {

						ign::feature::Feature feat;
						std::ostringstream ss;
						ss << mit1->first << " : " << mit1->second.getAttribute(linkedFeatIdName).toString() << "   " << mit2->first << " : " << idEdgeLinked2;
						// feat.setAttribute("message",ign::data::String( ss.str() ) );
						feat.setGeometry(lsEdge2);
						_shapeLogger->writeFeature("ClDebug", feat, ss.str());

						continue;//on ne garde que les CL dont les edges associés sont proches (sous un seuil)
					}		

					if (hausdorffDistEdges < distMin) {
						distMin = hausdorffDistEdges;
						cl2merge = std::make_pair(mit1->first, mit2->first);
					}
				}
			}

			if (cl2merge.first.empty()) {
				continue;//break
			}
			//_logger->log(epg::log::DEBUG, cl2merge.first);
			//_logger->log(epg::log::DEBUG, cl2merge.second);

			sEdgesMerged.insert(cl2merge.first);
			sEdgesMerged.insert(cl2merge.second);

			ign::feature::Feature fClNew = _fsCL->newFeature();
			fClNew = mIdClOriginsCountry1.find(cl2merge.first)->second;

			//_logger->log(epg::log::DEBUG, *_attrMergerOnBorder.getAttrNameW().begin());

			_attrMergerOnBorder.mergeFeatAttribute(fClNew, mIdClOriginsCountry2.find(cl2merge.second)->second, separator);
			//_addFeatAttributeMergingOnBorder(fClNew, mIdClOriginsCountry2.find(cl2merge.second)->second, separator);

			std::string idCLNew = _idGeneratorCL->next();
			fClNew.setId(idCLNew);
			lsCl.setFillZ(0);
			fClNew.setGeometry(lsCl);

			_fsCL->createFeature(fClNew, idCLNew);
		}
		++eit;
	}

	//suppression des CL sans #
	std::ostringstream ss;
	ss << "DELETE  FROM " << _fsCL->getTableName() << " WHERE " << countryCodeName << " not like '%#%'";
	context->getDataBaseManager().getConnection()->update(ss.str());

	//fusion des CL ayant les mêmes linkedFeatIdName et qui se touchent (ou presque)
	// A ce stade il y a des CL avec simple linkedFeatIdName et double linkedFeatIdName ...!?
	std::ostringstream ss1;
	ss1 << linkedFeatIdName << " LIKE '%#%'";
	ign::feature::FeatureFilter mergeFilter(ss1.str());
	ome2::calcul::detail::ClMerger::mergeAll(_fsCL, mergeFilter, _idGeneratorCL.get());

	_logger->log(epg::log::TITLE, "[ END FUSION CONNECTING LINES ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
///!!!! PLUS UTILISER
void app::calcul::CFeatGenerationOp::_mergeIntersectingCL2(
	double distMergeCL,
	double snapOnVertexBorder
) const {
	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	double const distMaxFromBorder = themeParameters->getValue(CL_BUFFER_DIST).toDouble();

	ign::feature::FeatureFilter filterCL;
	for (size_t i = 0; i < _vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCL, countryCodeName + " = '" + _vCountriesCodeName[i] + "'");
	}
	//std::string const natIdName = themeParameters->getValue(NATIONAL_IDENTIFIER).toString();
	ign::feature::FeatureIteratorPtr itCL = _fsCL->getFeatures(filterCL);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCL);
	boost::progress_display display(numFeatures, std::cout, "[ FUSION CONNECTING LINES WITH #]\n");

	std::set<std::string> sCL2Merged;
	std::string separator = "#";

	while (itCL->hasNext())
	{
		++display;
		ign::feature::Feature fCLCurr = itCL->next();
		std::string idCLCurr = fCLCurr.getId();
		ign::geometry::LineString lsCurr = fCLCurr.getGeometry().asLineString();
		std::string countryCodeCLCurr = fCLCurr.getAttribute(countryCodeName).toString();
		//debug	
		/*if (fCLCurr.getAttribute(natIdName).toString() == "TRONROUT0000000057243619")
			bool bStop = true;
		else if (fCLCurr.getAttribute(natIdName).toString() == "{2005FE06-ED7A-4CA6-A011-471ADD678B46}")
			bool bStop = true;*/
		if (countryCodeCLCurr.find("#") != std::string::npos)
			continue;
		if (sCL2Merged.find(idCLCurr) != sCL2Merged.end())
			continue;
		sCL2Merged.insert(idCLCurr);

		ign::geometry::LineString lsBorder;
		_getBorderFromEdge(lsCurr, lsBorder);
		epg::tools::MultiLineStringTool mslBorder(lsBorder);

		ign::feature::FeatureFilter filterArroundCL;
		filterArroundCL.setPropertyConditions(countryCodeName + " != '" + countryCodeCLCurr + "'");
		filterArroundCL.setExtent(lsCurr.getEnvelope());
		ign::feature::FeatureIteratorPtr itArroundCL = _fsCL->getFeatures(filterArroundCL);
		while (itArroundCL->hasNext())
		{
			ign::feature::Feature fCLArround = itArroundCL->next();
			std::string idClArround = fCLArround.getId();
			ign::geometry::LineString lsClArround = fCLArround.getGeometry().asLineString();
			//si CL deja traite on ne fait rien
			if (sCL2Merged.find(idClArround) != sCL2Merged.end())
				continue;
			std::string countryCodeCLArround = fCLArround.getAttribute(countryCodeName).toString();
			if (countryCodeCLArround.find("#") != std::string::npos)
				continue;
			//si pas d'intersection avec une CL d'un autre pays on ne crée pas de CL merged
			if (lsCurr.distance(lsClArround) > 0.1 )
				continue;

			//si intersection uniquement aux extremites des CLs on ne merge pas
			if (lsCurr.Intersection(lsClArround)->isPoint()) {
				ign::geometry::Point ptIntersect = lsCurr.Intersection(lsClArround)->asPoint();
				if ((ptIntersect == lsCurr.startPoint() && ptIntersect == lsClArround.startPoint())
					|| (ptIntersect == lsCurr.startPoint() && ptIntersect == lsClArround.endPoint())
					|| (ptIntersect == lsCurr.endPoint() && ptIntersect == lsClArround.startPoint())
					|| (ptIntersect == lsCurr.endPoint() && ptIntersect == lsClArround.endPoint())
					) {
					ign::math::Vec2d vecLsCurr, vecLsArround;
					if (ptIntersect == lsCurr.startPoint()) {
						vecLsCurr.x() = lsCurr.endPoint().x() - lsCurr.startPoint().x();
						vecLsCurr.y() = lsCurr.endPoint().y() - lsCurr.startPoint().y();
					}
					else {
						vecLsCurr.x() = lsCurr.startPoint().x() - lsCurr.endPoint().x();
						vecLsCurr.y() = lsCurr.startPoint().y() - lsCurr.endPoint().y();
					}
					if (ptIntersect == lsClArround.startPoint()) {
						vecLsArround.x() = lsClArround.endPoint().x() - lsClArround.startPoint().x();
						vecLsArround.y() = lsClArround.endPoint().y() - lsClArround.startPoint().y();
					}
					else {
						vecLsArround.x() = lsClArround.startPoint().x() - lsClArround.endPoint().x();
						vecLsArround.y() = lsClArround.startPoint().y() - lsClArround.endPoint().y();
					}
					double anglLs = epg::tools::geometry::angle(vecLsCurr, vecLsArround);
					if (anglLs > 0.01) //si l'angle est faible, ça peut deux tronçons superposés, s'intersectant en un seul point pour erreur d'arrondi
						continue;
				}	
			}

			ign::geometry::LineString lsIntersectedCL, lsSE, lsSS, lsES, lsEE;

			ign::geometry::LineString segmentSE(lsCurr.startPoint(), lsClArround.endPoint());
			_getGeomCL(lsSE, mslBorder, segmentSE, distMaxFromBorder, snapOnVertexBorder);

			ign::geometry::LineString segmentSS(lsCurr.startPoint(), lsClArround.startPoint());
			_getGeomCL(lsSS, mslBorder, segmentSS, distMaxFromBorder, snapOnVertexBorder);

			ign::geometry::LineString segmentES(lsCurr.endPoint(), lsClArround.startPoint());
			_getGeomCL(lsES, mslBorder, segmentES, distMaxFromBorder, snapOnVertexBorder);

			ign::geometry::LineString segmentEE(lsCurr.endPoint(), lsClArround.endPoint());
			_getGeomCL(lsEE, mslBorder, segmentEE, distMaxFromBorder, snapOnVertexBorder);
			
			//lsIntersectedCL = lsCurr;
			double lengthMin = 100000;
			//debug
			double testLength = lsSE.length();
			if (lsSE.length() < lengthMin && lsSE.length()> 0.1 ) {//on s'assure que la section de frontière n'est pas nulle et que les projections ne se sont pas snappés au même endroit sur la frontière
				int numSegSE = static_cast<int>(std::floor(lsSE.numSegments() / 2.));
				ign::geometry::Point ptIntLsSE = epg::tools::geometry::interpolate(lsSE, numSegSE, 0.5);
				if (ptIntLsSE.distance(lsCurr) < 0.01 //on s'assure que la section est bien l'intersection et non le complement
					&& ptIntLsSE.distance(lsClArround) < 0.01) {
					lsIntersectedCL = lsSE;
					lengthMin = lsIntersectedCL.length();
				}
			}
			if (lsSS.length() < lengthMin && lsSS.length() > 0.1 ) {
				int numSegSS = static_cast<int>(std::floor(lsSS.numSegments() / 2.));
				ign::geometry::Point ptIntLsSS = epg::tools::geometry::interpolate(lsSS, numSegSS, 0.5);
				if (ptIntLsSS.distance(lsCurr) < 0.01 && ptIntLsSS.distance(lsClArround) < 0.01) {
					lsIntersectedCL = lsSS;
					lengthMin = lsIntersectedCL.length();
				}
			}
			if (lsES.length() < lengthMin && lsES.length() > 0.1){
				int numSegES = static_cast<int>(std::floor(lsES.numSegments() / 2.));
				ign::geometry::Point ptIntLsES = epg::tools::geometry::interpolate(lsES, numSegES, 0.5);
				if (ptIntLsES.distance(lsCurr) < 0.01 && ptIntLsES.distance(lsClArround) < 0.01) {
					lsIntersectedCL = lsES;
					lengthMin = lsIntersectedCL.length();
				}
			}
			if (lsEE.length() < lengthMin && lsEE.length() > 0.1 ){
				int numSegEE = static_cast<int>(std::floor(lsEE.numSegments() / 2.));
				ign::geometry::Point ptIntLsEE = epg::tools::geometry::interpolate(lsEE, numSegEE, 0.5);
				if (ptIntLsEE.distance(lsCurr) < 0.01 && ptIntLsEE.distance(lsClArround) < 0.01) {
					lsIntersectedCL = lsEE;
					lengthMin = lsIntersectedCL.length();
				}
			}
			//ajouter contient?
			//NON
			if (lsClArround.length() < lengthMin ) {
				lsIntersectedCL = lsClArround;
				lengthMin = lsIntersectedCL.length();
			}
			if (lsCurr.length() < lengthMin) {
				lsIntersectedCL = lsCurr;
				lengthMin = lsIntersectedCL.length();
			}
			//verifier si empty?

			//si la proj de l'intersection sur la frontiere se fait sur le meme point avec le snap on ne crée pas de CL
			if (lsIntersectedCL.numPoints() < 2)
				continue;		

			ign::feature::Feature fCLNew = _fsCL->newFeature();

			if (countryCodeCLArround < countryCodeCLCurr) {
				fCLNew = fCLArround;
				_attrMergerOnBorder.mergeFeatAttribute(fCLNew, fCLCurr, separator);
				//_addFeatAttributeMergingOnBorder(fCLNew, fCLCurr, separator);
			}
			else {
				fCLNew = fCLCurr;
				_attrMergerOnBorder.mergeFeatAttribute(fCLNew, fCLArround, separator);
				//_addFeatAttributeMergingOnBorder(fCLNew, fCLArround, separator);
			}

			std::string idCLNew = _idGeneratorCL->next();
			fCLNew.setId(idCLNew);
			lsIntersectedCL.setFillZ(0);
			fCLNew.setGeometry(lsIntersectedCL);

			fCLNew.setAttribute(themeParameters->getValue(CF_STATUS).toString(), ign::data::String("edge_matched"));

			_fsCL->createFeature(fCLNew, idCLNew);
			{
				ign::feature::Feature fShaplog = fCLNew;
				ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
				lsSphaplog.clearZ();
				fShaplog.setGeometry(lsSphaplog);
				_shapeLogger->writeFeature("ClMergedBeforeUpdate", fShaplog);
			}

		}
	}

	for (std::set<std::string>::iterator sit = sCL2Merged.begin(); sit != sCL2Merged.end(); ++sit) {
		_fsCL->deleteFeature(*sit);
	}

}

///
///
///
void app::calcul::CFeatGenerationOp::_deleteClByAngleAndDistEdges(
	double angleMax,
	double distMax,
	double snapProjCl2edge
) const {
	
	_logger->log(epg::log::TITLE, "[ BEGIN DELETE CL BY ANGLE EDGES FOR : " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + _countryCodeDouble + "'");

	int numCl2delete = -1;
	while (numCl2delete != 0){
		std::set<std::string> sCl2delete;
		GraphType graphCl;
		_loadGraphCL(graphCl);
		ign::feature::FeatureIteratorPtr it = _fsCL->getFeatures(filterCLCountryCode);
		int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCLCountryCode);
		boost::progress_display display(numFeatures, std::cout, "[  DELETE CL BY ANGLE EDGES ]\n");
		while (it->hasNext()) {
			++display;
			ign::feature::Feature fCl = it->next();

			_logger->log(epg::log::DEBUG, fCl.getId());

			//verifier si la cl n'est pas liée a au moins une cl à chaque extremite
			edge_descriptor edCl = graphCl.getInducedEdges(fCl.getId()).second[0].descriptor;

			if (graphCl.degree(graphCl.source(edCl)) > 1 && graphCl.degree(graphCl.target(edCl)) > 1)
				continue;

			_logger->log(epg::log::DEBUG, fCl.getAttribute(linkedFeatIdName).toString());

			//on verifie l'angle entre la projection de la cl sur les portions d'edges
			std::vector<std::string> vEdgeslinked;
			epg::tools::StringTools::Split(fCl.getAttribute(linkedFeatIdName).toString(), "#", vEdgeslinked);
			std::string idEdgLinked1 = vEdgeslinked[0];
			std::string idEdgLinked2 = vEdgeslinked[1];
			ign::feature::Feature fEdg1, fEdg2;
			_fsEdge->getFeatureById(idEdgLinked1, fEdg1);
			_fsEdge->getFeatureById(idEdgLinked2, fEdg2);
			if (fEdg1.getId().empty()) {
				_logger->log(epg::log::WARN, "Suppression CL  " + fCl.getId() + " not matching linked edge : " + idEdgLinked1);
				ign::feature::Feature fShaplog = fCl;
				ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
				fShaplog.setGeometry(lsSphaplog);
				_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fShaplog);
				sCl2delete.insert(fCl.getId());
				continue;
			}
			if (fEdg2.getId().empty()) {
				_logger->log(epg::log::WARN, "Suppression CL  " + fCl.getId() + " not matching linked edge : " + idEdgLinked2);
				ign::feature::Feature fShaplog = fCl;
				ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
				fShaplog.setGeometry(lsSphaplog);
				_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fShaplog);
				sCl2delete.insert(fCl.getId());
				continue;
			}
			ign::geometry::LineString lsEdg1 = fEdg1.getGeometry().asLineString();
			ign::geometry::LineString lsEdg2 = fEdg2.getGeometry().asLineString();
			ign::geometry::LineString lsCl = fCl.getGeometry().asLineString();
			ign::geometry::LineString lsProjClEdg1, lsProjClEdg2;
			_getGeomProjClOnEdge(lsCl, lsEdg1, lsProjClEdg1, snapProjCl2edge);
			_getGeomProjClOnEdge(lsCl, lsEdg2, lsProjClEdg2, snapProjCl2edge);
			if (lsProjClEdg1.isEmpty()) {
				_logger->log(epg::log::WARN, "Suppression CL  " + fCl.getId() + " not projecting on matching linked edge : " + idEdgLinked1);
				ign::feature::Feature fShaplog = fCl;
				ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
				fShaplog.setGeometry(lsSphaplog);
				_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fShaplog);
				sCl2delete.insert(fCl.getId());
				continue;
			}
			if (lsProjClEdg2.isEmpty()) {
				_logger->log(epg::log::WARN, "Suppression CL  " + fCl.getId() + "  not projecting on matching linked edge : " + idEdgLinked2);
				ign::feature::Feature fShaplog = fCl;
				ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
				fShaplog.setGeometry(lsSphaplog);
				_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fShaplog);
				sCl2delete.insert(fCl.getId());
				continue;
			}
			ign::math::Vec2d vec1(lsProjClEdg1.endPoint().x() - lsProjClEdg1.startPoint().x(), lsProjClEdg1.endPoint().y() - lsProjClEdg1.startPoint().y());
			ign::math::Vec2d vec2(lsProjClEdg2.endPoint().x() - lsProjClEdg2.startPoint().x(), lsProjClEdg2.endPoint().y() - lsProjClEdg2.startPoint().y());
			double angleEdgesLinked = epg::tools::geometry::angle(vec1, vec2);

			double hausdorffDist = ign::geometry::algorithm::OptimizedHausdorffDistanceOp::distance(lsProjClEdg1, lsProjClEdg2);

			//si angle trop important entre les deux edges on ne crée pas de Cl
			if (angleEdgesLinked > angleMax && angleEdgesLinked < (M_PI - angleMax) || hausdorffDist > distMax) {
				sCl2delete.insert(fCl.getId());
				{
					ign::feature::Feature fShaplog = fCl;
					ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
					lsSphaplog.clearZ();
					fShaplog.setGeometry(lsSphaplog);
					_shapeLogger->writeFeature("ClDeleteByAngleDistEdges", fShaplog);
				}
			}
		}

		for (std::set<std::string>::iterator sit = sCl2delete.begin(); sit != sCl2delete.end(); ++sit) 
			_fsCL->deleteFeature(*sit);

		
		numCl2delete = sCl2delete.size();
	}


	//_logger->log(epg::log::INFO, "Nb CL supprimées par angle des edges superieur a un seuil et non utile à la continuité : " + ign::data::Integer(sCl2delete.size()).toString());
	_logger->log(epg::log::TITLE, "[ END DELETE CL BY ANGLE EDGES FOR :" + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
bool app::calcul::CFeatGenerationOp::_getCLToMerge(
	ign::feature::Feature fCL,
	double distMergeCL,
	std::map < std::string, ign::feature::Feature> & mCL2merge,
	std::set<std::string> & sCountryCode
) const {
	epg::Context* context = epg::ContextS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	mCL2merge[fCL.getId()] = fCL;
	sCountryCode.insert(fCL.getAttribute(countryCodeName).toString());

	ign::feature::FeatureFilter filterArroundCL;
	for (std::map<std::string, ign::feature::Feature>::iterator mit = mCL2merge.begin(); mit != mCL2merge.end(); ++mit) {
		epg::tools::FilterTools::addAndConditions(filterArroundCL, idName + " <> '" + mit->first + "'");
	}
	filterArroundCL.setExtent(fCL.getGeometry().getEnvelope());
	ign::feature::FeatureIteratorPtr itArroundCL = _fsCL->getFeatures(filterArroundCL);
	if (!itArroundCL->hasNext())
		return false;
	while (itArroundCL->hasNext())
	{
		ign::feature::Feature fCLArround = itArroundCL->next();
		if (fCL.getGeometry().distance(fCLArround.getGeometry()) > distMergeCL)
			continue;
		std::string countryCodeCL = fCLArround.getAttribute(countryCodeName).toString();
		sCountryCode.insert(countryCodeCL);
		_getCLToMerge(fCLArround, distMergeCL, mCL2merge, sCountryCode);
		mCL2merge[fCLArround.getId()] = fCLArround;
	}
	return true;
}

///
///
///
void app::calcul::CFeatGenerationOp::_getBorderFromEdge(
	ign::geometry::LineString & lsEdgeOnBorder,
	ign::geometry::LineString & lsBorder
) const {
	ign::feature::FeatureFilter filter;//(countryCodeName + " = 'be#fr'")
	filter.setExtent(lsEdgeOnBorder.getEnvelope());
	ign::feature::FeatureIteratorPtr itBoundary = _fsBoundary->getFeatures(filter);
	while (itBoundary->hasNext()) {
		ign::feature::Feature fBorder = itBoundary->next();
		ign::geometry::LineString lsBorderTemp = fBorder.getGeometry().asLineString();
		double distBorder = lsEdgeOnBorder.distance(lsBorderTemp);
		if (distBorder < 0.1) {
			lsBorder = lsBorderTemp;
			return;
		}
	}
}

///
///
///
bool app::calcul::CFeatGenerationOp::_isNextEdgeInAntennas(
	ign::feature::Feature & fEdgeCurr,
	ign::geometry::Point & ptCurr,
	ign::feature::Feature &  edgeNext,
	ign::geometry::Point & ptNext
) const {
	epg::Context* context = epg::ContextS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string countryCodeCurr = fEdgeCurr.getAttribute(countryCodeName).toString();


	if (fEdgeCurr.getGeometry().asLineString().startPoint() == ptCurr)
		ptNext = fEdgeCurr.getGeometry().asLineString().endPoint();
	else
		ptNext = fEdgeCurr.getGeometry().asLineString().startPoint();

	ign::feature::FeatureFilter filterArroundNdNext(idName + " <> '" + fEdgeCurr.getId() + "' and " + countryCodeName + " ='" + countryCodeCurr + "'");
	filterArroundNdNext.setExtent(ptNext.getEnvelope().expandBy(0.01));
	ign::feature::FeatureIteratorPtr itNextEdge= context->getFeatureStore(epg::EDGE)->getFeatures(filterArroundNdNext);
	
	if (!itNextEdge->hasNext()) //fin de l'antenne par un cul de sac
		return false;
	edgeNext = itNextEdge->next();
	if (itNextEdge->hasNext()) //fin de l'antenne par un noeud de valence strictement sup a 2
		return false;
	if (edgeNext.getAttribute("is_intersected_border").toBoolean())//fin de l'antenne par intersection de la frontière
		return false;
	
	return true;
}

///
///
///
void app::calcul::CFeatGenerationOp::_updateGeomCL(double snapProjCl2edge) const
{
	_logger->log(epg::log::TITLE, "[ BEGIN UPDATE GEOM CL " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const edgeTableName = _fsEdge->getTableName();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	std::string const fictitiousFieldName = themeParameters->getValue(EDGE_FICTITIOUS).toString();
	//std::string const natIdName = themeParameters->getValue(NATIONAL_IDENTIFIER).toString();

	if(_vCountriesCodeName.size() != 2)
		_logger->log(epg::log::WARN, "Attention, le countryCode " + _countryCodeDouble + " n'a pas deux country" );

	std::string countryCode1 = _vCountriesCodeName[0];
	std::string countryCode2 = _vCountriesCodeName[1];

	//ign::geometry::MultiPolygon mPolyCountry1, mPolyCountry2;
	//_getGeomCountry(countryCode1, mPolyCountry1);
	//_getGeomCountry(countryCode2, mPolyCountry2);

	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + _countryCodeDouble + "'");
	ign::feature::FeatureIteratorPtr it = _fsCL->getFeatures(filterCLCountryCode);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCLCountryCode);
	boost::progress_display display(numFeatures, std::cout, "[  UPDATE GEOM CL ]\n");

	std::set<std::string> sCL2delete;
	while (it->hasNext()) {
		++display;
		ign::feature::Feature fCL = it->next();
		ign::geometry::LineString lsCLCurr = fCL.getGeometry().asLineString();
		lsCLCurr.setFillZ(0);
		ign::geometry::LineString lsCLUpdated;

		std::vector<std::string> vEdgeslinked;
		epg::tools::StringTools::Split(fCL.getAttribute(linkedFeatIdName).toString(), "#", vEdgeslinked);
		std::string idEdgLinked1 = vEdgeslinked[0];
		std::string idEdgLinked2 = vEdgeslinked[1];	
		ign::feature::Feature fEdg1, fEdg2;

		_fsEdge->getFeatureById(idEdgLinked1, fEdg1);
		_fsEdge->getFeatureById(idEdgLinked2, fEdg2);

		if (fEdg1.getId().empty() || fEdg2.getId().empty()) //si on ne trouve pas l'un des troncons liés
			continue;


		ign::geometry::LineString lsEdgInit1 = fEdg1.getGeometry().asLineString();
		ign::geometry::LineString lsEdgInit2 = fEdg2.getGeometry().asLineString();
		ign::geometry::LineString lsEdg1, lsEdg2;

		_getGeomProjClOnEdge(lsCLCurr, lsEdgInit1, lsEdg1, snapProjCl2edge);
		_getGeomProjClOnEdge( lsCLCurr, lsEdgInit2,lsEdg2, snapProjCl2edge);

		//on coupe les edges au niveau de la projection des extremites des CLs sur ces edges, pour ne prendre que la portion d'edge que l'on apparie à la CL

		if (lsEdg1.isEmpty()) {
			_logger->log(epg::log::WARN, "Suppression CL  " + fCL.getId() + " not matching linked edge : " + idEdgLinked1);
			ign::feature::Feature fShaplog = fCL; 
			ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
			lsSphaplog.clearZ();
			fShaplog.setGeometry(lsSphaplog);
			_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fShaplog);
			sCL2delete.insert(fCL.getId());
			continue;
		}

		if (lsEdg2.isEmpty()) {
			_logger->log(epg::log::WARN, "Suppression CL  " + fCL.getId() + " not matching linked edge : " + idEdgLinked2);
			ign::feature::Feature fShaplog = fCL;
			ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
			lsSphaplog.clearZ();
			fShaplog.setGeometry(lsSphaplog);
			_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fShaplog);
			sCL2delete.insert(fCL.getId());
			continue;
		}


		std::string ficticiousValue1 = fEdg1.getAttribute(fictitiousFieldName).toString();
		bool isFictEdg1 = false;
		if (ficticiousValue1 == "true")
			isFictEdg1 = true;
		std::string ficticiousClValue2 = fEdg2.getAttribute(fictitiousFieldName).toString();
		bool isFictEdg2 = false;
		if (ficticiousClValue2 == "true")
			isFictEdg2 = true;

		if (isFictEdg1 && !isFictEdg2)
			lsCLUpdated = lsEdg1;
		else if (!isFictEdg1 && isFictEdg2)
				lsCLUpdated = lsEdg2;
		else {

			//si les 2 edges sont dans le même pays, on prend la geom de la portion de l'edge du pays
			/*bool isLs1InCountry1 = lsEdg1.intersects(mPolyCountry1);
			bool isLs2InCountry1 = lsEdg2.intersects(mPolyCountry1);
			bool isLs1InCountry2 = lsEdg1.intersects(mPolyCountry2);
			bool isLs2InCountry2 = lsEdg2.intersects(mPolyCountry2);
			if (isLs1InCountry1 && !isLs1InCountry2 && isLs2InCountry1 && !isLs2InCountry2)
				_getGeomCL(lsCLUpdated, lsEdg1, lsCLCurr.startPoint(), lsCLCurr.endPoint(), snapOnVertex);
			else if (isLs1InCountry2 && !isLs1InCountry1 && isLs2InCountry2 && !isLs2InCountry1)
				_getGeomCL(lsCLUpdated, lsEdg2, lsCLCurr.startPoint(), lsCLCurr.endPoint(), snapOnVertex);
			else {*/

			std::set<double> sAbsCurv;
			geometry::tools::LengthIndexedLineString lsIndex1(lsEdg1);
			geometry::tools::LengthIndexedLineString lsIndex2(lsEdg2);
			for (size_t i = 0; i < lsEdg1.numPoints() - 1; ++i) {
				double abscurv = lsIndex1.getPointLocation(i) / lsEdg1.length();
				sAbsCurv.insert(abscurv);
			}
			for (size_t i = 0; i < lsEdg2.numPoints() - 1; ++i) {
				double abscurv = lsIndex2.getPointLocation(i) / lsEdg2.length();
				sAbsCurv.insert(abscurv);
			}

			for (std::set<double>::iterator sit = sAbsCurv.begin(); sit != sAbsCurv.end(); ++sit) {
				ign::geometry::MultiPoint multiPt;
				double test = *sit*lsEdg1.length();
				ign::geometry::Point pt1 = lsIndex1.locateAlong(*sit*lsEdg1.length());
				ign::geometry::Point pt2 = lsIndex2.locateAlong(*sit*lsEdg2.length());
				multiPt.addGeometry(pt1);
				multiPt.addGeometry(pt2);
				ign::geometry::Point ptLsCentroid = multiPt.getCentroid();
				bool hasPtDistMin = false;
				if (sit != sAbsCurv.begin()) {
					if (lsCLUpdated.endPoint().distance(ptLsCentroid) < 0)
						hasPtDistMin = true;
				}
				if (!hasPtDistMin)
					lsCLUpdated.addPoint(ptLsCentroid);
			}
			ign::geometry::MultiPoint multiPtEnd;
			multiPtEnd.addGeometry(lsEdg1.endPoint());
			multiPtEnd.addGeometry(lsEdg2.endPoint());
			//multiPtEnd.addGeometry(endLsProj2);
			lsCLUpdated.addPoint(multiPtEnd.getCentroid());
			//lsCLUpdated.clearZ();
			lsCLUpdated.setFillZ(0);
		}

		fCL.setGeometry(lsCLUpdated);
		_fsCL->modifyFeature(fCL);
		
	}

	for( std::set<std::string>::iterator sit = sCL2delete.begin(); sit != sCL2delete.end();++sit)
		_fsCL->deleteFeature(*sit);

	_logger->log(epg::log::TITLE, "[ BEGIN UPDATE GEOM CL " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

}

///
///
///
void app::calcul::CFeatGenerationOp::_getGeomProjClOnEdge(
	ign::geometry::LineString & lsCl,
	ign::geometry::LineString & lsEdge,
	ign::geometry::LineString & lsprojClOnEdg,
	double snapProjCl2edge
) const {
	ign::geometry::Point startLsProj = epg::tools::geometry::project(lsEdge, lsCl.startPoint(), snapProjCl2edge);
	startLsProj.setZ(0);//interpolate avec projZ
	ign::geometry::Point endLsProj = epg::tools::geometry::project(lsEdge, lsCl.endPoint(), snapProjCl2edge);
	endLsProj.setZ(0);//interpolate

	app::geometry::tools::LineStringSplitter lsSplitter1(lsEdge);
	lsSplitter1.addCuttingGeometry(startLsProj);
	lsSplitter1.addCuttingGeometry(endLsProj);
	std::vector<ign::geometry::LineString> vLsCandidate1 = lsSplitter1.getSubLineStringsZ();
	for (size_t i = 0; i < vLsCandidate1.size(); ++i) {
		ign::geometry::LineString lsCandidate = vLsCandidate1[i];
		if ((lsCandidate.startPoint().egal2d(startLsProj) && lsCandidate.endPoint().egal2d(endLsProj))
			|| (lsCandidate.startPoint().egal2d(endLsProj) && lsCandidate.endPoint().egal2d(startLsProj))) {
			lsprojClOnEdg = lsCandidate;
			lsprojClOnEdg.startPoint().setZ(0);
			lsprojClOnEdg.endPoint().setZ(0);
			break;
		}
	}

	if (lsprojClOnEdg.isEmpty())
		return;
	//on s'assure du sens de ls
	if (!lsprojClOnEdg.startPoint().egal2d(startLsProj))
		lsprojClOnEdg = lsprojClOnEdg.reverse();
}

///
///
///
void app::calcul::CFeatGenerationOp::_deleteCLUnderThreshold() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN CLEAN CL UNDER THRESHOLD " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	epg::Context* context = epg::ContextS::getInstance();
	std::string const CLTableName = _fsCL->getTableName();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	std::ostringstream ssconditionDeleteCLUnderThreshold;

	ssconditionDeleteCLUnderThreshold << " ST_LENGTH(" << geomName << ") < 10 AND "<<countryCodeName << " ='" << _countryCodeDouble << "'";

	std::set<std::string> sCLToDelete;
	ign::feature::FeatureFilter filterCLInf10m(ssconditionDeleteCLUnderThreshold.str());
	ign::feature::FeatureIteratorPtr it = _fsCL->getFeatures(filterCLInf10m);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCLInf10m);
	boost::progress_display display(numFeatures, std::cout, "[ CLEAN CL UNDER THRESHOLD ]\n");

	while (it->hasNext()) {
		++display;
		ign::feature::Feature fCL10m = it->next();
		ign::feature::FeatureFilter filterNeighbor(idName + " <> '" + fCL10m.getId()+"'");
		filterNeighbor.setExtent(fCL10m.getGeometry().getEnvelope());
		ign::feature::FeatureIteratorPtr itArround = _fsCL->getFeatures(filterNeighbor);
		bool hasNeithbor = false;
		while (itArround->hasNext()) {
			ign::feature::Feature fCLNeighbor = itArround->next();
			if (fCLNeighbor.getGeometry().distance(fCL10m.getGeometry()) < 0.1) {
				hasNeithbor = true;
				break;
			}
		}
		if(!hasNeithbor)
			sCLToDelete.insert(fCL10m.getId());
	}

	for (std::set<std::string>::iterator sit = sCLToDelete.begin(); sit != sCLToDelete.end(); ++sit) {
		_fsCL->deleteFeature(*sit);
	}


	_logger->log(epg::log::INFO, "Nb CL isolées inférieures a un seuil supprimées  " + ign::data::Integer(sCLToDelete.size()).toString());
	_logger->log(epg::log::TITLE, "[ END CLEAN CL UNDER THRESHOLD " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

}

///
///
///
void app::calcul::CFeatGenerationOp::_getGeomCountry(
	std::string countryCodeSimple,
	ign::geometry::MultiPolygon & geomCountry
) const {
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string const landCoverTypeName = themeParameters->getValue(LAND_COVER_TYPE).toString();
	std::string const landAreaValue = themeParameters->getValue(TYPE_LAND_AREA).toString();


	ign::feature::FeatureIteratorPtr itLandmask = _fsLandmask->getFeatures(ign::feature::FeatureFilter(landCoverTypeName + " = '" + landAreaValue + "' AND " + countryCodeName + " = '" + countryCodeSimple + "'"));

	while (itLandmask->hasNext()) {
		ign::feature::Feature const& fLandmask = itLandmask->next();
		ign::geometry::MultiPolygon const& mp = fLandmask.getGeometry().asMultiPolygon();
		for (int i = 0; i < mp.numGeometries(); ++i) {
			geomCountry.addGeometry(mp.polygonN(i));
		}
	}
}

///
///
///
void app::calcul::CFeatGenerationOp::_loadGraphCL(GraphType & graphCL) const
{
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + _countryCodeDouble + "'");
	ign::feature::FeatureIteratorPtr it = _fsCL->getFeatures(filterCLCountryCode);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCLCountryCode);
	boost::progress_display display(numFeatures, std::cout, "[ LOAD GRAPH CL ]\n");
	ign::geometry::graph::builder::SimpleGraphBuilder<GraphType> graphBuilder(graphCL, 0.01);
	while (it->hasNext()) {
		++display;
		ign::feature::Feature fCL = it->next();
		graphBuilder.addEdge(fCL.getGeometry().asLineString(), fCL.getId());
	}
}

///
///
///
void app::calcul::CFeatGenerationOp::_loadGraphEdges(
	std::string countryCodeSimple,
	GraphType & graphEdges
) const {
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	ign::feature::FeatureFilter filterEdgeCountryCode(countryCodeName + " = '" + countryCodeSimple + "'");
	ign::feature::FeatureIteratorPtr it = _fsEdge->getFeatures(filterEdgeCountryCode);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsEdge, filterEdgeCountryCode);
	boost::progress_display display(numFeatures, std::cout, "[ LOAD GRAPH EDGE " + countryCodeSimple+" ]\n");
	ign::geometry::graph::builder::SimpleGraphBuilder<GraphType> graphBuilder(graphEdges, 0.01);
	while (it->hasNext()) {
		++display;
		ign::feature::Feature fedge= it->next();
		graphBuilder.addEdge(fedge.getGeometry().asLineString(), fedge.getId());
	}
}

///
///
///
bool app::calcul::CFeatGenerationOp::_isConnectedEdges(
	GraphType& graph,
	std::string const& idEdge1,
	std::string const& idEdge2
) const {
	// _logger->log(epg::log::DEBUG, "log100");
	_logger->log(epg::log::DEBUG, idEdge1);
	edge_descriptor edCl1 = graph.getInducedEdges(idEdge1).second[0].descriptor;
	// _logger->log(epg::log::DEBUG, "log101");
	_logger->log(epg::log::DEBUG, idEdge2);
	edge_descriptor edCl2 = graph.getInducedEdges(idEdge2).second[0].descriptor;
	// _logger->log(epg::log::DEBUG, "log102");

	if (graph.source(edCl1) == graph.source(edCl2) || graph.source(edCl1) == graph.target(edCl2) || graph.target(edCl1) == graph.source(edCl2) || graph.target(edCl1) == graph.target(edCl2))
		return true;

	return false;
}

///
///
///
std::pair<bool,std::pair<std::string, std::string>> app::calcul::CFeatGenerationOp::_getClLinkedEdges(
	std::string const& linkedFeatIdName,
	GraphType & graphCL,
	GraphType::edge_descriptor eCl
) const {
	std::string idCl = graphCL.origins(eCl)[0];
	_logger->log(epg::log::DEBUG, idCl);
	ign::feature::Feature clFeat;
	_fsCL->getFeatureById(idCl, clFeat);
	if (clFeat.getId().empty())
		return std::make_pair(false, std::make_pair("", ""));
	std::vector<std::string> vEdgeslinkedJ;
	std::string linkedFeat = clFeat.getAttribute(linkedFeatIdName).toString();
	_logger->log(epg::log::DEBUG, linkedFeat);
	epg::tools::StringTools::Split(linkedFeat, "#", vEdgeslinkedJ);
	return std::make_pair(true,std::make_pair(vEdgeslinkedJ[0], vEdgeslinkedJ[1]));
}

///
///
///
bool app::calcul::CFeatGenerationOp::_areParallelEdges(
	GraphType& graphCL,
	GraphType::edge_descriptor e1,
	GraphType::edge_descriptor e2 
) const {
	if (graphCL.source(e1) == graphCL.source(e2) && graphCL.target(e1) == graphCL.target(e2) ) return true;
	if (graphCL.target(e1) == graphCL.source(e2) && graphCL.source(e1) == graphCL.target(e2) ) return true;
	return false;
}

///
///
///
ign::geometry::Point app::calcul::CFeatGenerationOp::_getLinkedEdgesConnectingPoint(
	GraphType const& graph, 
	std::string const& idEdge1,
	std::string const& idEdge2
) const {
	edge_descriptor edCl1 = graph.getInducedEdges(idEdge1).second[0].descriptor;
	edge_descriptor edCl2 = graph.getInducedEdges(idEdge2).second[0].descriptor;

	ign::feature::Feature fEdge1;
	_fsEdge->getFeatureById(idEdge1, fEdge1);
	ign::geometry::LineString const& edgeGeom1 = fEdge1.getGeometry().asLineString();

	if ( graph.source(edCl1) == graph.source(edCl2) || graph.source(edCl1) == graph.target(edCl2) )
		return edgeGeom1.startPoint();

	return edgeGeom1.endPoint();
}

///
///
///
void app::calcul::CFeatGenerationOp::_setContinuityCl(GraphType& graphCL) const
{
	_logger->log(epg::log::TITLE, "[ BEGIN SET CL CONTINUITY ] : " + epg::tools::TimeTools::getTime());

	epg::Context* context = epg::ContextS::getInstance();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	GraphType graphEdges1, graphEdges2;
	_loadGraphEdges(_vCountriesCodeName[0], graphEdges1);
	_loadGraphEdges(_vCountriesCodeName[1], graphEdges2);

	GraphType::vertex_iterator vit, vitEnd;
	graphCL.vertices(vit, vitEnd);

	boost::progress_display display(graphCL.numVertices(), std::cout, "[ SET CL CONTINUITY ]\n");
	std::map<std::string, ign::feature::Feature> mClModified;

	while (vit != vitEnd) {
		++display;
		if (graphCL.degree(*vit) < 2 ) {
			++vit;
			continue;
		}
		std::vector< GraphType::oriented_edge_descriptor > vClsIncidentTemp;
		graphCL.incidentEdges(*vit, vClsIncidentTemp);

		//recal
		std::vector<std::vector< GraphType::oriented_edge_descriptor >> vVClsTrueIncident;
		std::set<size_t> sTreated;
		for (size_t i = 0; i < vClsIncidentTemp.size()-1; ++i) {
			if ( sTreated.find(i) != sTreated.end() ) continue;

			std::vector< GraphType::oriented_edge_descriptor > vClsIncidentTempConnectI;
			vClsIncidentTempConnectI.push_back(vClsIncidentTemp[i]);

			std::pair<bool,std::pair<std::string, std::string>> pLinkedEdgesI = _getClLinkedEdges(linkedFeatIdName, graphCL, vClsIncidentTemp[i].descriptor);

			//si CL n'existe plus
			if (!pLinkedEdgesI.first)
				continue;


			for (size_t j = i + 1; j < vClsIncidentTemp.size(); ++j) {
				if ( sTreated.find(j) != sTreated.end() ) continue;
				std::pair<bool, std::pair<std::string, std::string>> pLinkedEdgesJ = _getClLinkedEdges(linkedFeatIdName, graphCL, vClsIncidentTemp[j].descriptor);

				if (!pLinkedEdgesJ.first)
					continue;

				bool isConnected1 = pLinkedEdgesI.first == pLinkedEdgesJ.first || _isConnectedEdges(graphEdges1, pLinkedEdgesI.second.first, pLinkedEdgesJ.second.first);
				bool isConnected2 = pLinkedEdgesI.second == pLinkedEdgesJ.second || _isConnectedEdges(graphEdges2, pLinkedEdgesI.second.second, pLinkedEdgesJ.second.second);

				// tester si vClsIncidentTemp[i].descriptor et vClsIncidentTemp[j].descriptor sont paralleles
				if (_areParallelEdges(graphCL, vClsIncidentTemp[i].descriptor, vClsIncidentTemp[j].descriptor)) {

					// si oui regarder si *vit est plus proche que l'autre sommet des points d'intersections des linkedEdges
					std::string idClI = graphCL.origins(vClsIncidentTemp[i].descriptor)[0];
					ign::feature::Feature fClI;
					_fsCL->getFeatureById(idClI, fClI);
					ign::geometry::LineString const& lsI = fClI.getGeometry().asLineString();

					ign::geometry::Point const& vitGeom = vClsIncidentTemp[i].direction == ign::graph::DIRECT ? lsI.startPoint() : lsI.endPoint();
					ign::geometry::Point const& otherGeom = vClsIncidentTemp[i].direction == ign::graph::DIRECT ? lsI.endPoint() : lsI.startPoint();


					ign::geometry::MultiPoint mpLinkedEdgesConnectingPoints;
					mpLinkedEdgesConnectingPoints.addGeometry(_getLinkedEdgesConnectingPoint(graphEdges1, pLinkedEdgesI.second.first, pLinkedEdgesJ.second.first));
					mpLinkedEdgesConnectingPoints.addGeometry(_getLinkedEdgesConnectingPoint(graphEdges2, pLinkedEdgesI.second.second, pLinkedEdgesJ.second.second));
					ign::geometry::Point linkedEdgesConnectingPoint = mpLinkedEdgesConnectingPoints.getCentroid();

					// si il ne l'est pas (et que les linkeEdges ne sont pas egalement paralleles?) --> continue
					if (linkedEdgesConnectingPoint.distance(vitGeom) > linkedEdgesConnectingPoint.distance(otherGeom)) continue;
				}

				if (isConnected1 && isConnected2) {
					vClsIncidentTempConnectI.push_back(vClsIncidentTemp[j]);
					sTreated.insert(j);
				}

			}
			if (vClsIncidentTempConnectI.size() > 1)
				vVClsTrueIncident.push_back(vClsIncidentTempConnectI);
		}

		//on recalcule la nouvelle geometrie du point
		for (size_t j = 0; j < vVClsTrueIncident.size(); ++j) {
			std::vector< GraphType::oriented_edge_descriptor > vClsTrueIncident = vVClsTrueIncident[j];
			ign::geometry::MultiPoint multiPtToConnect;
			for (size_t i = 0; i < vClsTrueIncident.size(); ++i) {
				GraphType::edge_descriptor edCl = vClsTrueIncident[i].descriptor;
				std::string idClToModify = graphCL.origins(edCl)[0];
				ign::feature::Feature fClToModify;
				if (mClModified.find(idClToModify) != mClModified.end()) {
					fClToModify = mClModified.find(idClToModify)->second;
				} else {
					_fsCL->getFeatureById(idClToModify, fClToModify);
					//patch tant que les doublons ne sont pas suppr
					if (fClToModify.getId().empty())
						continue;
					mClModified[idClToModify] = fClToModify;
				}

				ign::geometry::LineString ls = fClToModify.getGeometry().asLineString();
				if (graphCL.source(edCl) == *vit)
					multiPtToConnect.addGeometry(ls.startPoint());
				else
					multiPtToConnect.addGeometry(ls.endPoint());

			}
			ign::geometry::Point ptUpdated = multiPtToConnect.asMultiPoint().getCentroid();
			//on modifie la geom des cl avec celle du nouveau point
			for (size_t i = 0; i < vClsTrueIncident.size(); ++i) {
				GraphType::edge_descriptor edCl = vClsTrueIncident[i].descriptor;
				std::string idClToModify = graphCL.origins(edCl)[0];
				// patch tant que les doublons ne sont pas suppr
				if (mClModified.find(idClToModify) == mClModified.end())
					continue;
				ign::feature::Feature fClToModify = mClModified.find(idClToModify)->second;
				ign::geometry::LineString lsClToModify = fClToModify.getGeometry().asLineString();
				if (graphCL.source(edCl) == *vit)
					lsClToModify.setPointN(ptUpdated, 0);
				else
					lsClToModify.setPointN(ptUpdated, lsClToModify.numPoints() - 1);

				lsClToModify.setFillZ(0);
				fClToModify.setGeometry(lsClToModify);
				mClModified[idClToModify] = fClToModify;
			}
		}
		++vit;
	}

	//on modifie les feature dans postgis et on log
	for(std::map<std::string, ign::feature::Feature>::iterator mit = mClModified.begin(); mit != mClModified.end();++mit ) {
		_fsCL->modifyFeature(mit->second);
	}
	_logger->log(epg::log::TITLE, "[ END SET CL CONTINUITY ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::_getClDoublonGeom() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN DELETE CL DOUBLON ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	std::set<std::string> sClDoublonToDelete;

	GraphType graphClDoublon;
	ign::geometry::graph::tools::SnapRoundPlanarizer< GraphType >  planarizerClDoublon(graphClDoublon,100);//scale =100 -> precision de 0.01
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + _countryCodeDouble + "'");
	ign::feature::FeatureIteratorPtr it = _fsCL->getFeatures(filterCLCountryCode);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCLCountryCode);
	boost::progress_display displayLoad(numFeatures, std::cout, "[ LOAD GRAPH PLANARIZE CL ]\n");
	while (it->hasNext()) {
		++displayLoad;
		ign::feature::Feature fCL = it->next();
		planarizerClDoublon.addEdge(fCL.getGeometry().asLineString(), fCL.getId());
	}
	planarizerClDoublon.planarize();

	GraphType::edge_iterator eit, eitEnd;
	graphClDoublon.edges(eit, eitEnd);
	boost::progress_display display(graphClDoublon.numEdges(), std::cout, "[ DELETE CL DOUBLON ]\n");
	while (eit != eitEnd) {
		++display;
		
		if (graphClDoublon.origins(*eit).size() == 1) {
			++eit;
			continue;
		}

		std::vector<std::string> vIdClDoublon = graphClDoublon.origins(*eit);
		std::map <std::string, double > mIdClDoublonsAngle;
		for (std::vector<std::string>::iterator vit = vIdClDoublon.begin(); vit != vIdClDoublon.end(); ++vit) {
			ign::feature::Feature fClDoublon;
			_fsCL->getFeatureById(*vit, fClDoublon);
			_logger->log(epg::log::WARN, "CL DOUBLON  : " + *vit);
			if (fClDoublon.getId().empty() ) { 
				//sClDoublonToDelete.insert(*vit);
				continue;
			}
			//ign::geometry::LineString lsClDoublon = fClDoublon.getGeometry().asLineString();
			ign::feature::Feature fShaplog = fClDoublon;
			ign::geometry::LineString lsSphaplog = fShaplog.getGeometry().asLineString();
			lsSphaplog.clearZ();
			fShaplog.setGeometry(lsSphaplog);
			_shapeLogger->writeFeature("ClDoublon", fShaplog);
			
			/*
			std::vector<std::string> vEdgeslinked;
			epg::tools::StringTools::Split(fClDoublon.getAttribute(linkedFeatIdName).toString(), "#", vEdgeslinked);
			std::string idEdgLinked1 = vEdgeslinked[0];
			std::string idEdgLinked2 = vEdgeslinked[1];
			ign::feature::Feature fEdg1, fEdg2;
			_fsEdge->getFeatureById(idEdgLinked1, fEdg1);
			_fsEdge->getFeatureById(idEdgLinked2, fEdg2);
			if (fEdg1.getId().empty() || fEdg2.getId().empty()) { //si on ne trouve pas l'un des troncons liés
				sClDoublonToDelete.insert(*vit);
				continue;
			}
			ign::geometry::LineString lsEdg1 = fEdg1.getGeometry().asLineString();
			ign::geometry::LineString lsEdg2 = fEdg2.getGeometry().asLineString();
			//calcul de l'angle
			ign::math::Vec2d vecLsDoublon(lsClDoublon.endPoint().x() - lsClDoublon.startPoint().x(), lsClDoublon.endPoint().y() - lsClDoublon.startPoint().y());
			ign::math::Vec2d vecEdg1(lsEdg1.endPoint().x() - lsEdg1.startPoint().x(), lsEdg1.endPoint().y() - lsEdg1.startPoint().y());
			ign::math::Vec2d vecEdg2(lsEdg2.endPoint().x() - lsEdg2.startPoint().x(), lsEdg2.endPoint().y() - lsEdg2.startPoint().y());
			//, vecLsArround;
			double anglClEdg1 = epg::tools::geometry::angle(vecLsDoublon, vecEdg1);
			double anglClEdg2 = epg::tools::geometry::angle(vecLsDoublon, vecEdg2);
			//if(anglClEdg1)*/
		}

		++eit;
	}
	_logger->log(epg::log::TITLE, "[ END DELETE CL DOUBLON ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::_mergingEdgesByOrigin(
	GraphType & graph
) const {
	std::map< typename GraphType::edge_descriptor, std::set< typename GraphType::edge_descriptor >  > mGatheredEdges;

	boost::progress_display display(graph.numEdges(), std::cout, "[ gathergByOrigin % complete ]\n");

	std::set< GraphType::edge_descriptor > sVisitedEdges;

	GraphType::edge_iterator eit, eend;
	for (graph.edges(eit, eend); eit != eend; ++eit) {
		++display;
		if (graph.target(*eit) == graph.source(*eit)) continue;
		if (sVisitedEdges.find(*eit) != sVisitedEdges.end()) continue;

		double lengthPivot = graph.getGeometry(*eit).length();

		GraphType::edge_descriptor edgePivot = *eit;

		sVisitedEdges.insert(*eit);

		GraphType::oriented_edge_descriptor tPivot[] = {
			GraphType::oriented_edge_descriptor(*eit, ign::graph::DIRECT),
			GraphType::oriented_edge_descriptor(*eit, ign::graph::REVERSE)
		};

		//liste des arcs a fusionner
		std::set< GraphType::edge_descriptor > sEdges;

		bool isLoop = false;
		for (size_t i = 0; i < 2; ++i)
		{
			GraphType::oriented_edge_descriptor nextEdge = tPivot[i];

			GraphType::vertex_descriptor vTarget = graph.target(nextEdge);
			GraphType::vertex_descriptor vStart = graph.source(nextEdge);

			if (graph.degree(vTarget) != 2) continue;

			while (true) {

				if (vTarget == vStart) {
					isLoop = true;
					break;
				}

				std::vector< GraphType::oriented_edge_descriptor > vIncidentEdges;

				graph.incidentEdges(vTarget, vIncidentEdges);
				nextEdge = (vIncidentEdges.front().descriptor == nextEdge.descriptor) ? vIncidentEdges.back() : vIncidentEdges.front();

				double lengthNextEdge = graph.getGeometry(nextEdge.descriptor).length();

				if (graph.origins(edgePivot)!= graph.origins(nextEdge.descriptor)) break;

				IGN_ASSERT(sVisitedEdges.find(nextEdge.descriptor) == sVisitedEdges.end());

				sVisitedEdges.insert(nextEdge.descriptor);

				if (lengthNextEdge > lengthPivot) {
					sEdges.insert(edgePivot);
					//on change le pivot
					edgePivot = nextEdge.descriptor;
					lengthPivot = lengthNextEdge;
				}
				else
					sEdges.insert(nextEdge.descriptor);

				vTarget = graph.target(nextEdge);

				if (graph.degree(vTarget) != 2) break;

			}
			if (isLoop) break;
		}

		if (!sEdges.empty() && !isLoop)
		{
			mGatheredEdges.insert(std::make_pair(edgePivot, sEdges));
		}
	}

	//patience
	boost::progress_display display2(mGatheredEdges.size(), std::cout, "[ mergingByOrigin % complete ]\n");

	typename std::map< GraphType::edge_descriptor, std::set<GraphType::edge_descriptor> >::iterator mit;
	for (mit = mGatheredEdges.begin(); mit != mGatheredEdges.end(); ++mit, ++display2)
	{
		//on ordonne les arcs a fusionner
		GraphType::edges_path path = epg::graph::tools::createPath(graph, mit->first, mit->second);

		if (graph.source(path.begin()->descriptor) == graph.target(path.rbegin()->descriptor)) 
			continue;
		
		//on fusionne
		//ign::geometry::LineString lsResult;
		std::vector<ign::geometry::Point> vGeomNewE;
		for (GraphType::edges_path_iterator it = path.begin(); it != path.end(); ++it) {

			ign::geometry::LineString lsTemp = graph.getGeometry(it->descriptor);
			if (it == path.begin()) vGeomNewE.push_back((it->direction == ign::graph::DIRECT) ? lsTemp.startPoint() : lsTemp.endPoint());
			if (it->direction == ign::graph::DIRECT)
			{
				ign::geometry::LineString::iterator itLs = lsTemp.begin();
				for (++itLs; itLs != lsTemp.end(); ++itLs) {
					vGeomNewE.push_back(*itLs);
				}
			}
			else {
				ign::geometry::LineString::reverse_iterator itLs = lsTemp.rbegin();
				for (++itLs; itLs != lsTemp.rend(); ++itLs) {
					vGeomNewE.push_back(*itLs);
				}
			}
		}

		GraphType::vertex_descriptor vSource = graph.source(*path.begin());
		GraphType::vertex_descriptor vTarget = graph.target(*path.rbegin());

		GraphType::oriented_edge_descriptor newD = graph.addEdge(vSource, vTarget, vGeomNewE, graph[mit->first]);

		for (GraphType::edges_path_iterator it = path.begin(); it != path.end(); ++it)
			if (it->descriptor == mit->first)
			{
				it->descriptor = newD.descriptor;
				break;
			}

		std::set< GraphType::vertex_descriptor > vVerticesToRemove;
		for (GraphType::edges_path_iterator it = path.begin(); it != path.end(); ) {
			typename GraphType::edges_path_iterator itTemp = it++;

			GraphType::vertex_descriptor vS = graph.source(*itTemp);
			GraphType::vertex_descriptor vT = graph.target(*itTemp);

			if (itTemp->descriptor == newD.descriptor) continue;

			graph.removeEdge(itTemp->descriptor);//persistent by default for Feature Graph

			if (itTemp != path.begin())
				vVerticesToRemove.insert(vS);
			if (itTemp->descriptor != path.rbegin()->descriptor)
				vVerticesToRemove.insert(vT);
		}

		typename std::set< GraphType::vertex_descriptor >::iterator sit;
		for (sit = vVerticesToRemove.begin(); sit != vVerticesToRemove.end(); ++sit)
			graph.removeVertex(*sit);//persistent by default for Feature Graph
		
	}
}
