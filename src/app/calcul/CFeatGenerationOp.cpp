// APP
#include <app/calcul/CFeatGenerationOp.h>
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LengthIndexedLineString.h>
#include <app/geometry/tools/LineStringSplitter.h>

// BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>

// SOCLE
#include <ign/geometry/algorithm/BufferOpGeos.h>
#include <ign/geometry/graph/builder/SimpleGraphBuilder.h>
#include <ign/geometry/graph/tools/SnapRoundPlanarizer.h>
#include <ign/geometry/algorithm/OptimizedHausdorffDistanceOp.h>

// EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/tools/FilterTools.h>
#include <epg/tools/StringTools.h>
#include <epg/utils/CopyTableUtils.h>
#include <epg/tools/geometry/interpolate.h>
#include <epg/tools/geometry/project.h>
#include <epg/tools/geometry/getSubLineString.h>
#include <epg/tools/geometry/angle.h>
#include <epg/tools/geometry/SegmentIndexedGeometry.h>
#include <epg/tools/geometry/LineIntersector.h>
#include<epg/graph/tools/merge.h>
#include<epg/graph/tools/createPath.h>

// OME2
#include <ome2/calcul/detail/ClMerger.h>
#include <ome2/feature/sql/NotDestroyedTools.h>


///
///
///
app::calcul::CFeatGenerationOp::CFeatGenerationOp(
	std::string const& countryCodeDouble,
	bool verbose
) {
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
	_shapeLogger->closeShape("clSnapclPoint");
	_shapeLogger->closeShape("clDeleteByLength");
	_shapeLogger->closeShape("clNoDeleteByLength");
	_shapeLogger->closeShape("CPSnappedOnCL");
	_shapeLogger->closeShape("CPFromInterWithCL");

	delete _mlsBorderSmoothed;	
}

///
///
///
void app::calcul::CFeatGenerationOp::GenerateConnectingLinesByCountry(
	std::string const& countryCodeDouble,
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
	double angleMaxToCutBorder = themeParameters->getValue(CFG_BOUNDARY_ANGLE_THRESHOLD).toDouble();
	angleMaxToCutBorder = angleMaxToCutBorder * M_PI / 180;

	//--
	ign::feature::FeatureIteratorPtr itBoundary = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsBoundary, ign::feature::FeatureFilter(countryCodeName + " = '" + _countryCodeDouble + "'"));
	while (itBoundary->hasNext())
	{
		ign::feature::Feature fBoundary = itBoundary->next();
		_logger->log(epg::log::INFO, "id boundary :"+ fBoundary.getId());

		std::string boundaryType = fBoundary.getAttribute("boundary_type").toString();
		// On ne traite que les frontières de type international_boundary ou coastline_sea_limit
		// On aurait pu filtrer en entrée mais le filtre semble très long, peut-être à cause de l'enum qui oblige à utiliser boundary_type::text like '%coastline_sea_limit%'
		if (boundaryType != "international_boundary" && boundaryType.find("coastline_sea_limit") == -1)
			continue;

		ign::geometry::LineString const& lsBoundary = fBoundary.getGeometry().asLineString();
		std::vector<ign::geometry::LineString> vLsBorderCutByAngle;
		_getBorderCutByAngle(lsBoundary, vLsBorderCutByAngle, angleMaxToCutBorder);
		for (size_t i = 0; i < vLsBorderCutByAngle.size(); ++i) {
			ign::geometry::LineString lsBoundaryCutByAngle = vLsBorderCutByAngle[i];

			ign::geometry::algorithm::BufferOpGeos buffOp;
			ign::geometry::GeometryPtr buffBorder(buffOp.buffer(lsBoundaryCutByAngle, distBuffer, 0, ign::geometry::algorithm::BufferOpGeos::CAP_FLAT));

			_getCLfromBorder(lsBoundaryCutByAngle, *buffBorder);
		}
	}

	_logger->log(epg::log::TITLE, "[ END CL GENERATION BY COUNTRY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::MergeConnectingLinesOnBorder(
	std::string const& countryCodeDouble,
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
void app::calcul::CFeatGenerationOp::SnapConnectingLines(
	std::string const& countryCodeDouble,
	bool verbose
) {
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

	_snapCl2Cl( distMaxClClosest );

	_logger->log(epg::log::TITLE, "[ END CL SNAPPING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::DeleteConnectingLines(
	std::string const& countryCodeDouble,
	bool verbose
) {
	CFeatGenerationOp op(countryCodeDouble, verbose);
    op._deleteConnectingLines();
}

///
///
///
void app::calcul::CFeatGenerationOp::_deleteConnectingLines() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN CL DELETING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	_deleteClByAngleAndDistEdges();

	_deleteCLUnderThreshold();

	_logger->log(epg::log::TITLE, "[ END CL DELETING FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::UpdateGeomConnectingLines(
	std::string const& countryCodeDouble,
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
	_logger->log(epg::log::TITLE, "[ BEGIN CL UPDATE GEOMETRY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	//--
	GraphType graphCL;
	_loadGraphCL(graphCL);

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const snapProjCl2edge = themeParameters->getValue(CL_SNAP_PROJ_CL_2_EDGE_DIST).toDouble();
	epg::Context* context = epg::ContextS::getInstance();
	std::string idName = context->getEpgParameters().getValue(ID).toString();
	std::string geomName = context->getEpgParameters().getValue(GEOM).toString();

	//--
	_updateGeomCL(snapProjCl2edge);

	//-- copie intermediaire pour debug
	epg::utils::CopyTableUtils::copyTable(
		themeParameters->getValue(CL_TABLE).toString(),
		idName,
		geomName,
		ign::geometry::Geometry::GeometryTypeLineString,
		themeParameters->getValue(CL_TABLE).toString()+"_no_cont",
		"",
		false,
		true
	);

	//--
	_setContinuityCl(graphCL);

	_logger->log(epg::log::TITLE, "[ END CL UPDATE GEOMETRY FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}



///
///
///
void app::calcul::CFeatGenerationOp::ComputeCL(
	std::string const& countryCodeDouble,
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
	std::string const& countryCodeDouble,
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

	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string boundaryTypeName = context->getEpgParameters().getValue(BOUNDARY_TYPE).toString();
	std::string typeInternationalBoundary = context->getEpgParameters().getValue(TYPE_INTERNATIONAL_BOUNDARY).toString();
	std::string typeCoastline = context->getEpgParameters().getValue(TYPE_COASTLINE).toString();

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distCLIntersected = themeParameters->getValue(CP_INTERSECTED_CL_DIST).toDouble();
	double const distUnderShoot = themeParameters->getValue(CP_UNDERSHOOT_DIST).toDouble();
	double const distCp2snapCl = themeParameters->getValue(CP_CP_2_CL_SNAP_DIST).toDouble();
	double const snapDistOnVertexFromCl = themeParameters->getValue(CP_VERTEX_CL_SNAP_DIST).toDouble();

	//--
    ign::feature::FeatureIteratorPtr itBoundary = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsBoundary, ign::feature::FeatureFilter(countryCodeName + " = '" + _countryCodeDouble + "'"));
	while (itBoundary->hasNext()) {

		ign::feature::Feature fBoundary = itBoundary->next();
		std::string const& boundaryType = fBoundary.getAttribute(boundaryTypeName).toString();
		ign::geometry::LineString const& lsBoundary = fBoundary.getGeometry().asLineString();

		_logger->log(epg::log::INFO, "id boundary :" + fBoundary.getId());
		
		if (boundaryType != typeInternationalBoundary && boundaryType.find(typeCoastline) == -1)
			continue;

		//--
		_getCPfromIntersectBorder(lsBoundary, distCLIntersected);

		//--
		ign::geometry::algorithm::BufferOpGeos buffOp;
		ign::geometry::GeometryPtr buffBorder(buffOp.buffer(lsBoundary, distUnderShoot, 0, ign::geometry::algorithm::BufferOpGeos::CAP_FLAT));
		_addToUndershootNearBorder(lsBoundary, *buffBorder, distUnderShoot);
	}

	//--
	_snapCPNearBy(0);

	//--
	std::map<std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint>> mClSplittedByCp;
	_snapCpOnClNearBy(distCp2snapCl, snapDistOnVertexFromCl, mClSplittedByCp);

	_cutClByCp(mClSplittedByCp);

	_logger->log(epg::log::TITLE, "[ END CP GENERATION FOR " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}
	
///
///
///
void app::calcul::CFeatGenerationOp::_init(
	std::string const& countryCodeDouble,
	bool verbose
) {
	//--
	_logger = epg::log::EpgLoggerS::getInstance();
	_logger->log(epg::log::INFO, "[START] initialization: " + epg::tools::TimeTools::getTime());

	//--
	epg::Context* context = epg::ContextS::getInstance();

	//--
	std::string const boundaryTableName = context->getEpgParameters().getValue(TARGET_BOUNDARY_TABLE).toString();

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const edgeTableName = context->getEpgParameters().getValue(EDGE_TABLE).toString();
	std::string const cpTableName = themeParameters->getValue(CP_TABLE).toString();
	std::string const clTableName = themeParameters->getValue(CL_TABLE).toString();
	std::string const listValueFormwayBigDist2merge = themeParameters->getValue(CP_VALUE_FORMWAY_BIGDIST2MERGE).toString();
	std::string smoothedBoundaryTableName = themeParameters->getValue(BOUNDARY_SMOOTHED_TABLE).toString();
	if (smoothedBoundaryTableName == "")
		smoothedBoundaryTableName = boundaryTableName;
	
	//--
	_shapeLogger = epg::log::ShapeLoggerS::getInstance();
	_shapeLogger->addShape("CLBeforeMerge", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClMergedBeforeUpdate", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDeletedNoCandidatefound", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDoublon", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDeleteByAngleDistEdges", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDebug", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("edgeClCutByCp", epg::log::ShapeLogger::LINESTRING); 
	_shapeLogger->addShape("lsBorderCutByAngle", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("clSnapclPoint", epg::log::ShapeLogger::POINT);
	_shapeLogger->addShape("clDeleteByLength", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("clNoDeleteByLength", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("CPSnappedOnCL", epg::log::ShapeLogger::POINT);
	_shapeLogger->addShape("CPFromInterWithCL", epg::log::ShapeLogger::POINT);

	//--
	_countryCodeDouble = countryCodeDouble;
	epg::tools::StringTools::Split(_countryCodeDouble, "#", _vCountriesCodeName);
	_verbose = verbose;

	///recuperation de la liste des attributs à concatener, de w et json dans la fusion des attributs
	std::string listAttrWName = themeParameters->getValue(LIST_ATTR_W).toString();
	std::string listAttrJsonName = themeParameters->getValue(LIST_ATTR_JSON).toString();
	_attrMergerOnBorder.setLists( listAttrWName, listAttrJsonName, "/");
	
	//--
	std::vector<std::string> vValueFormwayBigDist2merge;
	epg::tools::StringTools::Split(listValueFormwayBigDist2merge, "/", vValueFormwayBigDist2merge);
	for (size_t i = 0; i < vValueFormwayBigDist2merge.size(); ++i) {
		_sFormwayValues4BigDist2Merge.insert(vValueFormwayBigDist2merge[i]);
	}

	//--
	_reqFilterEdges2generateCF = themeParameters->getValue(SQL_FILTER_EDGES_2_GENERATE_CF).toString();

	//--
	_fsBoundary = context->getDataBaseManager().getFeatureStore(boundaryTableName, idName, geomName);
	_fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);
	_fsCP = context->getDataBaseManager().getFeatureStore(cpTableName, idName, geomName);
	_fsCL = context->getDataBaseManager().getFeatureStore(clTableName, idName, geomName);

	// id generator
	_idGeneratorCP = epg::sql::tools::IdGeneratorInterfacePtr(epg::sql::tools::IdGeneratorFactory::getNew(*_fsCP, "CONNECTINGPOINT"));
	_idGeneratorCL = epg::sql::tools::IdGeneratorInterfacePtr(epg::sql::tools::IdGeneratorFactory::getNew(*_fsCL, "CONNECTINGLINE"));

	//--
	ign::feature::FeatureFilter filter(countryCodeName + " = '" + _countryCodeDouble + "'");
	_mlsBorderSmoothed = new epg::tools::MultiLineStringTool(filter, *context->getDataBaseManager().getFeatureStore(smoothedBoundaryTableName, idName, geomName));
	
	//--
	_logger->log(epg::log::TITLE, "[ END INITIALIZATION ] : " + epg::tools::TimeTools::getTime());
}

///
///
///
void app::calcul::CFeatGenerationOp::_getBorderCutByAngle(
	ign::geometry::LineString const& lsBorder,
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
	ign::geometry::LineString const& lsBorder,
	ign::geometry::Geometry const& buffBorder
) const {
	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distBuffer = themeParameters->getValue(CL_BUFFER_DIST).toDouble();
	double const thresholdNoCL = themeParameters->getValue(CL_THRESHOLD_NO_CL).toDouble();
	double const ratioInBuff = themeParameters->getValue(CL_RATIO_IN_BUFFER).toDouble();
	double const snapOnVertexBorder = themeParameters->getValue(CL_SNAP_ON_VERTEX_BORDER_DIST).toDouble();
	double angleMax = themeParameters->getValue(CL_BORDER_MAX_ANGLE).toDouble();
	angleMax = angleMax * M_PI / 180;

	std::vector<ign::feature::FeatureAttributeType> listAttrEdge = _fsEdge->getFeatureType().attributes();

	ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + buffBorder.toString() + "'),3035))");
	epg::tools::FilterTools::addAndConditions(filter, countryCodeName +" NOT LIKE '%#%'");
	//ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_GeomFromText('" + buffBorder->toString() + "'))");
	if (_reqFilterEdges2generateCF != "")
		epg::tools::FilterTools::addAndConditions(filter, _reqFilterEdges2generateCF);

	ign::feature::FeatureIteratorPtr eit = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filter);
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filter);
	if (numFeatures == 0)
		return;
	boost::progress_display display(numFeatures, std::cout, "[ creating connecting lines % complete ]\n");

	while (eit->hasNext())
	{
		++display;
		//
		ign::feature::Feature fEdge = eit->next();
		ign::geometry::LineString const& lsEdge = fEdge.getGeometry().asLineString();
	
		//ign::geometry::Geometry* geomIntersect = lsEdge.Intersection(*buff);
		std::vector<ign::geometry::LineString> vLsProjectedOnBorder;
		app::geometry::tools::LineStringSplitter lsSplitter(lsEdge);
		lsSplitter.addCuttingGeometry(buffBorder);
		std::vector<ign::geometry::LineString> subEdgesBorder = lsSplitter.getSubLineStringsZ();

		//pas d'intersection par le buffer
		if (subEdgesBorder.size() == 1) {			
			double angleEdgBorder = _getAngleEdgeWithBorder(lsEdge, lsBorder);		
			//si l'edge est "proche" on considere qu'il est entierement dans le buffer et longe la frontiere
			if (lsEdge.distance(lsBorder) < distBuffer && (angleEdgBorder < angleMax || angleEdgBorder > (M_PI - angleMax) ) ) {
				ign::geometry::LineString lsCL = _getGeomCL(*_mlsBorderSmoothed, lsEdge, distBuffer, snapOnVertexBorder);
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
				ign::geometry::LineString const& lsSubEdgeCurr = subEdgesBorder[i];
				double currentSubEdgeLength = lsSubEdgeCurr.length();

 				double angleSubEdgBorder = _getAngleEdgeWithBorder(lsSubEdgeCurr, lsBorder);

				int numSeg = static_cast<int>(std::floor(lsSubEdgeCurr.numSegments() / 2.));
				ign::geometry::Point interiorPointSEC = epg::tools::geometry::interpolate(lsSubEdgeCurr, numSeg, 0.5);
				bool isSubSegInBuff = false;

				if (buffBorder.contains(interiorPointSEC) && (angleSubEdgBorder < angleMax || angleSubEdgBorder > (M_PI - angleMax) ) ) {
					isSubSegInBuff = true;
					numlastSubInBuff = i;

					lengthInBuff += currentSubEdgeLength;
					if (numfirstSubInBuff < 0)
						numfirstSubInBuff = i;
				}

				if (isSubSegInBuff || currentSubEdgeLength <= thresholdNoCL)
					lengthNearByBuff += currentSubEdgeLength;

				if ((currentSubEdgeLength > thresholdNoCL && !isSubSegInBuff) || i == subEdgesBorder.size() - 1) {
					if (lengthInBuff > ratioInBuff * lengthNearByBuff) {
						//recup ptStart, ptFin et proj des pt sur la border
						//recup de la border entre ces points pour recup de la geom CL
						ign::geometry::LineString lsRef = _getLinestring(subEdgesBorder, numfirstSubInBuff, numlastSubInBuff);
						ign::geometry::LineString lsCL = _getGeomCL(*_mlsBorderSmoothed, lsRef, distBuffer, snapOnVertexBorder);
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
			featCL.setAttribute(linkedFeatIdName, ign::data::String(fEdge.getId()));
			for (std::vector<ign::feature::FeatureAttributeType>::iterator vit = listAttrEdge.begin(); vit != listAttrEdge.end(); ++vit) {
				std::string attrName = vit->getName();
				if (attrName == geomName || attrName == idName || !_fsCL->getFeatureType().hasAttribute(attrName) )
					continue;
				featCL.setAttribute(attrName,fEdge.getAttribute(attrName));
			}

			std::string idCL = _idGeneratorCL->next();
			_fsCL->createFeature(featCL, idCL);
			
			//--
			_shapeLogger->writeFeature("CLBeforeMerge", featCL);
		}

	}	
}

///
///
///
ign::geometry::LineString app::calcul::CFeatGenerationOp::_getLinestring(
	std::vector<ign::geometry::LineString> const& subEdges,
	size_t idFirst,
	size_t idLast
) const {
	ign::geometry::LineString lsResult = subEdges[idFirst];
	for ( size_t i = idFirst + 1 ; i < idLast + 1 ; ++i ) {
		for ( size_t j = 1 ; j < subEdges[i].numPoints() ; ++j ) {
			lsResult.addPoint(subEdges[i].pointN(j));
		}
	}
	return lsResult;
}

///
///
///
void app::calcul::CFeatGenerationOp::_addToUndershootNearBorder(
	ign::geometry::LineString const& lsBorder,
	ign::geometry::Geometry const& buffBorder,
	double distUnderShoot
) const {
	//--
	epg::Context* context = epg::ContextS::getInstance();

	//--
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	//--
	std::vector<ign::feature::FeatureAttributeType> listAttrEdge = _fsEdge->getFeatureType().attributes();

	ign::feature::FeatureFilter filterBuffBorder("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + buffBorder.toString() + "'),3035))");
	if (_reqFilterEdges2generateCF != "")
		epg::tools::FilterTools::addAndConditions(filterBuffBorder, _reqFilterEdges2generateCF);
	//on ne prend que les edges ayant un cc simple pour ne pas créer de CP là où il y a des CLs
	epg::tools::FilterTools::addAndConditions(filterBuffBorder, "(" + countryCodeName + " = '" + _vCountriesCodeName[0] + "' or " + countryCodeName + " = '" + _vCountriesCodeName[1] + "')");

	ign::feature::FeatureIteratorPtr eit = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterBuffBorder);
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterBuffBorder);
	boost::progress_display display(numFeatures, std::cout, "[ computing connecting point by resolving undershoot % complete ]\n");

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
		ign::feature::FeatureIteratorPtr eitArroundPt = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterArroundPt);
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
		double distMin = std::numeric_limits<double>::infinity();
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
		fCF.setAttribute(linkedFeatIdName, ign::data::String(fEdge.getId()));
		for (std::vector<ign::feature::FeatureAttributeType>::iterator vit = listAttrEdge.begin();
			vit != listAttrEdge.end(); ++vit) {
			std::string attrName = vit->getName();
			if (attrName == geomName || attrName == idName || !_fsCP->getFeatureType().hasAttribute(attrName))
				continue;
			fCF.setAttribute(attrName, fEdge.getAttribute(attrName));
		}

		std::string idCP = _idGeneratorCP->next();
		_fsCP->createFeature(fCF, idCP);
		
		//--
		_shapeLogger->writeFeature("CPUndershoot", fCF);
	}
}

///
///
///
void app::calcul::CFeatGenerationOp::_getCPfromIntersectBorder(
	ign::geometry::LineString const& lsBorder,
	double distCLIntersected
) const {
	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	std::vector<ign::feature::FeatureAttributeType> const& listAttrEdge = _fsEdge->getFeatureType().attributes();

	//on prend les edges qui intersectent les frontières et ou ceux qui intersectent les Cls
	ign::feature::FeatureFilter filterFeaturesToMatch("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + lsBorder.toString() + "'),3035))");
	epg::tools::FilterTools::addOrConditions(filterFeaturesToMatch,"ST_INTERSECTS(" + geomName + ", (SELECT ST_Union(array(SELECT "+geomName+" FROM "+ _fsEdge->getTableName()+" WHERE "+ countryCodeName + " = '" + _countryCodeDouble + "'))))");

	if (_reqFilterEdges2generateCF != "")
		epg::tools::FilterTools::addAndConditions(filterFeaturesToMatch, _reqFilterEdges2generateCF);
	//on ne prend que les edges ayant un cc simple pour ne pas créer de CP là où il y a des CLs
	epg::tools::FilterTools::addAndConditions(filterFeaturesToMatch,"("+ countryCodeName+" = '" +_vCountriesCodeName[0]+"' OR "+countryCodeName + " = '" + _vCountriesCodeName[1] +"')");
	ign::feature::FeatureIteratorPtr itFeaturesToMatch = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterFeaturesToMatch);

	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterFeaturesToMatch);
	boost::progress_display display(numFeatures, std::cout, "[ creating connecting points % complete ]\n");

	while (itFeaturesToMatch->hasNext())
	{
		++display;

		ign::feature::Feature fToMatch = itFeaturesToMatch->next();
		ign::geometry::LineString const& lsFToMatch = fToMatch.getGeometry().asLineString();

		//DEBUG
		if( fToMatch.getId() == "8aee4cd4-6e72-416e-8187-cd8e8e6a1c93" ) {
			if ( lsBorder.distance(ign::geometry::Point(3901953.97, 3142259.44)) < 1) {
				bool test = true;
			}
		}

		ign::geometry::GeometryPtr geomPtr(lsFToMatch.Intersection(lsBorder));

		ign::feature::Feature fClArround;
		bool hasClConnected = _isEdgeConnected2cl(lsFToMatch, lsFToMatch.getEnvelope().expandBy(distCLIntersected), fClArround, distCLIntersected);
		
		bool debugInterWithCl = false;
		if (hasClConnected) {
			ign::geometry::GeometryPtr geomIntersectCl(fClArround.getGeometry().Intersection(lsFToMatch));

			//si il existe une intersection entre l'edge et une CL, et qu'elle est sous un seuil, on prend cette intersection à la place de celle avec la frontière
			if (!geomIntersectCl->isNull() && !geomIntersectCl->isEmpty()) {
				double distance = geomIntersectCl->distance(lsBorder);

				if(distance >= 0 && distance < distCLIntersected) {
					if (geomPtr->isNull() || geomPtr->distance(*geomIntersectCl) < distCLIntersected){ //utile ce test ou le precedent suffit?
						geomPtr.reset(geomIntersectCl.release());
						debugInterWithCl = true;
					}
				}
			}
		}

		if (geomPtr->isNull() || geomPtr->isEmpty())
			continue;

		//--
		ign::feature::Feature fCF = _fsCP->newFeature();
		fCF.setAttribute(linkedFeatIdName, ign::data::String(fToMatch.getId()));
		for (std::vector<ign::feature::FeatureAttributeType>::const_iterator vit = listAttrEdge.begin();
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

			//DEBUG
			if (debugInterWithCl)
				_shapeLogger->writeFeature("CPFromInterWithCL", fCF);
		}
		else if (geomPtr->isGeometryCollection())
		{
			ign::geometry::GeometryCollection geomCollect = geomPtr->asGeometryCollection();
			for (size_t i = 0; i < geomCollect.numGeometries(); ++i)
			{
				if (geomCollect.geometryN(i).isPoint())
				{
					ign::geometry::Point const& ptIntersect = geomCollect.geometryN(i).asPoint();		
					fCF.setGeometry(ptIntersect);
					std::string idCP = _idGeneratorCP->next();
					_fsCP->createFeature(fCF, idCP);

					//DEBUG
					if (debugInterWithCl)
						_shapeLogger->writeFeature("CPFromInterWithCL", fCF);
				}
			}
		}
	}
}

///
///
///
double app::calcul::CFeatGenerationOp::_getAngleEdgeWithBorder(
	ign::geometry::LineString const& lsEdge,
	ign::geometry::LineString const& lsBorder
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
ign::geometry::LineString app::calcul::CFeatGenerationOp::_getGeomCL(
	epg::tools::MultiLineStringTool & mslBorder,
	ign::geometry::LineString const& lsToProject,
	double distMaxBorder,
	double snapOnVertexBorder
) const {

	ign::geometry::Point ptStartToProject = lsToProject.startPoint();
	ign::geometry::Point ptEndToProject = lsToProject.endPoint();
	// std::pair< bool, ign::geometry::LineString > pathFound = mslBorder.getPathAlong(ptStartToProject, ptEndToProject, lsStart2EndToPrject, 2* distMaxBorder, distMaxBorder+1);
	std::pair< bool, ign::geometry::LineString > pathFound = mslBorder.getPathAlong(ptStartToProject, ptEndToProject, lsToProject, 1000, 1000, snapOnVertexBorder);

	if (pathFound.first)
		return pathFound.second;

	return ign::geometry::LineString();
}

///
///
///
void app::calcul::CFeatGenerationOp::_snapCPNearBy(
	double snapOnVertexBorder
) const {
	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double const distMergeCP = themeParameters->getValue(CP_MERGE_DIST_CP).toDouble();
	double const distMergeTractorCP = themeParameters->getValue(CP_MERGE_DIST_TRACTOR_CP).toDouble();
	double maxDistMerge = std::max(distMergeCP, distMergeTractorCP);

	ign::feature::FeatureFilter filterCP;
	for (size_t i = 0; i < _vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCP, countryCodeName + " = '" + _vCountriesCodeName[i] + "'");
	}
	ign::feature::FeatureIteratorPtr itCP = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCP, filterCP);

	std::set<std::string> sCP2Snap;
	std::string separator = "#";
	while (itCP->hasNext())
	{
		ign::feature::Feature fCPCurr = itCP->next();
		std::string idCP = fCPCurr.getId();

		//DEBUG
		_logger->log(epg::log::DEBUG, idCP);

		if (sCP2Snap.find(idCP) != sCP2Snap.end())
			continue;

		// on recupere recursivement les CP proches
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
			ign::feature::FeatureIteratorPtr fitBorder = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsBoundary, filterBorderNearCP);
			while (fitBorder->hasNext()) {
				ign::feature::Feature fBorder = fitBorder->next();
				ign::geometry::LineString const& lsBorder = fBorder.getGeometry().asLineString();

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
	std::string const typeName = themeParameters->getValue(FORM_OF_WAY_NAME).toString();
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
	ign::geometry::Geometry const& geomObjNearCl,
	ign::geometry::Envelope const& envArroundGeom,
	ign::feature::Feature & fCl2SnapOn,
	double distMinCl
) const {
	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	
	ign::feature::FeatureFilter filterArroundCp(countryCodeName + " = '" + _countryCodeDouble + "'");
	filterArroundCp.setExtent(envArroundGeom);

	ign::feature::FeatureIteratorPtr itClArround = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterArroundCp);
	bool hasEdgeConnected2Cl = false;
	while (itClArround->hasNext())
	{
		ign::feature::Feature fClArround = itClArround->next();
		ign::geometry::LineString const& lsClArround = fClArround.getGeometry().asLineString();

		double dist = lsClArround.distance(geomObjNearCl);

		if (dist < distMinCl) {
			hasEdgeConnected2Cl = true;
			distMinCl = dist;
			fCl2SnapOn = fClArround;
		}
	}
	return hasEdgeConnected2Cl;
}

///
///
///
void app::calcul::CFeatGenerationOp::_snapCpOnClNearBy(
	double distCp2snapCl,
	double snapDistOnVertexFromCl,
	std::map<std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint>> & mClSplittedByCp
) const {
	//--
	epg::Context* context = epg::ContextS::getInstance();

	//--
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	
	//--
	ign::feature::FeatureFilter filterCp;
	for (size_t i = 0; i < _vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCp, countryCodeName + " = '" + _vCountriesCodeName[i] + "'");
	}
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCP, filterCp);
	boost::progress_display display(numFeatures, std::cout, "[ snapping cp on cl % complete ]\n");

	//--
	std::map<std::string, ign::feature::Feature> mSnappedCpOnCl;
	ign::feature::FeatureIteratorPtr itCp = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCP, filterCp);
	while (itCp->hasNext())
	{
		++display;
		ign::feature::Feature fCp = itCp->next();
		ign::geometry::Point const& ptCp = fCp.getGeometry().asPoint();

		ign::feature::Feature fCl2SnapOn;
		double distMinCl = distCp2snapCl * 2;
		bool hasClNearBy = _isEdgeConnected2cl(ptCp, ptCp.getEnvelope().expandBy(distCp2snapCl), fCl2SnapOn, distMinCl);
			
		//pas de Cl proche pour snapper
		if (!hasClNearBy)
			continue;

		//on projette le Cp sur la Cl identifiée
		ign::geometry::LineString const& lsCl2SnapOn = fCl2SnapOn.getGeometry().asLineString();
		ign::geometry::Point ptSnap = epg::tools::geometry::project(lsCl2SnapOn, ptCp, snapDistOnVertexFromCl);
		ptSnap.setZ(0);

		//on modifie la geom du Cp
		fCp.setGeometry(ptSnap);
		mSnappedCpOnCl[fCp.getId()] = fCp;

		//on stocke la geom du Cp comme coupure de la Cl
		std::string idCl = fCl2SnapOn.getId();
		if (mClSplittedByCp.find(idCl) == mClSplittedByCp.end()) { //pas de decoupe pour cette Cl
			ign::geometry::MultiPoint mpt;
			mpt.addGeometry(ptSnap);
			std::pair<ign::feature::Feature, ign::geometry::MultiPoint> pairClSnappedOn = std::make_pair(fCl2SnapOn, mpt);
			mClSplittedByCp[idCl] = pairClSnappedOn;
		}
		else { //il y a déjà des points de coupure sur cette CL
			std::pair<ign::feature::Feature, ign::geometry::MultiPoint> pairClSnappedOn = mClSplittedByCp.find(idCl)->second;
			pairClSnappedOn.second.addGeometry(ptSnap);
			mClSplittedByCp[idCl] = pairClSnappedOn;
		}
	}

	for (std::map<std::string, ign::feature::Feature>::iterator mit = mSnappedCpOnCl.begin(); mit != mSnappedCpOnCl.end(); ++mit) {
		_fsCP->modifyFeature(mit->second);

		//DEBUG
		_shapeLogger->writeFeature("CPSnappedOnCL", mit->second);
	}
}

///
///
///
void app::calcul::CFeatGenerationOp::_cutClByCp(
	std::map<std::string,std::pair<ign::feature::Feature, ign::geometry::MultiPoint>> const& mClSplittedByCp
) const {
	std::set<std::string> sCl2delete;
	boost::progress_display display(mClSplittedByCp.size(), std::cout, "[ cutting cl by cp % complete ]\n");
	for (std::map<std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint>>::const_iterator mit = mClSplittedByCp.begin();
		mit != mClSplittedByCp.end(); ++mit) {
		++display;
		ign::feature::Feature const& fEdgDoubl2Cut = mit->second.first;
		ign::geometry::LineString const& lsCl = fEdgDoubl2Cut.getGeometry().asLineString();
		app::geometry::tools::LineStringSplitter lsSplitterClSnappedOn(lsCl, 0.1); //on snappe à 10cm si il y a plusieurs coupures
		//boucle sur mpt
		ign::geometry::MultiPoint const& mpt = mit->second.second;
		for (size_t i = 0; i < mpt.numGeometries(); ++i)
			lsSplitterClSnappedOn.addCuttingGeometry(mpt.geometryN(i).asPoint());

		std::vector<ign::geometry::LineString> subCl2cut = lsSplitterClSnappedOn.getSubLineStringsZ();

		if (subCl2cut.size() == 1)
			continue;
			
		sCl2delete.insert(fEdgDoubl2Cut.getId());

		for (size_t i = 0; i < subCl2cut.size(); ++i) {
			ign::feature::Feature fNewEdgeDouble = fEdgDoubl2Cut;
			fNewEdgeDouble.setGeometry(subCl2cut[i]);
			_fsEdge->createFeature(fNewEdgeDouble);

			//--
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
	std::map<std::string, ign::feature::Feature> & mCPNear
) const {
	//--
	epg::Context* context = epg::ContextS::getInstance();

	//--
	std::string const idName = context->getEpgParameters().getValue(ID).toString();

	//--
	mCPNear[fCP.getId()] = fCP;
	
	ign::feature::FeatureFilter filterArroundCP;
	for (std::map<std::string, ign::feature::Feature>::iterator mit = mCPNear.begin(); mit != mCPNear.end(); ++mit) {
		epg::tools::FilterTools::addAndConditions(filterArroundCP, idName + " <> '" + mit->first + "'");	//(idName + " <> '" + fCP.getId() + "'");
	}

	filterArroundCP.setExtent(fCP.getGeometry().getEnvelope().expandBy(distMergeCP));
	ign::feature::FeatureIteratorPtr itArroundCP = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCP, filterArroundCP);

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
	
	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	
	//--
	GraphType graphCl;
	_loadGraphCL(graphCl);

	//--
	ign::feature::FeatureFilter filterCL(countryCodeName + " = '" + _countryCodeDouble + "'");
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCL);
	boost::progress_display displayLoad(numFeatures, std::cout, "[ snapping cl to cl % complete ]\n");

	ign::feature::FeatureIteratorPtr itCL = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCL);
	std::map<std::string, std::pair<ENDING, ign::feature::Feature>> mFClModified;
	while (itCL->hasNext()) {
		++displayLoad;

		ign::feature::Feature fCL = itCL->next();
		ign::geometry::LineString lsCl = fCL.getGeometry().asLineString();

		edge_descriptor edCl = graphCl.getInducedEdges(fCL.getId()).second[0].descriptor;

		if(graphCl.degree(graphCl.source(edCl)) == 1)
			_snapTo(distMaxClClosest, START, fCL, lsCl, mFClModified);
			
		if(graphCl.degree(graphCl.target(edCl)) == 1)
			_snapTo(distMaxClClosest, END, fCL, lsCl, mFClModified);
	}

	for (std::map<std::string, std::pair<ENDING, ign::feature::Feature>>::iterator mit = mFClModified.begin(); mit != mFClModified.end(); ++mit) {
		_fsCL->modifyFeature(mit->second.second);
	}

}

///
///
///
void app::calcul::CFeatGenerationOp::_snapTo( 
	double distMaxClClosest,
	app::calcul::CFeatGenerationOp::ENDING ending,
	ign::feature::Feature & fCL,
	ign::geometry::LineString & newClGeom,
	std::map<std::string, std::pair<app::calcul::CFeatGenerationOp::ENDING, ign::feature::Feature>> & mFClModified
) const {
	std::vector<std::pair<ENDING, ign::feature::Feature>> vClExtremityClose = _getClExtremityClose( distMaxClClosest, ending, fCL);

	ign::geometry::LineString const& clGeom = fCL.getGeometry().asLineString();
	ign::geometry::Point const& endingGeom = ending == START ? clGeom.startPoint() : clGeom.endPoint();
	ign::geometry::Point const& otherEndingGeom = ending == START ? clGeom.endPoint() : clGeom.startPoint();

	ign::geometry::MultiPoint mp;
	mp.addGeometry(ending == START ? clGeom.startPoint() : clGeom.endPoint() );
	
	for ( int i = vClExtremityClose.size()-1 ; i >= 0 ; --i ) {
		
		std::map<std::string, std::pair<ENDING, ign::feature::Feature>>::const_iterator mit = mFClModified.find(vClExtremityClose[i].second.getId());
		if ( mit != mFClModified.end() )
			if ( mit->second.first == BOTH || mit->second.first == vClExtremityClose[i].first ) {
				//modif deja faite
				vClExtremityClose.erase(vClExtremityClose.begin() + i);
				continue;
			}
		// gestion des cl < distMaxClClosest = snapper l'extremité la plus proche 
		if ( clGeom.length() < distMaxClClosest ) {
			ign::geometry::LineString const& ClExtremityCloseGeom = vClExtremityClose[i].second.getGeometry().asLineString();
			ign::geometry::Point const& extremityGeom = vClExtremityClose[i].first == START ? ClExtremityCloseGeom.startPoint() : ClExtremityCloseGeom.endPoint();
			if ( otherEndingGeom.distance(extremityGeom) < endingGeom.distance(extremityGeom) ) {
				vClExtremityClose.erase(vClExtremityClose.begin() + i);
				continue;
			}
		}

		ign::geometry::LineString const& clCloseGem = vClExtremityClose[i].second.getGeometry().asLineString();

		if ( vClExtremityClose[i].first != END ) {
			mp.addGeometry(clCloseGem.startPoint());
		}
		if ( vClExtremityClose[i].first != START ) {
			mp.addGeometry(clCloseGem.endPoint());
		}
	}

	if (!vClExtremityClose.empty()) {
		ign::geometry::Point newPt = mp.asMultiPoint().getCentroid();
		newPt.setZ(0);
		
		newClGeom.setPointN(newPt, ending == START ? 0 : newClGeom.numPoints()-1);
		fCL.setGeometry(newClGeom);

		ENDING ending_ = ending;
		std::map<std::string, std::pair<app::calcul::CFeatGenerationOp::ENDING, ign::feature::Feature>>::iterator mit = mFClModified.find(fCL.getId());
		if ( mit != mFClModified.end() && mit->second.first != ending )
			ending_ = BOTH;

		mFClModified[fCL.getId()] = std::make_pair(ending_, fCL);

		for ( size_t i = 0 ; i < vClExtremityClose.size() ; ++i ) {
			ign::geometry::LineString clGeom = vClExtremityClose[i].second.getGeometry().asLineString();
			if ( vClExtremityClose[i].first != END )
				clGeom.setPointN(newPt, 0);
			if ( vClExtremityClose[i].first != START )
				clGeom.setPointN(newPt, clGeom.numPoints()-1);
			vClExtremityClose[i].second.setGeometry(clGeom);

			std::map<std::string, std::pair<ENDING, ign::feature::Feature>>::iterator mit = mFClModified.find(vClExtremityClose[i].second.getId());
			if ( mit != mFClModified.end() ) {
				mit->second.first = BOTH;
				mit->second.second = vClExtremityClose[i].second;
			} else {
				mFClModified.insert(std::make_pair(vClExtremityClose[i].second.getId(), vClExtremityClose[i]));
			}
		}

		ign::feature::Feature featConnPoint;
		featConnPoint.setGeometry(newPt);
		_shapeLogger->writeFeature("clSnapclPoint", featConnPoint);
	}
}

///
///
///
std::vector<std::pair<app::calcul::CFeatGenerationOp::ENDING, ign::feature::Feature>> app::calcul::CFeatGenerationOp::_getClExtremityClose(
	double distMaxClClosest,
	app::calcul::CFeatGenerationOp::ENDING ending,
	ign::feature::Feature const& fClCurr
) const {
	//recup des CL proches si dist = 0 rien faire, si dist inf distMaxClClosest snapper

	ign::geometry::LineString const& clGeom = fClCurr.getGeometry().asLineString();

	ign::geometry::Point const& ptClCurr = ending == START ? clGeom.startPoint() : clGeom.endPoint();

	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	//--
	ign::feature::FeatureFilter filterArroundCl (idName + " <> '" + fClCurr.getId() + "'");
	epg::tools::FilterTools::addAndConditions(filterArroundCl, countryCodeName + " = '" + _countryCodeDouble + "'");
	filterArroundCl.setExtent(ptClCurr.getEnvelope().expandBy(distMaxClClosest));

	ign::feature::FeatureIteratorPtr itClArround = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterArroundCl);
	std::vector<std::pair<ENDING, ign::feature::Feature>> vConnectedCl;
	while (itClArround->hasNext()) {

		ign::feature::Feature fClArround = itClArround->next();
		ign::geometry::LineString const& lsClArround = fClArround.getGeometry().asLineString();

		std::map<double, ENDING> mEnding;
		mEnding.insert(std::make_pair(lsClArround.startPoint().distance(ptClCurr), START));
		mEnding.insert(std::make_pair(lsClArround.endPoint().distance(ptClCurr), END));

		if ( mEnding.begin()->first == 0 || mEnding.begin()->first > distMaxClClosest )
			continue;
		
		ENDING endingClAround = NONE;

		if ( mEnding.begin()->first < distMaxClClosest )
			endingClAround = mEnding.begin()->second;

		if ( mEnding.rbegin()->first < distMaxClClosest && clGeom.length() > distMaxClClosest )
			endingClAround = endingClAround != NONE ? BOTH : mEnding.rbegin()->second;
		
		if ( endingClAround != NONE )
			vConnectedCl.push_back(std::make_pair(endingClAround, fClArround));
	}
	return vConnectedCl;
}

///
///
///
void app::calcul::CFeatGenerationOp::_mergeIntersectingClWithGraph(
	double distMaxEdges,
	double snapProjCl2edge
) const {
	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	std::string separator = "#";

	//--
	GraphType graphCl;
	ign::geometry::graph::tools::SnapRoundPlanarizer< GraphType >  planarizerCl(graphCl, 1e4);//scale =1e2 -> precision de 0.01

	//--
	ign::feature::FeatureFilter filterCL;
	for (size_t i = 0; i < _vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCL, countryCodeName + " = '" + _vCountriesCodeName[i] + "'");
	}

	//--
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCL);
	boost::progress_display displayLoad(numFeatures, std::cout, "[ loading cl planar graph % complete ]\n");

	ign::feature::FeatureIteratorPtr itCL = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCL);
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
	boost::progress_display display(graphCl.numEdges(), std::cout, "[ merging connecting lines % complete ]\n");
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

			double distMin = std::numeric_limits<double>::infinity();
			std::pair<std::string, std::string> cl2merge;

			for (std::map<std::string, ign::feature::Feature>::iterator mit1 = mIdClOriginsCountry1.begin(); mit1 != mIdClOriginsCountry1.end(); ++mit1) {

				std::string idEdgeLinked1 = mit1->second.getAttribute(linkedFeatIdName).toString();
				if (sEdgesMerged.find(mit1->first) != sEdgesMerged.end())//deja utilise pour merge
					continue;
				ign::feature::Feature fEdge1;
				_fsEdge->getFeatureById(idEdgeLinked1, fEdge1);
				ign::geometry::LineString lsEdge1 = fEdge1.getGeometry().asLineString();
				ign::geometry::GeometryPtr geomClEdge1(_getGeomProjClOnEdge(lsCl, lsEdge1, snapProjCl2edge));

				for (std::map<std::string, ign::feature::Feature>::iterator mit2 = mIdClOriginsCountry2.begin(); mit2 != mIdClOriginsCountry2.end(); ++mit2) {

					std::string idEdgeLinked2 = mit2->second.getAttribute(linkedFeatIdName).toString();
					if (sEdgesMerged.find(mit2->first) != sEdgesMerged.end())//deja utilise pour merge
						continue;
					ign::feature::Feature fEdge2;
					_fsEdge->getFeatureById(idEdgeLinked2, fEdge2);
					ign::geometry::LineString lsEdge2 = fEdge2.getGeometry().asLineString();
					ign::geometry::GeometryPtr geomClEdge2(_getGeomProjClOnEdge(lsCl, lsEdge2, snapProjCl2edge));

					double hausdorffDistEdges = ign::geometry::algorithm::OptimizedHausdorffDistanceOp::distance(*geomClEdge1, *geomClEdge2);

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

			lsCl.setFillZ(0);
			fClNew.setGeometry(lsCl);

			std::string idCLNew = _idGeneratorCL->next();
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
	ign::feature::FeatureIteratorPtr itCL = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCL);
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCL);
	boost::progress_display display(numFeatures, std::cout, "[ connecting lines merging with # % complete ]\n");

	std::set<std::string> sCL2Merged;
	std::string separator = "#";

	while (itCL->hasNext())
	{
		++display;
		ign::feature::Feature fCLCurr = itCL->next();
		std::string idCLCurr = fCLCurr.getId();
		ign::geometry::LineString const& lsCurr = fCLCurr.getGeometry().asLineString();
		std::string const& countryCodeCLCurr = fCLCurr.getAttribute(countryCodeName).toString();

		if (countryCodeCLCurr.find("#") != std::string::npos)
			continue;
		if (sCL2Merged.find(idCLCurr) != sCL2Merged.end())
			continue;
		sCL2Merged.insert(idCLCurr);

		ign::geometry::LineString lsBorder = _getBorderFromEdge(lsCurr);
		epg::tools::MultiLineStringTool mslBorder(lsBorder);

		ign::feature::FeatureFilter filterArroundCL;
		filterArroundCL.setPropertyConditions(countryCodeName + " != '" + countryCodeCLCurr + "'");
		filterArroundCL.setExtent(lsCurr.getEnvelope());
		ign::feature::FeatureIteratorPtr itArroundCL = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterArroundCL);
		while (itArroundCL->hasNext())
		{
			ign::feature::Feature fCLArround = itArroundCL->next();
			std::string idClArround = fCLArround.getId();
			ign::geometry::LineString const& lsClArround = fCLArround.getGeometry().asLineString();
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

			ign::geometry::LineString segmentSE(lsCurr.startPoint(), lsClArround.endPoint());
			ign::geometry::LineString lsSE = _getGeomCL(mslBorder, segmentSE, distMaxFromBorder, snapOnVertexBorder);

			ign::geometry::LineString segmentSS(lsCurr.startPoint(), lsClArround.startPoint());
			ign::geometry::LineString lsSS = _getGeomCL(mslBorder, segmentSS, distMaxFromBorder, snapOnVertexBorder);

			ign::geometry::LineString segmentES(lsCurr.endPoint(), lsClArround.startPoint());
			ign::geometry::LineString lsES = _getGeomCL(mslBorder, segmentES, distMaxFromBorder, snapOnVertexBorder);

			ign::geometry::LineString segmentEE(lsCurr.endPoint(), lsClArround.endPoint());
			ign::geometry::LineString lsEE = _getGeomCL(mslBorder, segmentEE, distMaxFromBorder, snapOnVertexBorder);

			ign::geometry::LineString lsIntersectedCL;
			
			//lsIntersectedCL = lsCurr;
			double lengthMin = std::numeric_limits<double>::infinity();

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
			}
			else {
				fCLNew = fCLCurr;
				_attrMergerOnBorder.mergeFeatAttribute(fCLNew, fCLArround, separator);
			}

			lsIntersectedCL.setFillZ(0);
			fCLNew.setGeometry(lsIntersectedCL);
			fCLNew.setAttribute(themeParameters->getValue(CF_STATUS).toString(), ign::data::String("edge_matched"));

			std::string idCLNew = _idGeneratorCL->next();
			_fsCL->createFeature(fCLNew, idCLNew);
			
			//--
			_shapeLogger->writeFeature("ClMergedBeforeUpdate", fCLNew);
		}
	}

	for (std::set<std::string>::iterator sit = sCL2Merged.begin(); sit != sCL2Merged.end(); ++sit) {
		_fsCL->deleteFeature(*sit);
	}

}

///
///
///
void app::calcul::CFeatGenerationOp::_deleteClByAngleAndDistEdges() const {
	
	_logger->log(epg::log::TITLE, "[ BEGIN DELETE CL BY ANGLE EDGES FOR : " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + _countryCodeDouble + "'");

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double angleMax = themeParameters->getValue(CL_EDGE_MAX_ANGLE).toDouble();
	double const distMax = themeParameters->getValue(CL_EDGE_MAX_DIST).toDouble();
	double const snapProjCl2edge = themeParameters->getValue(CL_SNAP_PROJ_CL_2_EDGE_DIST).toDouble();
	angleMax = angleMax * M_PI / 180;

	int numCl2delete = -1;
	while (numCl2delete != 0){
		//--
		GraphType graphCl;
		_loadGraphCL(graphCl);

		//--
		size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCLCountryCode);
		boost::progress_display display(numFeatures, std::cout, "[ deleting cl by angle between linked edges % complete ]\n");

		//--
		std::set<std::string> sCl2delete;
		ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCLCountryCode);
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
				_logger->log(epg::log::WARN, "Suppression CL " + fCl.getId() + " not matching linked edge : " + idEdgLinked1);
				_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fCl);

				sCl2delete.insert(fCl.getId());
				continue;
			}
			if (fEdg2.getId().empty()) {
				_logger->log(epg::log::WARN, "Suppression CL " + fCl.getId() + " not matching linked edge : " + idEdgLinked2);
				_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fCl);

				sCl2delete.insert(fCl.getId());
				continue;
			}
			ign::geometry::LineString lsEdg1 = fEdg1.getGeometry().asLineString();
			ign::geometry::LineString lsEdg2 = fEdg2.getGeometry().asLineString();
			ign::geometry::LineString lsCl = fCl.getGeometry().asLineString();

			ign::geometry::GeometryPtr geomProjClEdg1(_getGeomProjClOnEdge(lsCl, lsEdg1, snapProjCl2edge));
			ign::geometry::GeometryPtr geomProjClEdg2(_getGeomProjClOnEdge(lsCl, lsEdg2, snapProjCl2edge));
			if (geomProjClEdg1->isPoint()) {
				_logger->log(epg::log::WARN, "Suppression CL " + fCl.getId() + " not projecting on matching linked edge : " + idEdgLinked1);
				_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fCl);
				sCl2delete.insert(fCl.getId());
				continue;
			}
			if (geomProjClEdg2->isPoint()) {
				_logger->log(epg::log::WARN, "Suppression CL " + fCl.getId() + " not projecting on matching linked edge : " + idEdgLinked2);
				_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fCl);
				sCl2delete.insert(fCl.getId());
				continue;
			}
			ign::geometry::LineString const& lsProjClEdg1 = geomProjClEdg1->asLineString();
			ign::geometry::LineString const& lsProjClEdg2 = geomProjClEdg2->asLineString();

			ign::math::Vec2d vec1(lsProjClEdg1.endPoint().x() - lsProjClEdg1.startPoint().x(), lsProjClEdg1.endPoint().y() - lsProjClEdg1.startPoint().y());
			ign::math::Vec2d vec2(lsProjClEdg2.endPoint().x() - lsProjClEdg2.startPoint().x(), lsProjClEdg2.endPoint().y() - lsProjClEdg2.startPoint().y());
			double angleEdgesLinked = epg::tools::geometry::angle(vec1, vec2);

			double hausdorffDist = ign::geometry::algorithm::OptimizedHausdorffDistanceOp::distance(lsProjClEdg1, lsProjClEdg2);

			//si angle trop important entre les deux edges on ne crée pas de Cl
			if (angleEdgesLinked > angleMax && angleEdgesLinked < (M_PI - angleMax) || hausdorffDist > distMax) {
				sCl2delete.insert(fCl.getId());
				
				//--
				_shapeLogger->writeFeature("ClDeleteByAngleDistEdges", fCl);
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
ign::geometry::LineString app::calcul::CFeatGenerationOp::_getBorderFromEdge(
	ign::geometry::LineString const& lsEdgeOnBorder
) const {
	ign::feature::FeatureFilter filter;//(countryCodeName + " = 'be#fr'")
	filter.setExtent(lsEdgeOnBorder.getEnvelope());
	ign::feature::FeatureIteratorPtr itBoundary = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsBoundary, filter);
	while (itBoundary->hasNext()) {
		ign::feature::Feature fBorder = itBoundary->next();
		ign::geometry::LineString const& lsBorderTemp = fBorder.getGeometry().asLineString();
		double distBorder = lsEdgeOnBorder.distance(lsBorderTemp);
		if (distBorder < 0.1) {
			return lsBorderTemp;
		}
	}
	return ign::geometry::LineString();
}

///
///
///
bool app::calcul::CFeatGenerationOp::_isNextEdgeInAntennas(
	ign::feature::Feature const& fEdgeCurr,
	ign::geometry::Point const& ptCurr,
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
	ign::feature::FeatureIteratorPtr itNextEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*context->getFeatureStore(epg::EDGE), filterArroundNdNext);
	
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

	//--
	epg::Context* context = epg::ContextS::getInstance();

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	std::string const fictitiousFieldName = themeParameters->getValue(EDGE_FICTITIOUS_NAME).toString();

	if(_vCountriesCodeName.size() != 2)
		_logger->log(epg::log::WARN, "Attention, le code pays " + _countryCodeDouble + " n'est pas double" );

	std::string countryCode1 = _vCountriesCodeName[0];
	std::string countryCode2 = _vCountriesCodeName[1];

	//--
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + _countryCodeDouble + "'");
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCLCountryCode);
	boost::progress_display display(numFeatures, std::cout, "[ computing cl geometry by linked edges interpolation % complete ]\n");

	//--
	ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCLCountryCode);
	std::set<std::string> sCL2delete;
	while (it->hasNext()) {
		++display;
		ign::feature::Feature fCL = it->next();
		ign::geometry::LineString lsCLCurr = fCL.getGeometry().asLineString();
		lsCLCurr.setFillZ(0);
		
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

		ign::geometry::GeometryPtr geomEdg1(_getGeomProjClOnEdge(lsCLCurr, lsEdgInit1, snapProjCl2edge));
		ign::geometry::GeometryPtr geomEdg2(_getGeomProjClOnEdge(lsCLCurr, lsEdgInit2, snapProjCl2edge));

		//on coupe les edges au niveau de la projection des extremites des CLs sur ces edges, pour ne prendre que la portion d'edge que l'on apparie à la CL

		if (geomEdg1->isPoint() && geomEdg2->isPoint()) {
			_logger->log(epg::log::WARN, "Suppression CL " + fCL.getId());
			_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fCL);

			sCL2delete.insert(fCL.getId());
			continue;
		}

		bool isFictEdg1 = fEdg1.getAttribute(fictitiousFieldName).toString() == "true" ? true : false;
		bool isFictEdg2 = fEdg2.getAttribute(fictitiousFieldName).toString() == "true" ? true : false;

		if ( isFictEdg1 && geomEdg1->isPoint() || isFictEdg2 && geomEdg2->isPoint() ) {
			_logger->log(epg::log::WARN, "Suppression CL (géometrie cible ponctuelle)" + fCL.getId());
			_shapeLogger->writeFeature("ClDeletedNoCandidatefound", fCL);
			
			sCL2delete.insert(fCL.getId());
			continue;
		}

		ign::geometry::LineString lsCLUpdated;
		if (isFictEdg1 && !isFictEdg2)
			lsCLUpdated = geomEdg1->asLineString();
		else if (!isFictEdg1 && isFictEdg2)
			lsCLUpdated = geomEdg2->asLineString();
		else
			lsCLUpdated = _computeMeanGeom(*geomEdg1, *geomEdg2);

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
ign::geometry::LineString app::calcul::CFeatGenerationOp::_computeMeanGeom(
	ign::geometry::Geometry const& geom1,
	ign::geometry::Geometry const& geom2
) const {
	ign::geometry::LineString meanLs;

	std::set<double> sAbsCurv;

	geometry::tools::LengthIndexedLineString* lsIndex1Ptr = 0;
	geometry::tools::LengthIndexedLineString* lsIndex2Ptr = 0;

	double length1 = 0;
	double length2 = 0;

	if( geom1.isLineString() ) {
		ign::geometry::LineString const& ls = geom1.asLineString();
		length1 = ls.length();
		lsIndex1Ptr = new geometry::tools::LengthIndexedLineString(ls);
		for (size_t i = 0; i < ls.numPoints() - 1; ++i) {
			double abscurv = lsIndex1Ptr->getPointLocation(i) / length1;
			sAbsCurv.insert(abscurv);
		}
	}

	if( geom2.isLineString() ) {
		ign::geometry::LineString const& ls = geom2.asLineString();
		length2 = ls.length();
		lsIndex2Ptr = new geometry::tools::LengthIndexedLineString(ls);
		for (size_t i = 0; i < ls.numPoints() - 1; ++i) {
			double abscurv = lsIndex2Ptr->getPointLocation(i) / length2;
			sAbsCurv.insert(abscurv);
		}
	}

	for (std::set<double>::iterator sit = sAbsCurv.begin(); sit != sAbsCurv.end(); ++sit) {

		ign::geometry::Point pt1 = lsIndex1Ptr == 0 ? geom1.asPoint() : lsIndex1Ptr->locateAlong(*sit*length1);
		ign::geometry::Point pt2 = lsIndex2Ptr == 0 ? geom2.asPoint() : lsIndex2Ptr->locateAlong(*sit*length2);

		ign::geometry::MultiPoint multiPt;
		multiPt.addGeometry(pt1);
		multiPt.addGeometry(pt2);
		ign::geometry::Point meanPt = multiPt.getCentroid();
		
		if (sit != sAbsCurv.begin()) {
			if (meanLs.endPoint().distance(meanPt) >= 0) // WARNING : PAS DE SEUIL !!!
				meanLs.addPoint(meanPt);
		} else 
			meanLs.addPoint(meanPt);
	}

	//--
	ign::geometry::Point endPt1 = lsIndex1Ptr == 0 ? geom1.asPoint() : geom1.asLineString().endPoint();
	ign::geometry::Point endPt2 = lsIndex2Ptr == 0 ? geom2.asPoint() : geom2.asLineString().endPoint();

	ign::geometry::MultiPoint multiPtEnd;
	multiPtEnd.addGeometry(endPt1);
	multiPtEnd.addGeometry(endPt2);
	ign::geometry::Point meanEndPt = multiPtEnd.getCentroid();

	meanLs.addPoint(meanEndPt);
	meanLs.setFillZ(0);

	if (lsIndex1Ptr != 0 )
		delete lsIndex1Ptr;
	
	if (lsIndex2Ptr != 0 )
		delete lsIndex2Ptr;

	return meanLs;
}

///
///
///
ign::geometry::Geometry* app::calcul::CFeatGenerationOp::_getGeomProjClOnEdge(
	ign::geometry::LineString const& lsCl,
	ign::geometry::LineString const& lsEdge,
	double snapProjCl2edge
) const {
	ign::geometry::Point startLsProj = epg::tools::geometry::project(lsEdge, lsCl.startPoint(), snapProjCl2edge);
	startLsProj.setZ(0);//interpolate avec projZ
	ign::geometry::Point endLsProj = epg::tools::geometry::project(lsEdge, lsCl.endPoint(), snapProjCl2edge);
	endLsProj.setZ(0);//interpolate

	if( startLsProj.egal2d(endLsProj) )
		return startLsProj.clone();

	app::geometry::tools::LineStringSplitter lsSplitter(lsEdge);
	lsSplitter.addCuttingGeometry(startLsProj);
	lsSplitter.addCuttingGeometry(endLsProj);
	std::vector<ign::geometry::LineString> vLsCandidate = lsSplitter.getSubLineStringsZ();
	int idCandidate = -1;
	for (size_t i = 0; i < vLsCandidate.size(); ++i) {
		ign::geometry::LineString const& lsCandidate = vLsCandidate[i];
		if ((lsCandidate.startPoint().egal2d(startLsProj) && lsCandidate.endPoint().egal2d(endLsProj))
			|| (lsCandidate.startPoint().egal2d(endLsProj) && lsCandidate.endPoint().egal2d(startLsProj))) 
		{
			idCandidate = i;
			break;
		}
	}
	if (idCandidate < 0 || vLsCandidate[idCandidate].isEmpty())
		return startLsProj.clone();

	vLsCandidate[idCandidate].startPoint().setZ(0);
	vLsCandidate[idCandidate].endPoint().setZ(0);

	//on s'assure du sens de ls
	if (!vLsCandidate[idCandidate].startPoint().egal2d(startLsProj))
		return vLsCandidate[idCandidate].reverse().clone();
	else
		return vLsCandidate[idCandidate].clone();
}

// void app::calcul::CFeatGenerationOp::_deleteCLUnderThreshold() const
// {
// 	_logger->log(epg::log::TITLE, "[ BEGIN CLEAN CL UNDER THRESHOLD " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

// 	epg::Context* context = epg::ContextS::getInstance();
// 	std::string const CLTableName = _fsCL->getTableName();
// 	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
// 	std::string const idName = context->getEpgParameters().getValue(ID).toString();
// 	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

// 	std::ostringstream ssconditionDeleteCLUnderThreshold;

// 	ssconditionDeleteCLUnderThreshold << " ST_LENGTH(" << geomName << ") < 10 AND "<<countryCodeName << " ='" << _countryCodeDouble << "'";

// 	std::set<std::string> sCLToDelete;
// 	ign::feature::FeatureFilter filterCLInf10m(ssconditionDeleteCLUnderThreshold.str());
// 	ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCLInf10m);
// 	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCLInf10m);
// 	boost::progress_display display(numFeatures, std::cout, "[ cleaning cl under threshold % complete ]\n");

// 	while (it->hasNext()) {
// 		++display;
// 		ign::feature::Feature const& fCL10m = it->next();

// 		ign::feature::FeatureFilter filterNeighbor(idName + " <> '" + fCL10m.getId()+"'");
// 		filterNeighbor.setExtent(fCL10m.getGeometry().getEnvelope());

// 		ign::feature::FeatureIteratorPtr itArround = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterNeighbor);
// 		bool hasNeithbor = false;
// 		while (itArround->hasNext()) {
// 			ign::feature::Feature const& fCLNeighbor = itArround->next();
// 			if (fCLNeighbor.getGeometry().distance(fCL10m.getGeometry()) < 0.1) {
// 				hasNeithbor = true;
// 				break;
// 			}
// 		}
// 		if(!hasNeithbor) {
// 			sCLToDelete.insert(fCL10m.getId());
// 			_shapeLogger->writeFeature("clDeleteByLength", fCL10m);
// 		} else {
// 			_shapeLogger->writeFeature("clNoDeleteByLength", fCL10m);
// 		}
			
// 	}

// 	for (std::set<std::string>::iterator sit = sCLToDelete.begin(); sit != sCLToDelete.end(); ++sit) {
// 		_fsCL->deleteFeature(*sit);
// 	}

// 	_logger->log(epg::log::INFO, "Nb CL isolées inférieures a un seuil supprimées  " + ign::data::Integer(sCLToDelete.size()).toString());
// 	_logger->log(epg::log::TITLE, "[ END CLEAN CL UNDER THRESHOLD " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
// }

///
///
///
void app::calcul::CFeatGenerationOp::_deleteCLUnderThreshold() const
{
	_logger->log(epg::log::TITLE, "[ BEGIN CLEAN CL UNDER THRESHOLD " + _countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	//--
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	double minLength = themeParameters->getValue(CL_MIN_LENGTH).toDouble();

	//--
	GraphType graph;
	_loadGraphCL(graph);

	// patience
	boost::progress_display display(graph.numEdges(), std::cout, "[ deleting short cl % complete ]\n");

	//--
	std::set<std::string> sCLToDelete;
	std::set<edge_descriptor> sVisitedEdge;
	edge_iterator eit, eend;
	for( graph.edges( eit, eend ) ; eit != eend ; ++eit, ++display )
	{
		if( sVisitedEdge.find(*eit) != sVisitedEdge.end() )
			continue;

		edges_path path = _getPath(graph, *eit);

		if( graph.degree(graph.source(*path.begin())) > 1 || graph.degree(graph.target(*path.rbegin())) > 1 )
			continue;

		for (edges_path_const_iterator pit = path.begin() ; pit != path.end() ; ++pit)
			sVisitedEdge.insert(pit->descriptor);

		double length = _getLength(graph, path);

		if( length > minLength )
			continue;

		for (edges_path_const_iterator pit = path.begin() ; pit != path.end() ; ++pit) {
			sCLToDelete.insert(graph.origins(pit->descriptor)[0]);

			//--
			ign::feature::Feature feat;
			feat.setGeometry(graph.getGeometry(pit->descriptor));
			_shapeLogger->writeFeature("clDeleteByLength", feat);
		}
			
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
double app::calcul::CFeatGenerationOp::_getLength(
	GraphType const& graph,
	edges_path const& path
 ) const {
	double length = 0;
	for (edges_path_const_iterator pit = path.begin() ; pit != path.end() ; ++pit) {
		length += graph.getGeometry(pit->descriptor).length();
	}
	return length;
 }

///
///
///
app::calcul::CFeatGenerationOp::edges_path app::calcul::CFeatGenerationOp::_getPath(
	GraphType const& graph,
	edge_descriptor e
) const {
	edges_path path;

	oriented_edge_descriptor tPivot[] = { 
		oriented_edge_descriptor( e, ign::graph::DIRECT ), 
		oriented_edge_descriptor( e, ign::graph::REVERSE ) 
	};

	path.push_back(tPivot[0]);

	for( size_t  i = 0 ; i < 2 ; ++i) {
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

	return path;
}

///
///
///
void app::calcul::CFeatGenerationOp::_loadGraphCL(GraphType & graphCL) const
{
	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	//--
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + _countryCodeDouble + "'");
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCLCountryCode);
	boost::progress_display display(numFeatures, std::cout, "[ loading cl graph % complete ]\n");

	//--
	ign::geometry::graph::builder::SimpleGraphBuilder<GraphType> graphBuilder(graphCL, 0.01);
	ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCLCountryCode);
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
	std::string const& countryCodeSimple,
	GraphType & graphEdges
) const {
	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	//--
	ign::feature::FeatureFilter filterEdgeCountryCode(countryCodeName + " = '" + countryCodeSimple + "'");
	ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterEdgeCountryCode);
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterEdgeCountryCode);
	boost::progress_display display(numFeatures, std::cout, "[ loading "+ countryCodeSimple+" edges graph % complete ]\n");

	ign::geometry::graph::builder::SimpleGraphBuilder<GraphType> graphBuilder(graphEdges, 0.01);
	while (it->hasNext()) {
		++display;
		ign::feature::Feature fedge = it->next();
		graphBuilder.addEdge(fedge.getGeometry().asLineString(), fedge.getId());
	}
}

///
///
///
bool app::calcul::CFeatGenerationOp::_isConnectedEdges(
	GraphType const& graph,
	std::string const& idEdge1,
	std::string const& idEdge2
) const {
	edge_descriptor edCl1 = graph.getInducedEdges(idEdge1).second[0].descriptor;
	edge_descriptor edCl2 = graph.getInducedEdges(idEdge2).second[0].descriptor;

	if (graph.source(edCl1) == graph.source(edCl2) || graph.source(edCl1) == graph.target(edCl2) || graph.target(edCl1) == graph.source(edCl2) || graph.target(edCl1) == graph.target(edCl2))
		return true;

	return false;
}

///
///
///
std::pair<bool,std::pair<std::string, std::string>> app::calcul::CFeatGenerationOp::_getClLinkedEdges(
	std::string const& linkedFeatIdName,
	GraphType const& graphCL,
	GraphType::edge_descriptor eCl
) const {
	std::string idCl = graphCL.origins(eCl)[0];
	_logger->log(epg::log::DEBUG, idCl);

	ign::feature::Feature clFeat;
	_fsCL->getFeatureById(idCl, clFeat);

	if (clFeat.getId().empty())
		return std::make_pair(false, std::make_pair("", ""));

	std::string linkedFeat = clFeat.getAttribute(linkedFeatIdName).toString();
	_logger->log(epg::log::DEBUG, linkedFeat);

	std::vector<std::string> vEdgeslinkedJ;
	epg::tools::StringTools::Split(linkedFeat, "#", vEdgeslinkedJ);

	return std::make_pair(true,std::make_pair(vEdgeslinkedJ[0], vEdgeslinkedJ[1]));
}

///
///
///
bool app::calcul::CFeatGenerationOp::_areParallelEdges(
	GraphType const& graphCL,
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
void app::calcul::CFeatGenerationOp::_setContinuityCl(GraphType const& graphCL) const
{
	_logger->log(epg::log::TITLE, "[ BEGIN SET CL CONTINUITY ] : " + epg::tools::TimeTools::getTime());

	epg::Context* context = epg::ContextS::getInstance();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	GraphType graphEdges1, graphEdges2;
	_loadGraphEdges(_vCountriesCodeName[0], graphEdges1);
	_loadGraphEdges(_vCountriesCodeName[1], graphEdges2);

	GraphType::vertex_iterator vit, vitEnd;
	graphCL.vertices(vit, vitEnd);

	boost::progress_display display(graphCL.numVertices(), std::cout, "[ set cl continuity % complete ]\n");
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

	//--
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	//--
	GraphType graphClDoublon;
	ign::geometry::graph::tools::SnapRoundPlanarizer< GraphType >  planarizerClDoublon(graphClDoublon, 1e2);//scale =1e2 -> precision de 0.01

	//--
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + _countryCodeDouble + "'");
	size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCLCountryCode);
	boost::progress_display displayLoad(numFeatures, std::cout, "[ loading cl planar graph % complete ]\n");

	ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCLCountryCode);
	while (it->hasNext()) {
		++displayLoad;
		ign::feature::Feature fCL = it->next();
		planarizerClDoublon.addEdge(fCL.getGeometry().asLineString(), fCL.getId());
	}
	planarizerClDoublon.planarize();

	//--
	GraphType::edge_iterator eit, eitEnd;
	graphClDoublon.edges(eit, eitEnd);
	boost::progress_display display(graphClDoublon.numEdges(), std::cout, "[ deleting cl duplicates % complete ]\n");
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
			if (fClDoublon.getId().empty() )
				continue;

			_shapeLogger->writeFeature("ClDoublon", fClDoublon);
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

	boost::progress_display display(graph.numEdges(), std::cout, "[ gathering edges by origins % complete ]\n");

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
					//on c_mergingEdgesByOriginhange le pivot
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
	boost::progress_display display2(mGatheredEdges.size(), std::cout, "[ merging cl by origins % complete ]\n");

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
