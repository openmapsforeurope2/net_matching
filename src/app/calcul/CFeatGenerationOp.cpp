#include <app/calcul/CFeatGenerationOp.h>
#include <app/calcul/detail/ClMerger.h>
//#include <app/tools/CountryQueryErmTrans.h>
#include <app/params/ThemeParameters.h>
#include <app/geometry/tools/LengthIndexedLineString.h>

//BOOST
#include <boost/timer.hpp>
#include <boost/progress.hpp>

//SOCLE
#include <ign/geometry/algorithm/BufferOpGeos.h>
#include <ign/math/Line2T.h>
#include <ign/math/LineT.h>
#include <ign/geometry/graph/builder/SimpleGraphBuilder.h>
#include <ign/geometry/graph/tools/SnapRoundPlanarizer.h>
#include <ign/geometry/algorithm/OptimizedHausdorffDistanceOp.h>
//#include <ign/geometry/tools/LengthIndexedLineString.h>

//EPG
#include <epg/Context.h>
#include <epg/io/query/CountryQueriesReader.h>
#include <epg/io/query/utils/queryUtils.h>
#include <epg/tools/TimeTools.h>
#include <epg/tools/FilterTools.h>
#include <epg/tools/StringTools.h>
#include <epg/utils/CopyTableUtils.h>
#include <epg/utils/replaceTableName.h>
#include <epg/tools/geometry/LineStringSplitter.h>
#include <epg/tools/geometry/interpolate.h>
#include <epg/tools/geometry/project.h>
#include <epg/tools/geometry/getSubLineString.h>
#include <epg/tools/geometry/angle.h>
#include <epg/tools/geometry/SegmentIndexedGeometry.h>
#include <epg/tools/geometry/LineIntersector.h>

// BOOST
#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>






void app::calcul::CFeatGenerationOp::computeCL(std::string countryCodeDouble)
{
	_logger->log(epg::log::TITLE, "[ BEGIN CL GENERATION FOR " + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	epg::Context* context = epg::ContextS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const clTableName = _fsCL->getTableName();
	//CL
	double distBuffer = 20; // 5;
	double thresholdNoCL = 0.1;
	double ratioInBuff = 0.6;
	double snapOnVertexBorder = 2;
	double distMaxClClosest = 2;
	double angleMaxBorder = 25;
	double angleMaxEdges = 25;
	angleMaxBorder = angleMaxBorder * M_PI / 180;
	angleMaxEdges = angleMaxEdges * M_PI / 180;
	double distCLIntersected = 10;
	double distUnderShoot = 10;
	double distMergeCL = 1;
	double distMergeCP = 2;
	double distMaxEdges= 10; //5
	double snapProjCl2edge = 0.1;

	ign::feature::FeatureIteratorPtr itBoundary = _fsBoundary->getFeatures(ign::feature::FeatureFilter(countryCodeName + " = '" + countryCodeDouble + "'"));
	while (itBoundary->hasNext())
	{
		ign::feature::Feature fBoundary = itBoundary->next();

		std::string boundaryType = fBoundary.getAttribute("boundary_type").toString();

		// On ne traite que les frontières de type international_boundary ou coastline_sea_limit
		// On aurait pu filtrer en entrée mais le filtre semble très long, peut-être à cause de l'enum qui oblige à utiliser boundary_type::text like '%coastline_sea_limit%'
		if (boundaryType != "international_boundary" && boundaryType.find("coastline_sea_limit") == -1)
			continue;

		ign::geometry::LineString lsBoundary = fBoundary.getGeometry().asLineString();

		ign::geometry::algorithm::BufferOpGeos buffOp;
		ign::geometry::GeometryPtr buffBorder(buffOp.buffer(lsBoundary, distBuffer, 0, ign::geometry::algorithm::BufferOpGeos::CAP_FLAT));

		ign::feature::Feature fBuff;
		fBuff.setGeometry(buffBorder->clone());
		//_shapeLogger->writeFeature("getCLfromBorder_buffers", fBuff);

		_getCLfromBorder(lsBoundary, buffBorder, distBuffer, thresholdNoCL, angleMaxBorder, ratioInBuff, snapOnVertexBorder);


	}
	epg::utils::CopyTableUtils::copyTable(clTableName,idName,geomName, ign::geometry::Geometry::GeometryTypeLineString, clTableName+"_avtmerge", "", false,true);


	_mergeIntersectingClWithGraph(countryCodeDouble, distMaxEdges, snapProjCl2edge);
	epg::utils::CopyTableUtils::copyTable(clTableName, idName, geomName, ign::geometry::Geometry::GeometryTypeLineString, clTableName + "_apresmerge", "", false, true);

//	_snapCl2Cl(countryCodeDouble, distMaxClClosest);

	//_getClDoublonGeom(countryCodeDouble);
	_deleteClByAngleAndDistEdges(countryCodeDouble, angleMaxEdges, distMaxEdges, snapProjCl2edge);
	epg::utils::CopyTableUtils::copyTable(clTableName, idName, geomName, ign::geometry::Geometry::GeometryTypeLineString, clTableName + "_apresdeletebyangle", "", false, true);

	_deleteCLUnderThreshold(countryCodeDouble);
	epg::utils::CopyTableUtils::copyTable(clTableName, idName, geomName, ign::geometry::Geometry::GeometryTypeLineString, clTableName + "_apresdeletebyseuilmin", "", false, true);

	//boucle sur CC?
	GraphType graphCL;
	_loadGraphCL(countryCodeDouble,graphCL);

	_updateGeomCL(countryCodeDouble, snapProjCl2edge);
	epg::utils::CopyTableUtils::copyTable(clTableName, idName, geomName, ign::geometry::Geometry::GeometryTypeLineString, clTableName + "_apresupdate", "", false, true);

	_setContinuityCl(countryCodeDouble,graphCL);

	//_snapCl2Cl(countryCodeDouble, distMaxClClosest);

	_logger->log(epg::log::TITLE, "[ END CL GENERATION FOR " + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

}


void app::calcul::CFeatGenerationOp::computeCP(std::string countryCodeDouble)
{
	_logger->log(epg::log::TITLE, "[ BEGIN CP GENERATION FOR " + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	epg::Context* context = epg::ContextS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	//CL
	double distBuffer = 5;
	double thresholdNoCL = 0.1;
	double ratioInBuff = 0.6;
	//double snapOnVertexBorder = 1;
	double snapOnVertexBorder = 1;
	double angleMaxBorder = 25;
	angleMaxBorder = angleMaxBorder * M_PI / 180;
	double distCLIntersected = 10;
	double distUnderShoot = 10;
	double distMergeCL = 1;
	double distMergeCP = 5; //2


	ign::feature::FeatureIteratorPtr itBoundary = _fsBoundary->getFeatures(ign::feature::FeatureFilter(countryCodeName + " = '" + countryCodeDouble + "'"));
	while (itBoundary->hasNext()) {
		ign::feature::Feature fBoundary = itBoundary->next();
		std::string boundaryType = fBoundary.getAttribute("boundary_type").toString();
		if (boundaryType != "international_boundary" && boundaryType.find("coastline_sea_limit") == -1)
			continue;
		ign::geometry::LineString lsBoundary = fBoundary.getGeometry().asLineString();
		ign::geometry::algorithm::BufferOpGeos buffOp;
		ign::geometry::GeometryPtr buffBorder(buffOp.buffer(lsBoundary, distBuffer, 0, ign::geometry::algorithm::BufferOpGeos::CAP_FLAT));

		_getCPfromIntersectBorder(lsBoundary, distCLIntersected);

		_addToUndershootNearBorder(lsBoundary, buffBorder, distUnderShoot);
	}

	_snapCPNearBy(countryCodeDouble, distMergeCP, 0);

	_logger->log(epg::log::TITLE, "[ END CP GENERATION FOR " + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}

app::calcul::CFeatGenerationOp::CFeatGenerationOp( bool verbose)
{
	//debug//test
	_init( verbose);
}

app::calcul::CFeatGenerationOp::~CFeatGenerationOp()
{
	_shapeLogger->closeShape("CLBeforeMerge");
	_shapeLogger->closeShape("ClMergedBeforeUpdate");
	_shapeLogger->closeShape("ClDeletedNoCandidatefound");
	_shapeLogger->closeShape("ClDoublon");
	_shapeLogger->closeShape("ClDeleteByAngleDistEdges");
	_shapeLogger->closeShape("ClDebug");
}
	
///
///
///
void app::calcul::CFeatGenerationOp::_init( bool verbose)
{
	_logger = epg::log::EpgLoggerS::getInstance();
	_logger->log(epg::log::TITLE, "[ BEGIN INITIALIZATION ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	
	_shapeLogger = epg::log::ShapeLoggerS::getInstance();
	_shapeLogger->addShape("CLBeforeMerge", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClMergedBeforeUpdate", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDeletedNoCandidatefound", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDoublon", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDeleteByAngleDistEdges", epg::log::ShapeLogger::LINESTRING);
	_shapeLogger->addShape("ClDebug", epg::log::ShapeLogger::LINESTRING);

	_verbose = verbose;

	///recuperation de la liste des attributs à concatener dans la fusion des attributs
	std::string listAttr2concatName = themeParameters->getValue(LIST_ATTR_TO_CONCAT).toString();
	//on recup les attribut a concat et on les mets dans un vecteur en splitant
	std::vector<std::string> _vAttrNameToConcat;
	epg::tools::StringTools::Split(listAttr2concatName, "/", _vAttrNameToConcat);
	for (size_t i = 0; i < _vAttrNameToConcat.size(); ++i) {
		_sAttrNameToConcat.insert(_vAttrNameToConcat[i]);
	}

	///recuperation de la liste des attributs à concatener dans la fusion des attributs
	std::string listAttrWName = themeParameters->getValue(LIST_ATTR_W).toString();
	//on recup les attribut a concat et on les mets dans un vecteur en splitant
	std::vector<std::string> _vAttrNameW;
	epg::tools::StringTools::Split(listAttrWName, "/", _vAttrNameW);
	for (size_t i = 0; i < _vAttrNameW.size(); ++i) {
		_sAttrNameW.insert(_vAttrNameW[i]);
	}

	///recuperation des features
	std::string const boundaryTableName = epg::utils::replaceTableName(context->getEpgParameters().getValue(TARGET_BOUNDARY_TABLE).toString());
	_fsBoundary = context->getDataBaseManager().getFeatureStore(boundaryTableName, idName, geomName);

	std::string const landmaskTableName = epg::utils::replaceTableName(themeParameters->getValue(LANDMASK_TABLE).toString());
	_fsLandmask= context->getDataBaseManager().getFeatureStore(landmaskTableName, idName, geomName);

	_fsEdge = context->getFeatureStore(epg::EDGE);

	_reqFilterEdges2generateCF = themeParameters->getValue(SQL_FILTER_EDGES_2_GENERATE_CF).toString();

	std::string edgeTableName = context->getEpgParameters().getValue(EDGE_TABLE).toString();

	///Create tmp_cp table
	std::string cpTableName = epg::utils::replaceTableName(themeParameters->getValue(CP_TABLE).toString());
	// if (!context->getDataBaseManager().tableExists(cpTableName)) {
	{
		std::ostringstream ss;
		ss << "DROP TABLE IF EXISTS " << cpTableName << " ;";
		ss << "CREATE TABLE " << cpTableName
			<< " AS TABLE " << edgeTableName
			<< " WITH NO DATA;"
			<< "ALTER TABLE " << cpTableName << " ALTER COLUMN "
			<< geomName << " type geometry(PointZ, 0);"
			<< "ALTER TABLE " << cpTableName << " ALTER COLUMN "
			<< idName << " type varchar(255);"
			<< "ALTER TABLE " << cpTableName << " ALTER COLUMN "
			<< countryCodeName << " type varchar(255);"
			<< "ALTER TABLE " << cpTableName << " ALTER COLUMN "
			<< "w_national_identifier" << " type varchar(255);"
			<< "ALTER TABLE " << cpTableName << " ALTER COLUMN "
			<< "national_road_code" << " type varchar(255);"
			<< "ALTER TABLE " << cpTableName << " ALTER COLUMN "
			<< "european_route_number" << " type varchar(255);"
			<< "ALTER TABLE " << cpTableName << " ADD COLUMN " << context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString() << " character varying(255);";

		context->getDataBaseManager().getConnection()->update(ss.str());
	}
	// Create tmp_cl table
	std::string clTableName = epg::utils::replaceTableName(themeParameters->getValue(CL_TABLE).toString());
	// if (!context->getDataBaseManager().tableExists(clTableName)) {
	{
		std::ostringstream ss;
		ss << "DROP TABLE IF EXISTS " << clTableName << " ;";
		ss << "CREATE TABLE " << clTableName
			<< " AS TABLE " << edgeTableName
			<< " WITH NO DATA;"
			<< "ALTER TABLE " << clTableName << " ALTER COLUMN "
			<< geomName << " type geometry(LineString, 0);"
			<< "ALTER TABLE " << clTableName << " ALTER COLUMN "
			<< idName << " type varchar(255);"
			<< "ALTER TABLE " << clTableName << " ALTER COLUMN "
			<< countryCodeName << " type varchar(255);"
			<< "ALTER TABLE " << clTableName << " ALTER COLUMN "
			<< "w_national_identifier" << " type varchar(255);"
			<< "ALTER TABLE " << clTableName << " ALTER COLUMN "
			<< "national_road_code" << " type varchar(255);"
			<< "ALTER TABLE " << clTableName << " ALTER COLUMN "
			<< "european_route_number" << " type varchar(255);"
			<< "ALTER TABLE " << clTableName << " ADD COLUMN " << context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString() << " character varying(255);";
		//patch pour ne pas avoir des enums et eviter les soucis lors de la fusion des attributs	
		for (std::set<std::string>::iterator sit = _sAttrNameToConcat.begin(); sit != _sAttrNameToConcat.end(); ++sit)
			ss << "ALTER TABLE " << clTableName << " ALTER COLUMN " << *sit << " TYPE character varying(255);";
			
				
		context->getDataBaseManager().getConnection()->update(ss.str());
	}

	_fsCP = context->getDataBaseManager().getFeatureStore(cpTableName, idName, geomName);
	_fsCL = context->getDataBaseManager().getFeatureStore(clTableName, idName, geomName);

	// id generator
	_idGeneratorCP = epg::sql::tools::IdGeneratorInterfacePtr(epg::sql::tools::IdGeneratorFactory::getNew(*_fsCP, "CONNECTINGPOINT"));
	_idGeneratorCL = epg::sql::tools::IdGeneratorInterfacePtr(epg::sql::tools::IdGeneratorFactory::getNew(*_fsCL, "CONNECTINGLINE"));

	_logger->log(epg::log::TITLE, "[ END INITIALIZATION ] : " + epg::tools::TimeTools::getTime());
}



void app::calcul::CFeatGenerationOp::_getCLfromBorder(
	ign::geometry::LineString & lsBorder,
	ign::geometry::GeometryPtr& buffBorder,
	double distBuffer,
	double thresholdNoCL,
	double angleMax,
	double ratioInBuff,
	double snapOnVertexBorder
)
{
	epg::Context* context = epg::ContextS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	std::vector<ign::feature::FeatureAttributeType> listAttrEdge = _fsEdge->getFeatureType().attributes();

	ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + buffBorder->toString() + "'),3035))");
	//ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_GeomFromText('" + buffBorder->toString() + "'))");
	epg::tools::FilterTools::addAndConditions(filter, _reqFilterEdges2generateCF);

	ign::feature::FeatureIteratorPtr eit = _fsEdge->getFeatures(filter);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsEdge, filter);
	boost::progress_display display(numFeatures, std::cout, "[ CREATE CONNECTING LINES ]\n");

	ign::math::Vec2d vecBorder(lsBorder.endPoint().x() - lsBorder.startPoint().x(), lsBorder.endPoint().y() - lsBorder.startPoint().y());

	while (eit->hasNext())
	{
		++display;
		//
		ign::feature::Feature fEdge = eit->next();
		ign::geometry::LineString lsEdge = fEdge.getGeometry().asLineString();

		//debug
		if (fEdge.getId() == "7124d0d2-6e61-4321-b302-7d9ce6c40a1a")
			bool bStop = true;
	
		//ign::geometry::Geometry* geomIntersect = lsEdge.Intersection(*buff);
		std::vector<ign::geometry::LineString> vLsProjectedOnBorder;
		epg::tools::geometry::LineStringSplitter lsSplitter(lsEdge);
		lsSplitter.addCuttingGeometry(*buffBorder);
		std::vector<ign::geometry::LineString> subEdgesBorder = lsSplitter.getSubLineStrings();

		//pas d'intersection par le buffer
		if (subEdgesBorder.size() == 1) {			
			double angleEdgBorder = _getAngleEdgeWithBorder(lsEdge,lsBorder);		
			//si l'edge est "proche" on considere qu'il est entierement dans le buffer et longe la frontiere
			if (lsEdge.distance(lsBorder) < distBuffer && (angleEdgBorder < angleMax || angleEdgBorder > (M_PI - angleMax) ) ) {
				ign::geometry::LineString lsCL;
				_getGeomCL(lsCL, lsBorder, lsEdge.startPoint(), lsEdge.endPoint(), snapOnVertexBorder);
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
						_getGeomCL(lsCL, lsBorder, subEdgesBorder[numfirstSubInBuff].startPoint(), subEdgesBorder[numlastSubInBuff].endPoint(), snapOnVertexBorder);
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


void app::calcul::CFeatGenerationOp::_addToUndershootNearBorder(
	ign::geometry::LineString & lsBorder,
	ign::geometry::GeometryPtr& buffBorder,
	double distUnderShoot
)
{
	epg::Context* context = epg::ContextS::getInstance();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	std::vector<ign::feature::FeatureAttributeType> listAttrEdge = _fsEdge->getFeatureType().attributes();

	//ign::feature::FeatureFilter filterBuffBorder("ST_INTERSECTS(" + geomName + ", ST_GeomFromText('" + buffBorder->toString() + "'))");
	ign::feature::FeatureFilter filterBuffBorder("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + buffBorder->toString() + "'),3035))");
	epg::tools::FilterTools::addAndConditions(filterBuffBorder, _reqFilterEdges2generateCF);

	ign::feature::FeatureIteratorPtr eit = _fsEdge->getFeatures(filterBuffBorder);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsEdge, filterBuffBorder);
	boost::progress_display display(numFeatures, std::cout, "[ GET CONNECTING POINTS BY UNDERSHOOT LINES ]\n");

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
		std::vector< ign::geometry::Point > vPtIntersect;
		epg::tools::geometry::SegmentIndexedGeometry segIndexLsBorder(&lsBorder);
		////////////////////todo recup le bon segment de la frontiere
		epg::tools::geometry::LineIntersector::compute(ptClosestBorder, lsEdge.pointN(1), lsBorder, vPtIntersect);
		double distMin = 100000;
		for (std::vector< ign::geometry::Point >::iterator vit = vPtIntersect.begin(); vit != vPtIntersect.end(); ++vit) {
			double dist = ptClosestBorder.distance(*vit);
			if (dist < distMin) {
				projPt = *vit;
				distMin = dist;
			}
		}

		if (ptClosestBorder.distance(projPt) > distUnderShoot )
			continue;
		double distCLIntersected = distUnderShoot;
		if (_isEdgeIntersectedPtWithCL(fEdge, projPt, distCLIntersected))
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


void app::calcul::CFeatGenerationOp::_getCPfromIntersectBorder(
	ign::geometry::LineString & lsBorder,
	double distCLIntersected
)
{
	epg::Context* context = epg::ContextS::getInstance();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	std::vector<ign::feature::FeatureAttributeType> listAttrEdge = _fsEdge->getFeatureType().attributes();

	//ign::feature::FeatureFilter filterFeaturesToMatch("ST_INTERSECTS(" + geomName + ", ST_GeomFromText('" + lsBorder.toString() + "'))");
	ign::feature::FeatureFilter filterFeaturesToMatch("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + lsBorder.toString() + "'),3035))");
	epg::tools::FilterTools::addAndConditions(filterFeaturesToMatch, _reqFilterEdges2generateCF);
	ign::feature::FeatureIteratorPtr itFeaturesToMatch = _fsEdge->getFeatures(filterFeaturesToMatch);

	int numFeatures = context->getDataBaseManager().numFeatures(*_fsEdge, filterFeaturesToMatch);
	boost::progress_display display(numFeatures, std::cout, "[ CREATE CONNECTING POINTS ]\n");

	while (itFeaturesToMatch->hasNext())
	{
		++display;

		ign::feature::Feature fToMatch = itFeaturesToMatch->next();
		ign::geometry::LineString const& lsFToMatch = fToMatch.getGeometry().asLineString();
		ign::geometry::Geometry* geomPtr = lsFToMatch.Intersection(lsBorder);

		if (geomPtr->isNull())
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
			//si l'edge sert à une CL et, ne pas créer de CP?
			bool isConnectedToCL = _isEdgeIntersectedPtWithCL(fToMatch, geomPtr->asPoint(), distCLIntersected);
			if (!isConnectedToCL) {
				fCF.setGeometry(geomPtr->asPoint());
				std::string idCP = _idGeneratorCP->next();
				_fsCP->createFeature(fCF, idCP);
			}

		}

		else if (geomPtr->isGeometryCollection())
		{
			ign::geometry::GeometryCollection geomCollect = geomPtr->asGeometryCollection();
			for (size_t i = 0; i < geomCollect.numGeometries(); ++i)
			{
				if (geomCollect.geometryN(i).isPoint())
				{
					ign::geometry::Point ptIntersect = geomCollect.geometryN(i).asPoint();
					
					bool isConnectedToCL = _isEdgeIntersectedPtWithCL(fToMatch, ptIntersect, distCLIntersected);
					if (!isConnectedToCL) {
						fCF.setGeometry(ptIntersect);
						std::string idCP = _idGeneratorCP->next();
						_fsCP->createFeature(fCF, idCP);
					}
				}
			}
		}
	}

}

double app::calcul::CFeatGenerationOp::_getAngleEdgeWithBorder(
	ign::geometry::LineString& lsEdge,
	ign::geometry::LineString& lsBorder)
{
	double angleEdgWBorder;
	ign::math::Vec2d vecEdge(lsEdge.endPoint().x() - lsEdge.startPoint().x(), lsEdge.endPoint().y() - lsEdge.startPoint().y());

	ign::geometry::Point ptStartProjOnBorder= epg::tools::geometry::project(lsBorder, lsEdge.startPoint(), 0);
	ign::geometry::Point ptEndProjOnBorder = epg::tools::geometry::project(lsBorder, lsEdge.endPoint(), 0);
	
	ign::math::Vec2d vecBorder(ptEndProjOnBorder.x() - ptStartProjOnBorder.x(), ptEndProjOnBorder.y() - ptStartProjOnBorder.y());

	angleEdgWBorder = epg::tools::geometry::angle(vecBorder, vecEdge);

	return angleEdgWBorder;
	//	//double angleSubEdgBorder = epg::tools::geometry::angle(vecBorder, vecSubEdge);
}


void app::calcul::CFeatGenerationOp::_getGeomCL(
	ign::geometry::LineString& lsCL,
	ign::geometry::LineString& lsBorder,
	ign::geometry::Point ptStartToProject,
	ign::geometry::Point ptEndToProject,
	double snapOnVertexBorder
)
{
	std::pair< int, double > pairPtCurvStartCL = epg::tools::geometry::projectAlong(lsBorder, ptStartToProject, snapOnVertexBorder);
	std::pair< int, double > pairPtCurvEndCL = epg::tools::geometry::projectAlong(lsBorder, ptEndToProject, snapOnVertexBorder);
	//recup de la border entre ces points pour recup de la geom CL
	lsCL = epg::tools::geometry::getSubLineString(pairPtCurvStartCL, pairPtCurvEndCL, lsBorder, snapOnVertexBorder);
}




bool app::calcul::CFeatGenerationOp::_isEdgeIntersectedPtWithCL(
	ign::feature::Feature& fEdge,
	ign::geometry::Point ptIntersectBorder,
	double distCLIntersected 
)
{
	epg::Context* context = epg::ContextS::getInstance();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	std::string idLinkedEdge = fEdge.getId();
	ign::feature::FeatureFilter filterIntersectCL (linkedFeatIdName+" like '%" + idLinkedEdge +"%'");
	filterIntersectCL.setExtent(ptIntersectBorder.getEnvelope().expandBy(distCLIntersected));
	ign::feature::FeatureIteratorPtr itIntersectedCL = _fsCL->getFeatures(filterIntersectCL);

	if (itIntersectedCL->hasNext()) {
		//if(ptIntersectBorder.distance(itIntersectedCL->next().getGeometry()) < 0.1) //==0 ? ou rien, car proche de 1 on suppose que c'est de la CL l'intersection
		return true;
	}

	return false;
}

void app::calcul::CFeatGenerationOp::_snapCPNearBy(
	std::string countryCodeDouble,
	double distMergeCP,
	double snapOnVertexBorder
)
{

	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	ign::feature::FeatureFilter filterCP;
	std::vector<std::string> vCountriesCodeName;
	epg::tools::StringTools::Split(countryCodeDouble, "#", vCountriesCodeName);
	for (size_t i = 0; i < vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCP, countryCodeName + " = '" + vCountriesCodeName[i] + "'");
	}
	ign::feature::FeatureIteratorPtr itCP = _fsCP->getFeatures(filterCP);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCP, filterCP);
	boost::progress_display display(numFeatures, std::cout, "[ FUSION CONNECTING POINTS WITH #]\n");

	std::set<std::string> sCP2Snap;
	std::string separator = "#";

	while (itCP->hasNext())
	{
		++display;

		ign::feature::Feature fCPCurr = itCP->next();
		std::string idCP = fCPCurr.getId();

		if (sCP2Snap.find(idCP) != sCP2Snap.end())
			continue;

		std::map<std::string, ign::feature::Feature> mCPNear;
		bool hasNearestCP = _getNearestCP(fCPCurr, distMergeCP, mCPNear);

		std::list<std::string> lCp2Delete;

		if (!hasNearestCP) {
			_fsCP->deleteFeature(fCPCurr.getId());
			continue;
		}

		std::set<std::string> s1, s2;
		for(std::map<std::string, ign::feature::Feature>::iterator mit = mCPNear.begin(); mit != mCPNear.end(); ++mit) {
			if(mit->second.getAttribute(countryCodeName).toString() == vCountriesCodeName[0]) s1.insert(mit->first);
			else s2.insert(mit->first);
			sCP2Snap.insert(mit->first);
		}

		std::map<std::string, std::string> m1;
		for (std::set<std::string>::const_iterator sit1 = s1.begin() ; sit1 != s1.end() ; ++sit1) {
			ign::geometry::Point const& p1 = mCPNear[*sit1].getGeometry().asPoint();
			std::string samicopain = "";
			double distanceMax = distMergeCP;
			for (std::set<std::string>::const_iterator sit2 = s2.begin() ; sit2 != s2.end() ; ++sit2) {
				ign::geometry::Point const& p2 = mCPNear[*sit2].getGeometry().asPoint();
				double distance = p1.distance(p2);
				if (distance < distanceMax) {
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
			ign::geometry::Point const& p2 = mCPNear[*sit2].getGeometry().asPoint();
			std::string samicopain = "";
			double distanceMax = distMergeCP;
			for (std::set<std::string>::const_iterator sit1 = s1.begin() ; sit1 != s1.end() ; ++sit1) {
				ign::geometry::Point const& p1 = mCPNear[*sit1].getGeometry().asPoint();
				double distance = p2.distance(p1);
				if (distance < distanceMax) {
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
			filterBorderNearCP.setExtent(ptCentroidCP.getEnvelope().expandBy(distMergeCP));
			ign::geometry::LineString lsBorderClosest;
			double distMinBorder = 2 * distMergeCP;
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


void app::calcul::CFeatGenerationOp::_snapCPNearBy2(
	std::string countryCodeDouble,
	double distMergeCP,
	double snapOnVertexBorder
)
{

	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	ign::feature::FeatureFilter filterCP;
	std::vector<std::string> vCountriesCodeName;
	epg::tools::StringTools::Split(countryCodeDouble, "#", vCountriesCodeName);
	for (size_t i = 0; i < vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCP, countryCodeName + " = '" + vCountriesCodeName[i] + "'");
	}
	ign::feature::FeatureIteratorPtr itCP = _fsCP->getFeatures(filterCP);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCP, filterCP);
	boost::progress_display display(numFeatures, std::cout, "[ FUSION CONNECTING POINTS WITH #]\n");

	std::set<std::string> sCP2Snap;
	std::string separator = "#";

	while (itCP->hasNext())
	{
		++display;

		ign::feature::Feature fCPCurr = itCP->next();
		std::string idCP = fCPCurr.getId();

		if (sCP2Snap.find(idCP) != sCP2Snap.end())
			continue;

		std::map<std::string, ign::feature::Feature> mCPNear;
		bool hasNearestCP = _getNearestCP(fCPCurr, distMergeCP, mCPNear);
		if (hasNearestCP) {
			ign::geometry::MultiPoint multiPtCP;
			for (std::map<std::string, ign::feature::Feature>::iterator mit = mCPNear.begin(); mit != mCPNear.end(); ++mit) {
				sCP2Snap.insert(mit->first);
				multiPtCP.addGeometry(mit->second.getGeometry());
			}
			
			//geom
			ign::geometry::Point ptCentroidCP = multiPtCP.asMultiPoint().getCentroid();
			ign::feature::FeatureFilter filterBorderNearCP;// (countryCodeName + " = 'be#fr'");
			filterBorderNearCP.setExtent(ptCentroidCP.getEnvelope().expandBy(distMergeCP));
			ign::geometry::LineString lsBorderClosest;
			double distMinBorder = 2 * distMergeCP;
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

			//boucle sur mCPNear pour modif
			for (std::map<std::string, ign::feature::Feature>::iterator mit = mCPNear.begin(); mit != mCPNear.end(); ++mit) {
				ign::feature::Feature fCP2snap;
				_fsCP->getFeatureById(mit->first, fCP2snap);
				fCP2snap.setGeometry(ptCentroidOnBorderCP);
				_fsCP->modifyFeature(fCP2snap);
			}

		}
	}

}

bool app::calcul::CFeatGenerationOp::_getNearestCP(
	ign::feature::Feature fCP,
	double distMergeCP,
	std::map<std::string,ign::feature::Feature>& mCPNear
)
{
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

void app::calcul::CFeatGenerationOp::_addFeatAttributeMergingOnBorder(
	ign::feature::Feature& featMerged,
	ign::feature::Feature& featAttrToAdd,
	std::string separator
)
{
	ign::feature::FeatureType featTypMerged = featMerged.getFeatureType();
	//ign::feature::FeatureType featTyp2 = feat2.getFeatureType();
	//test si featTyp1 == featTyp2 sinon msg d'erreur
	std::vector<std::string> vAttrNames;
	featTypMerged.getAttributeNames(vAttrNames);

	for (size_t i = 0; i < vAttrNames.size(); ++i) {
		std::string attrValueMerged;
		std::string attrName = vAttrNames[i];
		if (_sAttrNameW.find(attrName) != _sAttrNameW.end()) //on ne fusionne pas les attributs de travail
			continue;
		std::string attrValueToMerge = featMerged.getAttribute(attrName).toString();
		std::string attrValueToAdd = featAttrToAdd.getAttribute(attrName).toString();
		if (attrName == "vertical_level")
			bool bt = true;
		if (attrValueToMerge != attrValueToAdd)
			attrValueMerged = attrValueToMerge + separator + attrValueToAdd;
		else
			attrValueMerged = attrValueToMerge;
	
		featMerged.setAttribute(attrName,ign::data::String(attrValueMerged));
	}
}

void  app::calcul::CFeatGenerationOp::_snapCl2Cl(
	std::string countryCodeDouble, 
	double distMaxClClosest
)
{
	_logger->log(epg::log::TITLE, "[ BEGIN SNAP CONNECTING LINES 2 CONNECTING LINES ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	ign::feature::FeatureFilter filterCL(countryCodeName + " = '" + countryCodeDouble + "'");

	ign::feature::FeatureIteratorPtr itCL = _fsCL->getFeatures(filterCL);
	int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCL);
	boost::progress_display displayLoad(numFeatures, std::cout, "[ SNAP CL 2 CL]\n");

	std::map<std::string,ign::feature::Feature> mFClModified;
	while (itCL->hasNext()) {
		++displayLoad;
		ign::feature::Feature fCL = itCL->next();
		ign::geometry::LineString lsCl = fCL.getGeometry().asLineString();
		
		ign::feature::Feature fCl2snapStart;
		bool isClosestStartCl2snapStart;
		bool hasClCandidate2snapStart = _hasClExtremityClose(countryCodeDouble, distMaxClClosest, fCL, lsCl.startPoint(), fCl2snapStart, isClosestStartCl2snapStart);

		if (hasClCandidate2snapStart) {
			ign::geometry::LineString lsCl2snapStart = fCl2snapStart.getGeometry().asLineString();
			ign::geometry::MultiPoint mp;
			mp.addGeometry(lsCl.startPoint());
			if(isClosestStartCl2snapStart)
				mp.addGeometry(lsCl2snapStart.startPoint());
			else 
				mp.addGeometry(lsCl2snapStart.endPoint());

			ign::geometry::Point pt2Modif = mp.asMultiPoint().getCentroid();

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
			
		ign::feature::Feature fCl2snapEnd;
		bool isClosestStartCl2snapEnd;
		bool hasClCandidate2snapEnd= _hasClExtremityClose(countryCodeDouble, distMaxClClosest, fCL, lsCl.endPoint(), fCl2snapEnd, isClosestStartCl2snapEnd);
		if (hasClCandidate2snapEnd) {
			ign::geometry::LineString lsCl2snapEnd = fCl2snapEnd.getGeometry().asLineString();
			ign::geometry::MultiPoint mp;
			mp.addGeometry(lsCl.endPoint());
			if (isClosestStartCl2snapEnd)
				mp.addGeometry(lsCl2snapEnd.startPoint());
			else
				mp.addGeometry(lsCl2snapEnd.endPoint());

			ign::geometry::Point pt2Modif = mp.asMultiPoint().getCentroid();

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

	for (std::map<std::string, ign::feature::Feature>::iterator mit = mFClModified.begin(); mit != mFClModified.end(); ++mit) {
		_fsCL->modifyFeature(mit->second);
	}

}


bool  app::calcul::CFeatGenerationOp::_hasClExtremityClose(
	std::string countryCodeDouble,
	double distMaxClClosest,
	ign::feature::Feature fClCurr,
	ign::geometry::Point ptClCurr,
	ign::feature::Feature& fCl2snap,
	bool& isClosestStartCl2snap 
)
{
	//recup des CL proches si dist = 0 rien faire, si dist inf distMaxClClosest snapper
	epg::Context* context = epg::ContextS::getInstance();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	ign::feature::FeatureFilter filterArroundCl (idName + " <> '" + fClCurr.getId() + "'");
	epg::tools::FilterTools::addAndConditions(filterArroundCl, countryCodeName + " = '" + countryCodeDouble + "'");
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



void app::calcul::CFeatGenerationOp::_mergeIntersectingClWithGraph(
	std::string countryCodeDouble,
	double distMaxEdges,
	double snapProjCl2edge
)
{
	_logger->log(epg::log::TITLE, "[ BEGIN FUSION CONNECTING LINES ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	std::string separator = "#";

	GraphType graphCl;
	ign::geometry::graph::tools::SnapRoundPlanarizer< GraphType >  planarizerCl(graphCl, 10000);//scale =100 -> precision de 0.01
	ign::feature::FeatureFilter filterCL;
	std::vector<std::string> vCountriesCodeName;
	epg::tools::StringTools::Split(countryCodeDouble, "#", vCountriesCodeName);
	for (size_t i = 0; i < vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCL, countryCodeName + " = '" + vCountriesCodeName[i] + "'");
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

		//debug
		if (vClOrigins[0] == "CONNECTINGLINE2" && vClOrigins[1] == "CONNECTINGLINE5")
			bool bStop = true;

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
			if (countryCodeCl == vCountriesCodeName[0])
				mIdClOriginsCountry1[fCl.getId()] = fCl;
			else if (countryCodeCl == vCountriesCodeName[1])
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
			sEdgesMerged.insert(cl2merge.first);
			sEdgesMerged.insert(cl2merge.second);
			ign::feature::Feature fClNew = _fsCL->newFeature();
			fClNew = mIdClOriginsCountry1.find(cl2merge.first)->second;
			_addFeatAttributeMergingOnBorder(fClNew, mIdClOriginsCountry2.find(cl2merge.second)->second, separator);

			std::string idCLNew = _idGeneratorCL->next();
			fClNew.setId(idCLNew);
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
	detail::ClMerger::mergeAll(_idGeneratorCL.get(), _fsCL, mergeFilter);

	_logger->log(epg::log::TITLE, "[ END FUSION CONNECTING LINES ] : " + epg::tools::TimeTools::getTime());
}


void app::calcul::CFeatGenerationOp::_mergeIntersectingCL(
	std::string countryCodeDouble,
	double distMergeCL,
	double snapOnVertexBorder
)
{
	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	ign::feature::FeatureFilter filterCL;
	std::vector<std::string> vCountriesCodeName;
	epg::tools::StringTools::Split(countryCodeDouble, "#", vCountriesCodeName);
	for (size_t i = 0; i < vCountriesCodeName.size(); ++i) {
		epg::tools::FilterTools::addOrConditions(filterCL, countryCodeName + " = '" + vCountriesCodeName[i] + "'");
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

			ign::geometry::LineString lsIntersectedCL, lsSE, lsSS, lsES, lsEE;;

			_getGeomCL(lsSE, lsBorder, lsCurr.startPoint(), lsClArround.endPoint(), snapOnVertexBorder);
			_getGeomCL(lsSS, lsBorder, lsCurr.startPoint(), lsClArround.startPoint(), snapOnVertexBorder);
			_getGeomCL(lsES, lsBorder, lsCurr.endPoint(), lsClArround.startPoint(), snapOnVertexBorder);
			_getGeomCL(lsEE, lsBorder, lsCurr.endPoint(), lsClArround.endPoint(), snapOnVertexBorder);
			
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
				_addFeatAttributeMergingOnBorder(fCLNew, fCLCurr, separator);
			}
			else {
				fCLNew = fCLCurr;
				_addFeatAttributeMergingOnBorder(fCLNew, fCLArround, separator);
			}

			std::string idCLNew = _idGeneratorCL->next();
			fCLNew.setId(idCLNew);
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



void app::calcul::CFeatGenerationOp::_deleteClByAngleAndDistEdges(std::string countryCodeDouble, double angleMax, double distMax, double snapProjCl2edge)
{
	
	_logger->log(epg::log::TITLE, "[ BEGIN DELETE CL BY ANGLE EDGES FOR : " + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + countryCodeDouble + "'");

	int numCl2delete = -1;
	while (numCl2delete != 0){
		std::set<std::string> sCl2delete;
		GraphType graphCl;
		_loadGraphCL(countryCodeDouble, graphCl);
		ign::feature::FeatureIteratorPtr it = _fsCL->getFeatures(filterCLCountryCode);
		int numFeatures = context->getDataBaseManager().numFeatures(*_fsCL, filterCLCountryCode);
		boost::progress_display display(numFeatures, std::cout, "[  DELETE CL BY ANGLE EDGES ]\n");
		while (it->hasNext()) {
			++display;
			ign::feature::Feature fCl = it->next();
			//verifier si la cl n'est pas liée a au moins une cl à chaque extremite
			edge_descriptor edCl = graphCl.getInducedEdges(fCl.getId()).second[0].descriptor;
			if (graphCl.degree(graphCl.source(edCl)) > 1 && graphCl.degree(graphCl.target(edCl)) > 1)
				continue;

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
	_logger->log(epg::log::TITLE, "[ END DELETE CL BY ANGLE EDGES FOR :" + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
}


bool app::calcul::CFeatGenerationOp::_getCLToMerge(
	ign::feature::Feature fCL,
	double distMergeCL,
	std::map < std::string, ign::feature::Feature>& mCL2merge,
	std::set<std::string>& sCountryCode
)
{
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


void app::calcul::CFeatGenerationOp::_getBorderFromEdge(
	ign::geometry::LineString& lsEdgeOnBorder,
	ign::geometry::LineString& lsBorder
)
{
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



bool app::calcul::CFeatGenerationOp::_isNextEdgeInAntennas(ign::feature::Feature& fEdgeCurr, ign::geometry::Point& ptCurr, ign::feature::Feature&  edgeNext, ign::geometry::Point& ptNext)
{
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


void app::calcul::CFeatGenerationOp::_updateGeomCL(std::string countryCodeDouble, double snapProjCl2edge)
{
	_logger->log(epg::log::TITLE, "[ BEGIN UPDATE GEOM CL " + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const edgeTableName = _fsEdge->getTableName();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
	//std::string const natIdName = themeParameters->getValue(NATIONAL_IDENTIFIER).toString();

	std::vector<std::string> vCountriesCodeName;
	epg::tools::StringTools::Split(countryCodeDouble, "#", vCountriesCodeName);

	if(vCountriesCodeName.size() != 2)
		_logger->log(epg::log::WARN, "Attention, le countryCode " + countryCodeDouble + " n'a pas deux country" );

	std::string countryCode1 = vCountriesCodeName[0];
	std::string countryCode2 = vCountriesCodeName[1];

	//ign::geometry::MultiPolygon mPolyCountry1, mPolyCountry2;
	//_getGeomCountry(countryCode1, mPolyCountry1);
	//_getGeomCountry(countryCode2, mPolyCountry2);

	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + countryCodeDouble + "'");
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
			for (size_t i = 0; i < lsEdg1.numPoints()-1; ++i) {
				double abscurv = lsIndex1.getPointLocation(i)/lsEdg1.length();
				sAbsCurv.insert(abscurv);
			}
			for (size_t i = 0; i < lsEdg2.numPoints()-1; ++i) {
				double abscurv = lsIndex2.getPointLocation(i)/lsEdg2.length();
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
				if(!hasPtDistMin)
					lsCLUpdated.addPoint(ptLsCentroid);
			}
			ign::geometry::MultiPoint multiPtEnd;
			multiPtEnd.addGeometry(lsEdg1.endPoint());
			multiPtEnd.addGeometry(lsEdg2.endPoint());
			//multiPtEnd.addGeometry(endLsProj2);
			lsCLUpdated.addPoint(multiPtEnd.getCentroid());
//		}
		lsCLUpdated.clearZ();
		fCL.setGeometry(lsCLUpdated);
		_fsCL->modifyFeature(fCL);
		
	}

	for( std::set<std::string>::iterator sit = sCL2delete.begin(); sit != sCL2delete.end();++sit)
		_fsCL->deleteFeature(*sit);

	_logger->log(epg::log::TITLE, "[ BEGIN UPDATE GEOM CL " + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

}


void app::calcul::CFeatGenerationOp::_getGeomProjClOnEdge(ign::geometry::LineString& lsCl, ign::geometry::LineString& lsEdge, ign::geometry::LineString& lsprojClOnEdg, double snapProjCl2edge)
{
	ign::geometry::Point startLsProj = epg::tools::geometry::project(lsEdge, lsCl.startPoint(), snapProjCl2edge);
	startLsProj.setZ(0);//interpolate avec projZ
	ign::geometry::Point endLsProj = epg::tools::geometry::project(lsEdge, lsCl.endPoint(), snapProjCl2edge);
	endLsProj.setZ(0);//interpolate

	epg::tools::geometry::LineStringSplitter lsSplitter1(lsEdge);
	lsSplitter1.addCuttingGeometry(startLsProj);
	lsSplitter1.addCuttingGeometry(endLsProj);
	std::vector<ign::geometry::LineString> vLsCandidate1 = lsSplitter1.getSubLineStrings();
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

void app::calcul::CFeatGenerationOp::_deleteCLUnderThreshold(std::string countryCodeDouble)
{

	_logger->log(epg::log::TITLE, "[ BEGIN CLEAN CL UNDER THRESHOLD " + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

	epg::Context* context = epg::ContextS::getInstance();
	std::string const CLTableName = _fsCL->getTableName();
	std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
	std::string const idName = context->getEpgParameters().getValue(ID).toString();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

	std::ostringstream ssconditionDeleteCLUnderThreshold;

	ssconditionDeleteCLUnderThreshold << " ST_LENGTH(" << geomName << ") < 10 AND "<<countryCodeName << " ='" << countryCodeDouble << "'";

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
	_logger->log(epg::log::TITLE, "[ END CLEAN CL UNDER THRESHOLD " + countryCodeDouble + " ] : " + epg::tools::TimeTools::getTime());

}

void app::calcul::CFeatGenerationOp::_getGeomCountry(std::string countryCodeSimple, ign::geometry::MultiPolygon& geomCountry)
{
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


void app::calcul::CFeatGenerationOp::_loadGraphCL(std::string countryCodeDouble, GraphType& graphCL)
{
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + countryCodeDouble + "'");
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


void app::calcul::CFeatGenerationOp::_loadGraphEdges(std::string countryCodeSimple, GraphType& graphEdges)
{
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

bool app::calcul::CFeatGenerationOp::_isConnectedEdges(GraphType& graph, std::string idEdge1, std::string idEdge2)
{
	edge_descriptor edCl1 = graph.getInducedEdges(idEdge1).second[0].descriptor;
	edge_descriptor edCl2 = graph.getInducedEdges(idEdge2).second[0].descriptor;

	if (graph.source(edCl1) == graph.source(edCl2) || graph.source(edCl1) == graph.target(edCl2) || graph.target(edCl1) == graph.source(edCl2) || graph.target(edCl1) == graph.target(edCl2))
		return true;

	return false;
}

std::pair<std::string, std::string> app::calcul::CFeatGenerationOp::_getClLinkedEdges( std::string const& linkedFeatIdName, GraphType& graphCL, GraphType::edge_descriptor eCl ) {
	std::string idCl = graphCL.origins(eCl)[0];
	ign::feature::Feature clFeat;
	_fsCL->getFeatureById(idCl, clFeat);
	std::vector<std::string> vEdgeslinkedJ;
	epg::tools::StringTools::Split(clFeat.getAttribute(linkedFeatIdName).toString(), "#", vEdgeslinkedJ);

	return std::make_pair(vEdgeslinkedJ[0], vEdgeslinkedJ[1]);
}

bool app::calcul::CFeatGenerationOp::_areParallelEdges( GraphType& graphCL, GraphType::edge_descriptor e1,  GraphType::edge_descriptor e2 ) {
	if (graphCL.source(e1) == graphCL.source(e2) && graphCL.target(e1) == graphCL.target(e2) ) return true;
	if (graphCL.target(e1) == graphCL.source(e2) && graphCL.source(e1) == graphCL.target(e2) ) return true;
	return false;
}

ign::geometry::Point app::calcul::CFeatGenerationOp::_getLinkedEdgesConnectingPoint( GraphType const& graph, std::string const& idEdge1, std::string const& idEdge2 ){
	edge_descriptor edCl1 = graph.getInducedEdges(idEdge1).second[0].descriptor;
	edge_descriptor edCl2 = graph.getInducedEdges(idEdge2).second[0].descriptor;

	ign::feature::Feature fEdge1;
	_fsEdge->getFeatureById(idEdge1, fEdge1);
	ign::geometry::LineString const& edgeGeom1 = fEdge1.getGeometry().asLineString();

	if ( graph.source(edCl1) == graph.source(edCl2) || graph.source(edCl1) == graph.target(edCl2) )
		return edgeGeom1.startPoint();

	return edgeGeom1.endPoint();
}

void app::calcul::CFeatGenerationOp::_setContinuityCl(std::string countryCodeDouble, GraphType& graphCL)
{
	_logger->log(epg::log::TITLE, "[ BEGIN SET CL CONTINUITY ] : " + epg::tools::TimeTools::getTime());

	epg::Context* context = epg::ContextS::getInstance();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	std::vector<std::string> vCountriesCodeName;
	epg::tools::StringTools::Split(countryCodeDouble, "#", vCountriesCodeName);

	GraphType graphEdges1, graphEdges2;
	_loadGraphEdges(vCountriesCodeName[0], graphEdges1);
	_loadGraphEdges(vCountriesCodeName[1], graphEdges2);

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

			std::pair<std::string, std::string> pLinkedEdgesI = _getClLinkedEdges(linkedFeatIdName, graphCL, vClsIncidentTemp[i].descriptor);

			for (size_t j = i + 1; j < vClsIncidentTemp.size(); ++j) {
				if ( sTreated.find(j) != sTreated.end() ) continue;

				std::pair<std::string, std::string> pLinkedEdgesJ = _getClLinkedEdges(linkedFeatIdName, graphCL, vClsIncidentTemp[j].descriptor);

				bool isConnected1 = pLinkedEdgesI.first == pLinkedEdgesJ.first || _isConnectedEdges(graphEdges1, pLinkedEdgesI.first, pLinkedEdgesJ.first);
				bool isConnected2 = pLinkedEdgesI.second == pLinkedEdgesJ.second || _isConnectedEdges(graphEdges2, pLinkedEdgesI.second, pLinkedEdgesJ.second);

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
					mpLinkedEdgesConnectingPoints.addGeometry(_getLinkedEdgesConnectingPoint(graphEdges1, pLinkedEdgesI.first, pLinkedEdgesJ.first));
					mpLinkedEdgesConnectingPoints.addGeometry(_getLinkedEdgesConnectingPoint(graphEdges2, pLinkedEdgesI.second, pLinkedEdgesJ.second));
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
				if (mClModified.find(idClToModify) != mClModified.end())
					fClToModify = mClModified.find(idClToModify)->second;
				else {
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



void app::calcul::CFeatGenerationOp::_getClDoublonGeom(std::string countryCodeDouble)
{
	_logger->log(epg::log::TITLE, "[ BEGIN DELETE CL DOUBLON ] : " + epg::tools::TimeTools::getTime());
	epg::Context* context = epg::ContextS::getInstance();
	std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
	std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

	std::set<std::string> sClDoublonToDelete;

	GraphType graphClDoublon;
	ign::geometry::graph::tools::SnapRoundPlanarizer< GraphType >  planarizerClDoublon(graphClDoublon,100);//scale =100 -> precision de 0.01
	ign::feature::FeatureFilter filterCLCountryCode(countryCodeName + " = '" + countryCodeDouble + "'");
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