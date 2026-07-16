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
#include <epg/tools/geometry/interpolate.h>
#include <epg/tools/geometry/project.h>
#include <epg/tools/geometry/getSubLineString.h>
#include <epg/tools/geometry/angle.h>
#include <epg/tools/geometry/LineIntersector.h>
#include <epg/utils/CopyTableUtils.h>
#include<epg/graph/tools/merge.h>
#include<epg/graph/tools/reverse.h>

// OME2
#include <ome2/calcul/detail/ClMerger.h>
#include <ome2/feature/sql/NotDestroyedTools.h>

namespace app
{
    namespace calcul
    {
		///
		///
		///
		CFeatGenerationOp::CFeatGenerationOp(
			std::string const& borderCode,
			bool verbose
		) : _borderCode(borderCode),
			_verbose(verbose),
			_tagFromCl("from_cl"),
			_mlsBorderSmoothed(0),
			_attrMerger(0)
		{
			_init();
		}

		///
		///
		///
		CFeatGenerationOp::~CFeatGenerationOp()
		{
			_shapeLogger->closeShape("CLBeforeMerge");
			_shapeLogger->closeShape("ClDeletedNoCandidatefound");
			_shapeLogger->closeShape("ClDeleteByAngleDistEdges");
			_shapeLogger->closeShape("ClDebug");
			_shapeLogger->closeShape("edgeClCutByCp");
			_shapeLogger->closeShape("lsBorderCutByAngle");
			_shapeLogger->closeShape("clSnapclPoint");
			_shapeLogger->closeShape("clDeleteByLength");
			_shapeLogger->closeShape("debug_merging_of_2_cp_on_cl");
			_shapeLogger->closeShape("debug_merging_of_cp_on_cl_and_cp_border");
			_shapeLogger->closeShape("debug_duplicate");

			if (_mlsBorderSmoothed)
				delete _mlsBorderSmoothed;

			if (_attrMerger)
				delete _attrMerger;
		}

		///
		///
		///
		void CFeatGenerationOp::GenerateConnectingLinesByCountry(
			std::string const& borderCode,
			bool verbose
		) {
			CFeatGenerationOp op(borderCode, verbose);
			op._generateConnectingLinesByCountry();
		}

		///
		///
		///
		void CFeatGenerationOp::_init() {
			//--
			_logger = epg::log::EpgLoggerS::getInstance();
			_logger->log(epg::log::INFO, "[START] initialization: " + epg::tools::TimeTools::getTime());

			//--
			epg::Context* context = epg::ContextS::getInstance();

			//--
			epg::params::EpgParameters const& epgParams = context->getEpgParameters();
			std::string const idName = epgParams.getValue(ID).toString();
			std::string const geomName = epgParams.getValue(GEOM).toString();
			std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
			std::string const edgeTableName = epgParams.getValue(EDGE_TABLE).toString();
			std::string const boundaryTableName = epgParams.getValue(TARGET_BOUNDARY_TABLE).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const cpTableName = themeParameters->getValue(CP_TABLE).toString();
			std::string const clTableName = themeParameters->getValue(CL_TABLE).toString();
			std::string const fromOfWayException = themeParameters->getValue(CP_FORM_OF_WAY_EXCEPTION).toString();
			std::string smoothedBoundaryTableName = themeParameters->getValue(BOUNDARY_SMOOTHED_TABLE).toString();
			if (smoothedBoundaryTableName == "")
				smoothedBoundaryTableName = boundaryTableName;
			std::string const listAttrSeparator = themeParameters->getValue(LIST_ATTR_SEPARATOR).toString();
			std::string const listAttrJsonName = themeParameters->getValue(LIST_ATTR_JSON).toString();
			std::string const listAttrWName = themeParameters->getValue(LIST_ATTR_W).toString();
			
			//--
			_shapeLogger = epg::log::ShapeLoggerS::getInstance();
			_shapeLogger->addShape("CLBeforeMerge", epg::log::ShapeLogger::LINESTRING);
			_shapeLogger->addShape("ClDeletedNoCandidatefound", epg::log::ShapeLogger::LINESTRING);
			_shapeLogger->addShape("ClDeleteByAngleDistEdges", epg::log::ShapeLogger::LINESTRING);
			_shapeLogger->addShape("ClDebug", epg::log::ShapeLogger::LINESTRING);
			_shapeLogger->addShape("edgeClCutByCp", epg::log::ShapeLogger::LINESTRING); 
			_shapeLogger->addShape("lsBorderCutByAngle", epg::log::ShapeLogger::LINESTRING);
			_shapeLogger->addShape("clSnapclPoint", epg::log::ShapeLogger::POINT);
			_shapeLogger->addShape("clDeleteByLength", epg::log::ShapeLogger::LINESTRING);
			_shapeLogger->addShape("debug_merging_of_2_cp_on_cl", epg::log::ShapeLogger::POINT);
			_shapeLogger->addShape("debug_merging_of_cp_on_cl_and_cp_border", epg::log::ShapeLogger::POINT);
			_shapeLogger->addShape("debug_duplicate", epg::log::ShapeLogger::LINESTRING);

			//--
			epg::tools::StringTools::Split(_borderCode, "#", _vCountry);

			//--
			_isClStatement = countryCodeName + " LIKE '%" + _vCountry.front() + "%' AND "+ countryCodeName +" LIKE '%" + _vCountry.back() + "%'";

			//--
			_mIsCountryStatement.insert(std::make_pair(_vCountry.front(), countryCodeName + " LIKE '%" + _vCountry.front() + "%' AND "+ countryCodeName +" NOT LIKE '%" + _vCountry.back() + "%'"));
			_mIsCountryStatement.insert(std::make_pair(_vCountry.back(), countryCodeName + " LIKE '%" + _vCountry.back() + "%' AND "+ countryCodeName +" NOT LIKE '%" + _vCountry.front() + "%'"));

			//--
			_attrMerger = new ome2::calcul::utils::AttributeMerger(
				listAttrWName,
				listAttrJsonName,
				countryCodeName,
				listAttrSeparator,
				"#"
			);
			
			//--
			std::vector<std::string> vValueFormOfWay;
			epg::tools::StringTools::Split(fromOfWayException, "/", vValueFormOfWay);
			for (size_t i = 0; i < vValueFormOfWay.size(); ++i) {
				_sFormOfWayException.insert(vValueFormOfWay[i]);
			}

			//--
			_sqlFilterForCpGeneration = themeParameters->getValue(SQL_FILTER_EDGES_2_GENERATE_CF).toString();

			//--
			_fsBoundary = context->getDataBaseManager().getFeatureStore(boundaryTableName, idName, geomName);
			_fsEdge = context->getDataBaseManager().getFeatureStore(edgeTableName, idName, geomName);
			_fsCP = context->getDataBaseManager().getFeatureStore(cpTableName, idName, geomName);
			_fsCL = context->getDataBaseManager().getFeatureStore(clTableName, idName, geomName);

			// id generator
			_idGeneratorCP = epg::sql::tools::IdGeneratorInterfacePtr(epg::sql::tools::IdGeneratorFactory::getNew(*_fsCP, "CONNECTINGPOINT"));
			_idGeneratorCL = epg::sql::tools::IdGeneratorInterfacePtr(epg::sql::tools::IdGeneratorFactory::getNew(*_fsCL, "CONNECTINGLINE"));

			//--
			ign::feature::FeatureFilter filter(countryCodeName + " = '" + _borderCode + "'");
			_mlsBorderSmoothed = new epg::tools::MultiLineStringTool(filter, *context->getDataBaseManager().getFeatureStore(smoothedBoundaryTableName, idName, geomName));
			
			//--
			_logger->log(epg::log::TITLE, "[ END INITIALIZATION ] : " + epg::tools::TimeTools::getTime());
		}

		///
		///
		///
		void CFeatGenerationOp::_generateConnectingLinesByCountry() const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN CL GENERATION BY COUNTRY FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			//--
			epg::Context* context = epg::ContextS::getInstance();

			//--
			std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			double const distBuffer = themeParameters->getValue(CL_BUFFER_DIST).toDouble();
			double angleMaxToCutBorder = themeParameters->getValue(CFG_BOUNDARY_ANGLE_THRESHOLD).toDouble();
			angleMaxToCutBorder = angleMaxToCutBorder * M_PI / 180;

			//--
			ign::feature::FeatureIteratorPtr itBoundary = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsBoundary, ign::feature::FeatureFilter(countryCodeName + " = '" + _borderCode + "'"));
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
				for (size_t i = 0; i < vLsBorderCutByAngle.size(); ++i)
				{
					ign::geometry::LineString lsBoundaryCutByAngle = vLsBorderCutByAngle[i];

					ign::geometry::algorithm::BufferOpGeos buffOp;
					ign::geometry::GeometryPtr buffBorder(buffOp.buffer(lsBoundaryCutByAngle, distBuffer, 0, ign::geometry::algorithm::BufferOpGeos::CAP_FLAT));

					_getCLfromBorder(lsBoundaryCutByAngle, *buffBorder);
				}
			}

			_logger->log(epg::log::TITLE, "[ END CL GENERATION BY COUNTRY FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
		}

		///
		///
		///
		void CFeatGenerationOp::MergeConnectingLinesOnBorder(
			std::string const& borderCode,
			bool verbose
		) {
			CFeatGenerationOp op(borderCode, verbose);
			op._mergeConnectingLinesOnBorder();
		}

		///
		///
		///
		void CFeatGenerationOp::_mergeConnectingLinesOnBorder() const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN CL MERGING FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			double const distMaxEdges = themeParameters->getValue(CL_EDGE_MAX_DIST).toDouble();
			double const snapProjCl2edge = themeParameters->getValue(CL_SNAP_PROJ_CL_2_EDGE_DIST).toDouble();

			//--
			_mergeIntersectingClWithGraph(distMaxEdges, snapProjCl2edge);

			_logger->log(epg::log::TITLE, "[ END CL MERGING FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
		}

		///
		///
		///
		void CFeatGenerationOp::SnapConnectingLines(
			std::string const& borderCode,
			bool verbose
		) {
			CFeatGenerationOp op(borderCode, verbose);
			op._snapConnectingLines();
		}

		///
		///
		///
		void CFeatGenerationOp::_snapConnectingLines() const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN CL SNAPPING FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			double const distMaxClClosest = themeParameters->getValue(CL_CL_CLOSEST_MAX_DIST).toDouble();

			_snapCl2Cl( distMaxClClosest );

			_logger->log(epg::log::TITLE, "[ END CL SNAPPING FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
		}

		///
		///
		///
		void CFeatGenerationOp::DeleteConnectingLines(
			std::string const& borderCode,
			bool verbose
		) {
			CFeatGenerationOp op(borderCode, verbose);
			op._deleteConnectingLines();
		}

		///
		///
		///
		void CFeatGenerationOp::_deleteConnectingLines() const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN CL DELETING FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			_deleteClByAngleAndDistEdges();

			_deleteCLUnderThreshold();

			_logger->log(epg::log::TITLE, "[ END CL DELETING FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
		}

		///
		///
		///
		void CFeatGenerationOp::UpdateGeomConnectingLines(
			std::string const& borderCode,
			bool verbose
		) {
			CFeatGenerationOp op(borderCode, verbose);
			op._updateGeomConnectingLines();
		}

		///
		///
		///
		void CFeatGenerationOp::_updateGeomConnectingLines() const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN CL UPDATE GEOMETRY FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string idName = context->getEpgParameters().getValue(ID).toString();
			std::string geomName = context->getEpgParameters().getValue(GEOM).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			double const snapProjCl2edge = themeParameters->getValue(CL_SNAP_PROJ_CL_2_EDGE_DIST).toDouble();

			//--
			GraphType graphCL;
			_loadGraphCL(graphCL);

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

			_logger->log(epg::log::TITLE, "[ END CL UPDATE GEOMETRY FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
		}



		///
		///
		///
		void CFeatGenerationOp::ComputeCL(
			std::string const& borderCode,
			bool verbose
		) {
			CFeatGenerationOp op(borderCode, verbose);
			op._computeCL();
		}

		///
		///
		///
		void CFeatGenerationOp::_computeCL() const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN CL GENERATION FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			_generateConnectingLinesByCountry();
			_mergeConnectingLinesOnBorder();
			_snapConnectingLines();
			_deleteConnectingLines();
			_updateGeomConnectingLines();

			_logger->log(epg::log::TITLE, "[ END CL GENERATION FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

		}

		///
		///
		///
		void CFeatGenerationOp::ComputeCP(
			std::string const& borderCode,
			bool verbose
		) {
			CFeatGenerationOp op(borderCode, verbose);
			op._computeCP();
		}

		///
		///
		///
		void CFeatGenerationOp::_computeCP() const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN CP GENERATION FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			//--
			ign::geometry::MultiLineString mlsBorder = _getBorderGeom();
			epg::tools::geometry::SegmentIndexedGeometry segIndexBorder(&mlsBorder);

			//--
			_getCPfromCl();
			_getCPfromBorder(mlsBorder);
			_removeDuplicateCP();
			_getCPfromBorderUnderShoot(mlsBorder, segIndexBorder);
			_snapCPNearBy(segIndexBorder);
			_cutClByCp();

			_logger->log(epg::log::TITLE, "[ END CP GENERATION FOR " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
		}

		///
		///
		///
		void CFeatGenerationOp::_getBorderCutByAngle(
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
		void CFeatGenerationOp::_getCLfromBorder(
			ign::geometry::LineString const& lsBorder,
			ign::geometry::Geometry const& buffBorder
		) const {
			//--
			epg::Context* context = epg::ContextS::getInstance();

			//--
			epg::params::EpgParameters const& epgParams = context->getEpgParameters();
			std::string const idName = epgParams.getValue(ID).toString();
			std::string const geomName = epgParams.getValue(GEOM).toString();
			std::string const linkedFeatIdName = epgParams.getValue(LINKED_FEATURE_ID).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			double const distBuffer = themeParameters->getValue(CL_BUFFER_DIST).toDouble();
			double const thresholdNoCL = themeParameters->getValue(CL_THRESHOLD_NO_CL).toDouble();
			double const ratioInBuff = themeParameters->getValue(CL_RATIO_IN_BUFFER).toDouble();
			double const snapOnVertexBorder = themeParameters->getValue(CL_SNAP_ON_VERTEX_BORDER_DIST).toDouble();
			double angleMax = themeParameters->getValue(CL_BORDER_MAX_ANGLE).toDouble();
			angleMax = angleMax * M_PI / 180;

			//--
			std::vector<ign::feature::FeatureAttributeType> listAttrEdge = _fsEdge->getFeatureType().attributes();

			//--
			ign::feature::FeatureFilter filter("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + buffBorder.toString() + "'),3035))");
			// exclusion des arcs fusionnant déjà des arcs des deux pays
			epg::tools::FilterTools::addAndConditions(filter, " NOT ("+_isClStatement+")");
			if (_sqlFilterForCpGeneration != "")
				epg::tools::FilterTools::addAndConditions(filter, _sqlFilterForCpGeneration);

			ign::feature::FeatureIteratorPtr eit = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filter);
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filter);
			if (numFeatures == 0)
				return;
			boost::progress_display display(numFeatures, std::cout, "[ creating CL % complete ]\n");

			while (eit->hasNext())
			{
				++display;

				//--
				ign::feature::Feature fEdge = eit->next();
				ign::geometry::LineString const& lsEdge = fEdge.getGeometry().asLineString();
			
				//--
				app::geometry::tools::LineStringSplitter lsSplitter(lsEdge);
				lsSplitter.addCuttingGeometry(buffBorder);
				std::vector<ign::geometry::LineString> subEdgesBorder = lsSplitter.getSubLineStringsZ();

				//--
				std::vector<ign::geometry::LineString> vLsProjectedOnBorder;
				//pas d'intersection avec le contour du buffer
				if (subEdgesBorder.size() == 1)
				{			
					double angleEdgeBorder = _getAngleEdgeWithBorder(lsEdge, lsBorder);		
					//si l'edge est "proche" on considere qu'il est entierement dans le buffer et longe la frontiere
					if (lsEdge.distance(lsBorder) < distBuffer && (angleEdgeBorder < angleMax || angleEdgeBorder > (M_PI - angleMax) ) )
					{
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

					for (size_t i = 0; i < subEdgesBorder.size(); ++i)
					{
						ign::geometry::LineString const& lsSubEdgeCurr = subEdgesBorder[i];
						double currentSubEdgeLength = lsSubEdgeCurr.length();

						double angleSubEdgBorder = _getAngleEdgeWithBorder(lsSubEdgeCurr, lsBorder);

						int numSeg = static_cast<int>(std::floor(lsSubEdgeCurr.numSegments() / 2.));
						ign::geometry::Point interiorPointSEC = epg::tools::geometry::interpolate(lsSubEdgeCurr, numSeg, 0.5);
						bool isSubSegInBuff = false;

						if (buffBorder.contains(interiorPointSEC) && (angleSubEdgBorder < angleMax || angleSubEdgBorder > (M_PI - angleMax) ) )
						{
							isSubSegInBuff = true;
							numlastSubInBuff = i;

							lengthInBuff += currentSubEdgeLength;
							if (numfirstSubInBuff < 0)
								numfirstSubInBuff = i;
						}

						if (isSubSegInBuff || currentSubEdgeLength <= thresholdNoCL)
							lengthNearByBuff += currentSubEdgeLength;

						if ((currentSubEdgeLength > thresholdNoCL && !isSubSegInBuff) || i == subEdgesBorder.size() - 1)
						{
							if (lengthInBuff > ratioInBuff * lengthNearByBuff)
							{
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

				for (size_t i = 0; i < vLsProjectedOnBorder.size(); ++i)
				{
					//create CL
					//generation de l'id CL
					ign::feature::Feature featCL = _fsCL->newFeature();
					vLsProjectedOnBorder[i].setFillZ(0);
					featCL.setGeometry(vLsProjectedOnBorder[i]);
					featCL.setAttribute(linkedFeatIdName, ign::data::String(fEdge.getId()));
					for (std::vector<ign::feature::FeatureAttributeType>::iterator vit = listAttrEdge.begin(); vit != listAttrEdge.end(); ++vit)
					{
						std::string attrName = vit->getName();
						if (attrName == geomName || attrName == idName || !_fsCL->getFeatureType().hasAttribute(attrName) )
							continue;
						featCL.setAttribute(attrName, fEdge.getAttribute(attrName));
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
		ign::geometry::LineString CFeatGenerationOp::_getLinestring(
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
		void CFeatGenerationOp::_getCPfromCl() const {
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

			//--
			ign::feature::FeatureFilter filterEdge("ST_INTERSECTS(" + geomName + ", (SELECT ST_Union(array(SELECT "+geomName+" FROM "+ _fsEdge->getTableName()+" WHERE "+ _isClStatement + "))))");
			if (_sqlFilterForCpGeneration != "")
				epg::tools::FilterTools::addAndConditions(filterEdge, _sqlFilterForCpGeneration);
			epg::tools::FilterTools::addAndConditions(filterEdge, "((" + _mIsCountryStatement.begin()->second + ") OR (" + _mIsCountryStatement.rbegin()->second + "))");

			//--
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterEdge);
			boost::progress_display display(numFeatures, std::cout, "[ creating CP from CL % complete ]\n");

			ign::feature::FeatureIteratorPtr itEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterEdge);

			while (itEdge->hasNext())
			{
				++display;

				ign::feature::Feature fEdge = itEdge->next();
				ign::geometry::LineString const& edgeGeom = fEdge.getGeometry().asLineString();

				std::vector<ign::geometry::LineString> vCl = _getIntersectingCls(edgeGeom);
				ign::geometry::MultiLineString mlsCl(vCl);

				ign::geometry::GeometryPtr intersectionGeomPtr(edgeGeom.Intersection(mlsCl));

				if (intersectionGeomPtr->isNull() || intersectionGeomPtr->isEmpty())
					continue;

				_recordCp(*intersectionGeomPtr, fEdge, _tagFromCl);
			}
		}

		///
		///
		///
		ign::geometry::MultiLineString CFeatGenerationOp::_getBorderGeom() const {
			ign::geometry::MultiLineString mlsBorder;

			//--
			epg::Context* context = epg::ContextS::getInstance();

			//--
			epg::params::EpgParameters const& epgParams = context->getEpgParameters();
			std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();
			std::string const boundaryTypeName = epgParams.getValue(BOUNDARY_TYPE).toString();
			std::string const typeInternationalBoundary = epgParams.getValue(TYPE_INTERNATIONAL_BOUNDARY).toString();
			std::string const typeCoastline = epgParams.getValue(TYPE_COASTLINE).toString();

			//--
			ign::feature::FeatureIteratorPtr itBoundary = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsBoundary, ign::feature::FeatureFilter(countryCodeName + " = '" + _borderCode + "'"));
			while (itBoundary->hasNext()) {

				ign::feature::Feature fBoundary = itBoundary->next();
				std::string const& boundaryType = fBoundary.getAttribute(boundaryTypeName).toString();
				ign::geometry::LineString const& lsBoundary = fBoundary.getGeometry().asLineString();
				
				if (boundaryType != typeInternationalBoundary && boundaryType.find(typeCoastline) == -1)
					continue;

				mlsBorder.addGeometry(lsBoundary);
			}

			return mlsBorder;
		}

		///
		///
		///
		bool CFeatGenerationOp::_isDangle(
			ign::feature::Feature const& fEdge,
			CFeatGenerationOp::ENDING ending
		) const {
			ign::geometry::LineString const& edgeGeom = fEdge.getGeometry().asLineString();
			ign::geometry::Point const& endingPoint = ending == START ? edgeGeom.startPoint() : edgeGeom.endPoint();

			//--
			ign::feature::FeatureFilter filterEdge;
			if (_sqlFilterForCpGeneration != "")
				filterEdge.setPropertyConditions(_sqlFilterForCpGeneration);
			filterEdge.setExtent(endingPoint.getEnvelope().expandBy(1));

			//--
			ign::feature::FeatureIteratorPtr itEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterEdge);
			while (itEdge->hasNext()) {
				ign::feature::Feature fEdgeAround = itEdge->next();
				ign::geometry::LineString const& edgeAroundGeom = fEdgeAround.getGeometry().asLineString();

				if (fEdgeAround.getId() == fEdge.getId())
					continue;

				double distStart = edgeAroundGeom.startPoint().distance(endingPoint);
				double distEnd = edgeAroundGeom.endPoint().distance(endingPoint);

				if( std::min(distStart, distEnd) < 1e-5 )
					return false;
			}
			return true;
		}

		///
		///
		///
		void CFeatGenerationOp::_getCPfromBorderUnderShoot(
			ign::geometry::Geometry const& borderGeom,
			epg::tools::geometry::SegmentIndexedGeometry const& segIndexBorder
		) const {
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			double const distUnderShoot = themeParameters->getValue(CP_UNDERSHOOT_DIST).toDouble();
			double const maxLengthOverShoot = themeParameters->getValue(CP_OVERSHOOT_MAX_LENTH).toDouble();

			//--
			ign::geometry::algorithm::BufferOpGeos buffOp;
			ign::geometry::GeometryPtr buffBorder(buffOp.buffer(borderGeom, distUnderShoot, 0, ign::geometry::algorithm::BufferOpGeos::CAP_FLAT));

			std::ostringstream oss;
			oss << "(";
			oss << "SELECT ST_Intersects(ST_StartPoint(" << geomName << "), p.buff) OR ST_Intersects(ST_EndPoint(" << geomName << "), p.buff)";
			oss << "FROM (";
			oss << "SELECT ST_GeomFromText('"<< buffBorder->toString() << "', 3035) AS buff";
			oss << ") p";
			oss << ")";

			//--
			ign::feature::FeatureFilter filterEdge(oss.str());
			if (_sqlFilterForCpGeneration != "")
				epg::tools::FilterTools::addAndConditions(filterEdge, _sqlFilterForCpGeneration);
			epg::tools::FilterTools::addAndConditions(filterEdge, "((" + _mIsCountryStatement.begin()->second + ") OR (" + _mIsCountryStatement.rbegin()->second + "))");

			//--
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterEdge);
			boost::progress_display display(numFeatures, std::cout, "[ computing CP by resolving undershoot % complete ]\n");

			ign::feature::FeatureIteratorPtr itEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterEdge);
			while (itEdge->hasNext())
			{
				++display;

				ign::feature::Feature fEdge = itEdge->next();
				ign::geometry::LineString const& lsEdge = fEdge.getGeometry().asLineString();

				// est ce qu'un overshoot précède l'undershoot ?
				std::vector<ign::geometry::LineString> vBorderCuttingParts;
				segIndexBorder.getSegments( lsEdge.getEnvelope(), vBorderCuttingParts );

				app::geometry::tools::LineStringSplitter splitter(lsEdge);
				ign::geometry::MultiLineString mlsCuttingParts(vBorderCuttingParts);
				splitter.addCuttingGeometry(mlsCuttingParts);
				std::vector< ign::geometry::LineString > vEdgePart = splitter.getSubLineStrings();

				std::map<double, ENDING> mEnding;

				double distanceStart = lsEdge.startPoint().distance(borderGeom);
				if( distanceStart < distUnderShoot && _isDangle(fEdge, START) && vEdgePart.front().length() > maxLengthOverShoot )
					mEnding.insert(std::make_pair(distanceStart, START));

				double distanceEnd = lsEdge.endPoint().distance(borderGeom);
				if( distanceEnd < distUnderShoot && _isDangle(fEdge, END) && vEdgePart.back().length() > maxLengthOverShoot )
					mEnding.insert(std::make_pair(distanceEnd, END));

				//DEBUG
				//VOIR les cas:
				// - undershoot aux deux extremites d'un petit arc
				// - intersection a une meme extremite proche frontiere (dist < distUnderShoot)
				// - intersection et undershoot sur petit arc

				//--
				for( std::map<double, ENDING>::const_iterator mit = mEnding.begin() ; mit != mEnding.end() ; ++mit )
				{
					ign::geometry::Point const& endingPoint = mit->second == START ? lsEdge.startPoint() : lsEdge.endPoint();
					ign::geometry::Point const& previousPoint = mit->second == START ? lsEdge.pointN(1) : lsEdge.pointN(lsEdge.numPoints()-2);

					ign::geometry::Point projPt;
					double distMin = std::numeric_limits<double>::max();

					//proj axiale
					std::vector<ign::geometry::LineString> vBorderSegments;
					segIndexBorder.getSegments( endingPoint.getEnvelope().expandBy(distUnderShoot*1.5), vBorderSegments );
					if (vBorderSegments.empty())
						continue;
					ign::geometry::MultiLineString mlsBorderSegments(vBorderSegments);
					std::vector< ign::geometry::Point > vPtIntersect = epg::tools::geometry::LineIntersector::compute(previousPoint, endingPoint, mlsBorderSegments);
					ign::math::Line2d line(previousPoint.toVec2d(), endingPoint.toVec2d());
					
					for (std::vector< ign::geometry::Point >::iterator vit = vPtIntersect.begin(); vit != vPtIntersect.end(); ++vit) {
						double abs = line.project(vit->toVec2d());
						if (abs < 1) continue; /* contient les cas abs < 0  et 0 < abs < 1 */

						double dist = endingPoint.distance(*vit);
						if (dist < distMin) {
							projPt = *vit;
							distMin = dist;
						}
					}

					//proj ortho
					ign::geometry::Point projPt2 = epg::tools::geometry::project(mlsBorderSegments, endingPoint, 0);
					//on verifie que la projection n'est pas un rebroussement
					ign::geometry::Point projPt2Reverse = epg::tools::geometry::project(lsEdge, projPt2, 0);
					double angle = _getAngle(endingPoint, projPt2, endingPoint, projPt2Reverse);
					if( angle > 1.3963 ) {/*80 degrees*/
						double distance2 = endingPoint.distance(projPt2);
						if (distance2 < distMin/3 && distance2 < distUnderShoot/3) {
							distMin = distance2;
							projPt = projPt2;
						}
					}

					//--
					if (distMin > distUnderShoot)
						continue;

					//DEBUG voir pour garder le Z de l'extremité projetée
					projPt.setZ(0);
					_recordCp(projPt, fEdge);
				}
			}
		}

		///
		///
		///
		void CFeatGenerationOp::_getCPfromBorder(
			ign::geometry::Geometry const& borderGeom
		) const {
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
			std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();

			//--
			ign::feature::FeatureFilter filterEdge("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + borderGeom.toString() + "'),3035))");
			if (_sqlFilterForCpGeneration != "")
				epg::tools::FilterTools::addAndConditions(filterEdge, _sqlFilterForCpGeneration);
			epg::tools::FilterTools::addAndConditions(filterEdge, "((" + _mIsCountryStatement.begin()->second + ") OR (" + _mIsCountryStatement.rbegin()->second + "))");
			
			//--
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterEdge);
			boost::progress_display display(numFeatures, std::cout, "[ creating CP from border % complete ]\n");

			ign::feature::FeatureIteratorPtr itEdge = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterEdge);
			while (itEdge->hasNext())
			{
				++display;

				ign::feature::Feature fEdge = itEdge->next();
				ign::geometry::LineString const& edgeGeom = fEdge.getGeometry().asLineString();
				
				//TODO : optimisation possible en utilisant segIndexBorder
				ign::geometry::GeometryPtr intersectionGeomPtr(edgeGeom.Intersection(borderGeom));

				_recordCp(*intersectionGeomPtr, fEdge);
			}
		}

		///
		///
		///
		void CFeatGenerationOp::_removeDuplicateCP() const {
			std::list<std::string> lCpToDelete;

			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
			
			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getValue(W_TAG_NAME).toString();

			//--
			ign::feature::FeatureFilter filterCp(wTagName + " IS DISTINCT FROM '" + _tagFromCl +"'");

			//--
			ign::feature::FeatureIteratorPtr itCp = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCP, filterCp);

			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCP, filterCp);
			boost::progress_display display(numFeatures, std::cout, "[ deleting duplicate CP % complete ]\n");

			while (itCp->hasNext())
			{
				++display;

				ign::feature::Feature fCp = itCp->next();
				std::string cpId = fCp.getId();

				std::pair<bool, ign::feature::Feature> foundFromCLCandidate = _hasDuplicateCandidate(fCp, true);

				if( !foundFromCLCandidate.first )
					continue;

				std::pair<bool, ign::feature::Feature> foundFromBorderCandidate = _hasDuplicateCandidate(foundFromCLCandidate.second, false);

				if( !foundFromBorderCandidate.first || foundFromBorderCandidate.second.getId() != cpId ) {
					_logger->log(epg::log::WARN, "candidate association failure for CP : "+cpId);
					continue;
				}

				//DEBUG
				{
					ign::feature::Feature fDuplicate;
					fDuplicate.setGeometry(ign::geometry::LineString(fCp.getGeometry().asPoint(), foundFromCLCandidate.second.getGeometry().asPoint()));
					_shapeLogger->writeFeature("debug_duplicate", fDuplicate);
				}

				lCpToDelete.push_back(cpId);
			}

			for( std::list<std::string>::const_iterator lit = lCpToDelete.begin() ; lit != lCpToDelete.end() ; ++lit )
				_fsCP->deleteFeature(*lit);
		}

		///
		///
		///
		std::pair<bool, ign::feature::Feature> CFeatGenerationOp::_hasDuplicateCandidate(
			ign::feature::Feature const& fCp,
			bool fromCl
		) const {
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
			std::string const idName = context->getEpgParameters().getValue(ID).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			double const pairingDist = themeParameters->getValue(CP_FROM_CL_FROM_BORDER_PAIRING_DIST).toDouble();
			std::string const wTagName = themeParameters->getValue(W_TAG_NAME).toString();

			//--
			ign::geometry::Point const& cpGeom = fCp.getGeometry().asPoint();
			std::string const linkedEdgeId = fCp.getAttribute(linkedFeatIdName).toString();
			
			//--
			ign::feature::FeatureFilter filterCp;
			filterCp.setExtent(cpGeom.getEnvelope().expandBy(pairingDist));
			epg::tools::FilterTools::addAndConditions(filterCp, linkedFeatIdName + " = '" + linkedEdgeId + "'");
			epg::tools::FilterTools::addAndConditions(filterCp, idName + " != '" + fCp.getId() + "'");
			if(fromCl)
				epg::tools::FilterTools::addAndConditions(filterCp, wTagName + " = '" + _tagFromCl + "'");

			//--
			ign::feature::FeatureIteratorPtr itCp = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCP, filterCp);

			double minDist = std::numeric_limits<double>::infinity();
			ign::feature::Feature minFeat;
			while (itCp->hasNext())
			{
				ign::feature::Feature fCandidate = itCp->next();
				ign::geometry::Point const& candidateGeom = fCandidate.getGeometry().asPoint();

				double distance = candidateGeom.distance(cpGeom);

				if( distance < minDist ) {
					minDist = distance;
					minFeat = fCandidate;
				}
			}

			if( minDist > pairingDist )
				return std::make_pair(false, ign::feature::Feature());

			return std::make_pair(true, minFeat);
		}

		///
		///
		///
		void CFeatGenerationOp::_recordCp(
			ign::geometry::Geometry const& cpGeom,
			ign::feature::Feature const& linkedEdgeFeat,
			std::string tag
		) const {
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
			std::string const idName = context->getEpgParameters().getValue(ID).toString();
			std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getValue(W_TAG_NAME).toString();

			//--
			std::vector<ign::feature::FeatureAttributeType> const& listAttrEdge = _fsEdge->getFeatureType().attributes();

			//--
			ign::feature::Feature fCF = _fsCP->newFeature();
			fCF.setAttribute(linkedFeatIdName, ign::data::String(linkedEdgeFeat.getId()));
			for (std::vector<ign::feature::FeatureAttributeType>::const_iterator vit = listAttrEdge.begin(); vit != listAttrEdge.end(); ++vit)
			{
				std::string attrName = vit->getName();
				if (attrName == geomName || attrName == idName || !_fsCP->getFeatureType().hasAttribute(attrName))
					continue;
				fCF.setAttribute(attrName, linkedEdgeFeat.getAttribute(attrName));
			}
			if( tag != "" ) 
				fCF.setAttribute(wTagName, ign::data::String(tag));

			if (cpGeom.isPoint())
			{
				fCF.setGeometry(cpGeom.asPoint());
				std::string idCP = _idGeneratorCP->next();
				_fsCP->createFeature(fCF, idCP);
			}
			else if (cpGeom.isGeometryCollection())
			{
				ign::geometry::GeometryCollection const& geomCollect = cpGeom.asGeometryCollection();
				for (size_t i = 0; i < geomCollect.numGeometries(); ++i)
				{
					if (geomCollect.geometryN(i).isPoint())
					{
						ign::geometry::Point const& ptIntersect = geomCollect.geometryN(i).asPoint();		
						fCF.setGeometry(ptIntersect);
						std::string idCP = _idGeneratorCP->next();
						_fsCP->createFeature(fCF, idCP);
					}
				}
			}
		}

		///
		///
		///
		double CFeatGenerationOp::_getAngleEdgeWithBorder(
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
		double CFeatGenerationOp::_getAngle(
			ign::geometry::Point const& ptSource1,
			ign::geometry::Point const& ptTarget1,
			ign::geometry::Point const& ptSource2,
			ign::geometry::Point const& ptTarget2
		) const {
			if( ptSource1.distance(ptTarget1) < 1e-5 || ptSource2.distance(ptTarget2) < 1e-5 )
				return std::numeric_limits<double>::infinity();
			ign::math::Vec2d v1(ptTarget1.x() - ptSource1.x(), ptTarget1.y() - ptSource1.y());
			ign::math::Vec2d v2(ptTarget2.x() - ptSource2.x(), ptTarget2.y() - ptSource2.y());

			return epg::tools::geometry::angle(v1, v2);
		}

		///
		///
		///
		ign::geometry::LineString CFeatGenerationOp::_getGeomCL(
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
		void CFeatGenerationOp::_snapCPNearBy(
			epg::tools::geometry::SegmentIndexedGeometry const& segIndexBorder
		) const {
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
			std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

			//--
			ign::feature::FeatureFilter filterCP("((" + _mIsCountryStatement.begin()->second + ") OR (" + _mIsCountryStatement.rbegin()->second + "))");

			//--
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCP, filterCP);
			boost::progress_display display(numFeatures, std::cout, "[ snapping CP near by % complete ]\n");

			ign::feature::FeatureIteratorPtr itCP = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCP, filterCP);

			std::set<std::string> sCP2Snap;
			std::string separator = "#";
			while (itCP->hasNext())
			{
				++display;
				
				ign::feature::Feature fCPCurr = itCP->next();
				std::string idCP = fCPCurr.getId();

				//DEBUG
				// if ( fCPCurr.getGeometry().distance(ign::geometry::Point(3803179.324, 3104038.058)) < 0.5 ) {
				// 	 bool test = true;
				// }
				// if ( fCPCurr.getGeometry().distance(ign::geometry::Point(3803180.104, 3104037.315)) < 0.5 ) {
				// 	 bool test = true;
				// }
				// if ( fCPCurr.getGeometry().distance(ign::geometry::Point(3803183.32, 3104039.89)) < 1 ) {
				// 	 bool test = true;
				// }

				if (sCP2Snap.find(idCP) != sCP2Snap.end())
					continue;

				// on recupere recursivement les CP proches
				//SUPPRIMER LE retour BOOL (tester si mCPNear.size() > 1)
				std::map<std::string, ign::feature::Feature> mCPNear;
				bool hasNearestCP = _getNearestCP(fCPCurr, mCPNear);

				//DEBUG BEGIN
				//que ce passe t il en cas de merging d'un CP sur CL avec CP sur border ?
				//dans quels cas et comment les CP sur CL sont ils mergés ?
				size_t countCpOnCl = 0;
				size_t countCpOnBorder = 0;
				ign::geometry::MultiPoint debugMultiPtCP;
				for ( std::map<std::string, ign::feature::Feature>::const_iterator debut_it = mCPNear.begin() ; debut_it != mCPNear.end() ; ++debut_it )
				{
					debugMultiPtCP.addGeometry(debut_it->second.getGeometry().asPoint());
					std::vector<ign::geometry::LineString> vCls = _getIntersectingCls(debut_it->second.getGeometry());
					if ( vCls.size() > 0 ) {
						++countCpOnCl;
					} else {
						++countCpOnBorder;
					}
				}
				{
					ign::feature::Feature feat;
					feat.setGeometry(debugMultiPtCP);
					if(countCpOnCl > 1)
						_shapeLogger->writeFeature("debug_merging_of_2_cp_on_cl", feat);
					if(countCpOnCl > 0 && countCpOnBorder > 0)
						_shapeLogger->writeFeature("debug_merging_of_cp_on_cl_and_cp_border", feat);
				}
				//DEBUG END

				std::list<std::string> lCp2Delete;

				if (!hasNearestCP) {
					_fsCP->deleteFeature(fCPCurr.getId());
					continue;
				}

				std::set<std::string> s1, s2;
				std::map<std::string, ign::feature::Feature> mLinkedEdgeFeature;
				for(std::map<std::string, ign::feature::Feature>::iterator mit = mCPNear.begin(); mit != mCPNear.end(); ++mit)
				{
					std::string cpCountryCode = mit->second.getAttribute(countryCodeName).toString();
					if(cpCountryCode.find(_vCountry.front()) != std::string::npos) 
						s1.insert(mit->first);
					else 
						s2.insert(mit->first);

					sCP2Snap.insert(mit->first);

					mLinkedEdgeFeature.insert(std::make_pair(mit->first, ign::feature::Feature()));

					_fsEdge->getFeatureById(mit->second.getAttribute(linkedFeatIdName).toString(), mLinkedEdgeFeature[mit->first]);
				}

				// map pour optimisation
				std::map<std::string, std::map<std::string, std::pair<bool, double>>> mmAreMergeable;
				for (std::set<std::string>::const_iterator sit1 = s1.begin() ; sit1 != s1.end() ; ++sit1)
				{
					std::map<std::string, std::pair<bool, double>> mAreMergeable;
					for (std::set<std::string>::const_iterator sit2 = s2.begin() ; sit2 != s2.end() ; ++sit2)
					{
						double distance = mCPNear[*sit1].getGeometry().asPoint().distance(mCPNear[*sit2].getGeometry().asPoint());
						bool areMergeable = _areMergeable(mLinkedEdgeFeature[*sit1], mLinkedEdgeFeature[*sit2], distance);
						mAreMergeable.insert(std::make_pair(*sit2, std::make_pair(areMergeable, distance)));
					}
					mmAreMergeable.insert(std::make_pair(*sit1, mAreMergeable));
				}

				std::map<std::string, std::string> m1;
				for (std::set<std::string>::const_iterator sit1 = s1.begin() ; sit1 != s1.end() ; ++sit1)
				{
					std::string samicopain = "";
					double distanceMax = std::numeric_limits<double>::max();
					for (std::set<std::string>::const_iterator sit2 = s2.begin() ; sit2 != s2.end() ; ++sit2)
					{
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
				for (std::set<std::string>::const_iterator sit2 = s2.begin() ; sit2 != s2.end() ; ++sit2)
				{
					std::string samicopain = "";
					double distanceMax = std::numeric_limits<double>::max();
					for (std::set<std::string>::const_iterator sit1 = s1.begin() ; sit1 != s1.end() ; ++sit1)
					{
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
				for (std::map<std::string, std::string>::const_iterator mit1 = m1.begin() , next_mit1 = mit1 ; mit1 != m1.end() ; mit1 = next_mit1)
				{
					++next_mit1;
					for (std::map<std::string, std::string>::const_iterator mit2 = m2.begin() ; mit2 != m2.end() ; ++mit2)
					{
						if( mit1->second == mit2->first && mit2->second == mit1->first)
						{
							mapCpGroup.insert(value_type(mit1->first, ++group));
							mapCpGroup.insert(value_type(mit2->first, group));
							m1.erase(mit1);
							m2.erase(mit2);
							break;
						}
					}
				}

				for (std::map<std::string, std::string>::const_iterator mit1 = m1.begin() , next_mit1 = mit1 ; mit1 != m1.end() ; mit1 = next_mit1)
				{
					++next_mit1;
					std::map<std::string, std::string>::const_iterator mit2 = m2.find(mit1->second);
					if( mit2 == m2.end()) continue; /*cp de m2 déjà affecté à un groupe*/

					mapCpGroup.insert(value_type(mit1->first,++group));
					mapCpGroup.insert(value_type(mit2->first,group));
					m1.erase(mit1);
					m2.erase(mit2);
				}

				for (std::map<std::string, std::string>::const_iterator mit2 = m2.begin() , next_mit2 = mit2 ; mit2 != m2.end() ; mit2 = next_mit2)
				{
					++next_mit2;
					std::map<std::string, std::string>::const_iterator mit1 = m1.find(mit2->second);
					if( mit1 == m1.end()) continue; /*cp de m1 déjà affecté à un groupe*/

					mapCpGroup.insert(value_type(mit2->first,++group));
					mapCpGroup.insert(value_type(mit1->first,group));
					m2.erase(mit2);
					m1.erase(mit1);
				}

				// il reste les cp dont le cp associé est déjà relié à un groupe
				for (std::map<std::string, std::string>::const_iterator mit2 = m2.begin() ; mit2 != m2.end() ; ++mit2)
				{
					auto l_mit = mapCpGroup.left.find(mit2->second);
					IGN_ASSERT(l_mit != mapCpGroup.left.end());
					mapCpGroup.insert(value_type(mit2->first, l_mit->second));
				}
				for (std::map<std::string, std::string>::const_iterator mit1 = m1.begin() ; mit1 != m1.end() ; ++mit1)
				{
					auto l_mit = mapCpGroup.left.find(mit1->second);
					IGN_ASSERT(l_mit != mapCpGroup.left.end());
					mapCpGroup.insert(value_type(mit1->first, l_mit->second));
				}

				//TODO ajouter une fonction pour separer un group en plusieur s'il contient plusieur CP from_cl non superposé
				// faire un test pour voir si ça existe (3073514,03  3845868,64)

				for (size_t i = 1 ; i <= group ; ++i)
				{
					ign::geometry::Point newGeom;

					//--
					auto range = mapCpGroup.right.equal_range(i);
					std::list<std::string> lCp;
					for (auto r_it = range.first; r_it != range.second; ++r_it)
						lCp.push_back(r_it->second);
					
					//construire une map de CP mCp a partir de mCPNear et lCp
					std::list<ign::feature::Feature> lCPfromCL = _pickUpCPfromCL(lCp, mCPNear);
					if( lCPfromCL.size() > 0 ) {
						newGeom = _getBestCpFromClGeom(lCPfromCL);
					} else {
						ign::geometry::MultiPoint mlpCP;
						for (std::list<std::string>::const_iterator lit = lCp.begin() ; lit != lCp.end(); ++lit)
							_addNoDuplicate(mlpCP, mCPNear[*lit].getGeometry().asPoint());

						ign::geometry::Point centroidPoint = mlpCP.getCentroid();
						
						//--
						std::vector<ign::geometry::LineString> vBorderSegments;
						segIndexBorder.getSegments( mlpCP.getEnvelope(), vBorderSegments );
						ign::geometry::MultiLineString mlsBorderSegments(vBorderSegments);

						newGeom = epg::tools::geometry::project(mlsBorderSegments, centroidPoint, 0 /*border vertex snap dist*/);

						//DEBUG voir pour calculer le Z
						newGeom.setZ(0);
					}

					//--
					for (std::list<std::string>::const_iterator lit = lCp.begin() ; lit != lCp.end(); ++lit)
					{
						mCPNear[*lit].setGeometry(newGeom);
						_fsCP->modifyFeature(mCPNear[*lit]);
					}
				}

				for (std::list<std::string>::const_iterator lit = lCp2Delete.begin() ; lit != lCp2Delete.end() ; ++lit)
					_fsCP->deleteFeature(*lit);
			}
		}

		///
		///
		///
		bool CFeatGenerationOp::_addNoDuplicate(
			ign::geometry::MultiPoint & mpt,
			ign::geometry::Point const& pt,
			double precision
		) const {
			if( !mpt.isEmpty() && mpt.distance(pt) < precision )
				return false;

			mpt.addGeometry(pt);
			return true;
		}

		///
		///
		///
		ign::geometry::Point CFeatGenerationOp::_getBestCpFromClGeom(
			std::list<ign::feature::Feature> const& lCpFromCl
		) const {
			//--
			if ( lCpFromCl.size() == 1 )
				return lCpFromCl.begin()->getGeometry().asPoint();

			//--
			// si une CL commune : calcul du point cible = moyenne des abscisses
			ign::geometry::MultiPoint mlptNoDuplicate;
			_addNoDuplicate(mlptNoDuplicate, lCpFromCl.begin()->getGeometry().asPoint());
			std::map<std::string, ign::feature::Feature> mConnectedClRef = _getConnectecCl(lCpFromCl.begin()->getGeometry().asPoint());
			for ( std::list<ign::feature::Feature>::const_iterator lit = std::next(lCpFromCl.begin()) ; lit != lCpFromCl.end() ; ++lit )
			{
				_addNoDuplicate(mlptNoDuplicate, lit->getGeometry().asPoint());
				std::map<std::string, ign::feature::Feature> mConnectedCl = _getConnectecCl(lit->getGeometry().asPoint());
				for ( std::map<std::string, ign::feature::Feature>::const_iterator mit = mConnectedClRef.begin() ; mit != mConnectedClRef.end() ; )
				{
					if( mConnectedCl.find(mit->first) == mConnectedCl.end() ) {
						mit = mConnectedClRef.erase(mit);
					} else {
						++mit;
					}
				}
				if( mConnectedClRef.size() == 0 )
					break;
			}

			if( mConnectedClRef.size() == 1 && mlptNoDuplicate.numGeometries() > 1 )
			{
				ign::geometry::LineString const& clGeom = mConnectedClRef.begin()->second.getGeometry().asLineString();
				geometry::tools::LengthIndexedLineString indexedClGeom(clGeom);
				double meanAbs = 0;
				for ( size_t i = 0 ; i < mlptNoDuplicate.numGeometries() ; ++i )
						meanAbs += indexedClGeom.project(mlptNoDuplicate.pointN(i));

				meanAbs /= mlptNoDuplicate.numGeometries();
				return indexedClGeom.locateAlong(meanAbs);
			}

			//--
			ign::geometry::Point centroid = mlptNoDuplicate.getCentroid();

			//--
			double distanceMin = std::numeric_limits<double>::infinity();
			ign::geometry::Point ptMin;
			for ( std::list<ign::feature::Feature>::const_iterator lit = lCpFromCl.begin() ; lit != lCpFromCl.end() ; ++lit )
			{
				double distance = centroid.distance(lit->getGeometry().asPoint());
				if ( distance < distanceMin ) {
					distanceMin = distance;
					ptMin = lit->getGeometry().asPoint();
				}
			}
			return ptMin;
		}

		///
		///
		///
		std::map<std::string, ign::feature::Feature> CFeatGenerationOp::_getConnectecCl(
			ign::geometry::Point const& pt,
			double precision
		) const {
			std::map<std::string, ign::feature::Feature> mConnectedCl;

			//--
			ign::feature::FeatureFilter filterCl(_isClStatement);
			filterCl.setExtent(pt.getEnvelope().expandBy(precision));

			//--
			ign::feature::FeatureIteratorPtr itCl = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterCl);
			while (itCl->hasNext())
			{
				ign::feature::Feature fCl = itCl->next();

				if( fCl.getGeometry().distance(pt) < precision )
					mConnectedCl.insert(std::make_pair(fCl.getId(), fCl));
			}
			return mConnectedCl;
		}

		///
		///
		///
		std::list<ign::feature::Feature> CFeatGenerationOp::_pickUpCPfromCL(
			std::list<std::string> const& lCp,
			std::map<std::string, ign::feature::Feature> const& mCPNear
		) const {
			std::list<ign::feature::Feature> lCpFromCl;

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getValue(W_TAG_NAME).toString();

			for( std::list<std::string>::const_iterator lit = lCp.begin() ; lit != lCp.end() ; ++lit )
			{
				std::map<std::string, ign::feature::Feature>::const_iterator mit = mCPNear.find(*lit);
				if( mit != mCPNear.end() ) {
					if( mit->second.getAttribute(wTagName).toString() == _tagFromCl )
						lCpFromCl.push_back(mit->second);
				}
			}
			return lCpFromCl;
		}

		///
		///
		///
		bool CFeatGenerationOp::_areMergeable(
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
		bool CFeatGenerationOp::_areDistanceTypeCompatible(
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

			bool isWalkwayOrTractor1 = _sFormOfWayException.find(type1) != _sFormOfWayException.end();
			bool isWalkwayOrTractor2 = _sFormOfWayException.find(type2) != _sFormOfWayException.end();

			return isWalkwayOrTractor1 && isWalkwayOrTractor2 ? distance < distMergeTractorCP : distance < distMergeCP;
		}

		///
		///
		///
		bool CFeatGenerationOp::_areCollinear(
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
		std::vector<ign::geometry::LineString> CFeatGenerationOp::_getIntersectingCls(
			ign::geometry::Geometry const& geom
		) const {
			std::vector<ign::geometry::LineString> vCl;

			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
			std::string geomName = context->getEpgParameters().getValue(GEOM).toString();
			
			//--
			ign::feature::FeatureFilter filterCl("ST_INTERSECTS(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + geom.toString() + "'),3035))");
			epg::tools::FilterTools::addAndConditions(filterCl, _isClStatement);

			//--
			ign::feature::FeatureIteratorPtr itCl = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterCl);
			while (itCl->hasNext())
			{
				vCl.push_back(itCl->next().getGeometry().asLineString());
			}

			return vCl;
		}

		///
		///
		///
		void CFeatGenerationOp::_cutClByCp() const 
		{
			std::list<std::string> lCl2Delete;

			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
			std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const wTagName = themeParameters->getValue(W_TAG_NAME).toString();

			//--
			ign::feature::FeatureFilter filterCL(_isClStatement);
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filterCL);
			boost::progress_display display(numFeatures, std::cout, "[ cutting CL by CP % complete ]\n");

			ign::feature::FeatureIteratorPtr itCL = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filterCL);
			while (itCL->hasNext())
			{
				++display;

				ign::feature::Feature fCL = itCL->next();
				ign::geometry::LineString const lsCl = fCL.getGeometry().asLineString();
				
				//on cherche les CP a proximité
				ign::geometry::MultiPoint mpCp;
				ign::feature::FeatureFilter filterCp("ST_DISTANCE(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + lsCl.toString() + "'),3035)) < 0.1");
				epg::tools::FilterTools::addAndConditions(filterCp, wTagName + " = '" + _tagFromCl +"'");
				ign::feature::FeatureIteratorPtr itCL = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCP, filterCp);
				while (itCL->hasNext()) {
					mpCp.addGeometry(itCL->next().getGeometry().asPoint());
				}

				if( mpCp.numGeometries() == 0 )
					continue;

				app::geometry::tools::LineStringSplitter clSplitter(lsCl, 0.1); //on snappe à 10cm si il y a plusieurs coupures
				clSplitter.addCuttingGeometry(mpCp);

				std::vector<ign::geometry::LineString> subCl = clSplitter.getSubLineStringsZ();

				if (subCl.size() == 1)
					continue;

				for( size_t i = 0 ; i < subCl.size() ; ++i )
				{
					if ( i > 0 )
						subCl[i].startPoint() = _getClosestGeometry( mpCp, lsCl, subCl[i].startPoint() );
					if ( i < subCl.size()-1 )
						subCl[i].endPoint() = _getClosestGeometry( mpCp, lsCl, subCl[i].endPoint() );
				}

				lCl2Delete.push_back(fCL.getId());

				for (size_t i = 0; i < subCl.size(); ++i)
				{
					fCL.setGeometry(subCl[i]);
					_fsEdge->createFeature(fCL);

					//--
					_shapeLogger->writeFeature("edgeClCutByCp", fCL);
					_logger->log(epg::log::INFO, "edge create by _cutClByCp: " + fCL.getId());
				}
			}

			for (std::list<std::string>::iterator lit = lCl2Delete.begin(); lit != lCl2Delete.end(); ++lit)
				_fsEdge->deleteFeature(*lit);
		}

		///
		///
		///
		ign::geometry::Point CFeatGenerationOp::_getClosestGeometry(
			ign::geometry::MultiPoint const& mlp,
			ign::geometry::LineString const& ls,
			ign::geometry::Point const& pt
		) const {
			double distMin = std::numeric_limits<double>::infinity();
			int indexMin = -1;
			for( size_t i = 0 ; i < mlp.numGeometries() ; ++i ) {
				ign::geometry::Point projPt = epg::tools::geometry::project(ls, mlp.pointN(i), 0);
				double dist = projPt.distance(pt);
				if( dist < distMin ) {
					distMin = dist;
					indexMin = i;
				}
			}
			return mlp.pointN(indexMin);
		}

		///
		///
		///
		// void CFeatGenerationOp::_cutClByCp(
		// 	std::map<std::string,std::pair<ign::feature::Feature, ign::geometry::MultiPoint>> const& mClSplittedByCp
		// ) const {
		// 	std::set<std::string> sCl2delete;
		// 	boost::progress_display display(mClSplittedByCp.size(), std::cout, "[ cutting cl by cp % complete ]\n");
		// 	for (std::map<std::string, std::pair<ign::feature::Feature, ign::geometry::MultiPoint>>::const_iterator mit = mClSplittedByCp.begin();
		// 		mit != mClSplittedByCp.end(); ++mit) {
		// 		++display;
		// 		ign::feature::Feature const& fEdgDoubl2Cut = mit->second.first;
		// 		ign::geometry::LineString const& lsCl = fEdgDoubl2Cut.getGeometry().asLineString();
		// 		app::geometry::tools::LineStringSplitter lsSplitterClSnappedOn(lsCl, 0.1); //on snappe à 10cm si il y a plusieurs coupures
		// 		//boucle sur mpt
		// 		ign::geometry::MultiPoint const& mpt = mit->second.second;
		// 		for (size_t i = 0; i < mpt.numGeometries(); ++i)
		// 			lsSplitterClSnappedOn.addCuttingGeometry(mpt.geometryN(i).asPoint());

		// 		std::vector<ign::geometry::LineString> subCl2cut = lsSplitterClSnappedOn.getSubLineStringsZ();

		// 		if (subCl2cut.size() == 1)
		// 			continue;

		// 		sCl2delete.insert(fEdgDoubl2Cut.getId());

		// 		for (size_t i = 0; i < subCl2cut.size(); ++i) {
		// 			ign::feature::Feature fNewEdgeDouble = fEdgDoubl2Cut;
		// 			fNewEdgeDouble.setGeometry(subCl2cut[i]);
		// 			_fsEdge->createFeature(fNewEdgeDouble);

		// 			//--
		// 			_shapeLogger->writeFeature("edgeClCutByCp", fEdgDoubl2Cut);
		// 			_logger->log(epg::log::INFO, "edge create by _cutClByCp: " + fNewEdgeDouble.getId());
		// 		}
		// 	}

		// 	for (std::set<std::string>::iterator sit = sCl2delete.begin(); sit != sCl2delete.end(); ++sit)
		// 		_fsEdge->deleteFeature(*sit);
		// }

		///
		///
		///
		bool CFeatGenerationOp::_getNearestCP(
			ign::feature::Feature const& fCP,
			std::map<std::string, ign::feature::Feature> & mCPNear
		) const {
			//--
			mCPNear[fCP.getId()] = fCP;

			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
			std::string const idName = context->getEpgParameters().getValue(ID).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const formOfWayName = themeParameters->getValue(FORM_OF_WAY_NAME).toString();
			double const distMergeCP = themeParameters->getValue(CP_MERGE_DIST_CP).toDouble();
			double const distMergeTractorCP = themeParameters->getValue(CP_MERGE_DIST_TRACTOR_CP).toDouble();

			//--
			std::string const& formOfWay = fCP.getAttribute(formOfWayName).toString();
			double mergeDist = _sFormOfWayException.find(formOfWay) != _sFormOfWayException.end() ? distMergeTractorCP : distMergeCP;

			//--
			ign::feature::FeatureFilter filterArroundCP("ST_DISTANCE(" + geomName + ", ST_SetSRID(ST_GeomFromText('" + fCP.getGeometry().toString() + "'),3035)) < "+ign::data::Double(mergeDist).toString());
			for (std::map<std::string, ign::feature::Feature>::iterator mit = mCPNear.begin(); mit != mCPNear.end(); ++mit) {
				epg::tools::FilterTools::addAndConditions(filterArroundCP, idName + " <> '" + mit->first + "'");	//(idName + " <> '" + fCP.getId() + "'");
			}

			ign::feature::FeatureIteratorPtr itArroundCP = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCP, filterArroundCP);

			if (!itArroundCP->hasNext())
				return false;

			while (itArroundCP->hasNext())
			{
				ign::feature::Feature const& fCPArround = itArroundCP->next();
				std::string const& formOfWayArround = fCPArround.getAttribute(formOfWayName).toString();

				if (formOfWayArround != formOfWay) {
					double mergeDistArround = _sFormOfWayException.find(formOfWayArround) != _sFormOfWayException.end() ? distMergeTractorCP : distMergeCP;
					if (mergeDistArround < mergeDist) {
						double dist = fCPArround.getGeometry().asPoint().distance(fCP.getGeometry().asPoint());

						if( dist > mergeDistArround )
							continue;
					}
				}

				_getNearestCP(fCPArround, mCPNear);
			}
			return true;
		}

		///
		///
		///
		void  CFeatGenerationOp::_snapCl2Cl(
			double distMaxClClosest
		) const {
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
			
			//--
			GraphType graphCl;
			_loadGraphCL(graphCl);

			//--
			ign::feature::FeatureFilter filterCL(_isClStatement);
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCL);

			boost::progress_display displayLoad(numFeatures, std::cout, "[ snapping CL to CL % complete ]\n");

			ign::feature::FeatureIteratorPtr itCL = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCL);
			std::map<std::string, std::pair<ENDING, ign::feature::Feature>> mFClModified;
			while (itCL->hasNext()) 
			{
				++displayLoad;

				ign::feature::Feature fCL = itCL->next();
				ign::geometry::LineString lsCl = fCL.getGeometry().asLineString();

				edge_descriptor edCl = graphCl.getInducedEdges(fCL.getId()).second[0].descriptor;

				if(graphCl.degree(graphCl.source(edCl)) == 1)
					_snapTo(distMaxClClosest, START, fCL, lsCl, mFClModified);
					
				if(graphCl.degree(graphCl.target(edCl)) == 1)
					_snapTo(distMaxClClosest, END, fCL, lsCl, mFClModified);
			}

			for (std::map<std::string, std::pair<ENDING, ign::feature::Feature>>::iterator mit = mFClModified.begin(); mit != mFClModified.end(); ++mit)
				_fsCL->modifyFeature(mit->second.second);
		}

		///
		///
		///
		void CFeatGenerationOp::_snapTo( 
			double distMaxClClosest,
			CFeatGenerationOp::ENDING ending,
			ign::feature::Feature & fCL,
			ign::geometry::LineString & newClGeom,
			std::map<std::string, std::pair<CFeatGenerationOp::ENDING, ign::feature::Feature>> & mFClModified
		) const {
			std::vector<std::pair<ENDING, ign::feature::Feature>> vClExtremityClose = _getClExtremityClose( distMaxClClosest, ending, fCL);

			ign::geometry::LineString const& clGeom = fCL.getGeometry().asLineString();
			ign::geometry::Point const& endingGeom = ending == START ? clGeom.startPoint() : clGeom.endPoint();
			ign::geometry::Point const& otherEndingGeom = ending == START ? clGeom.endPoint() : clGeom.startPoint();

			ign::geometry::MultiPoint mp;
			mp.addGeometry(ending == START ? clGeom.startPoint() : clGeom.endPoint() );
			
			for ( int i = vClExtremityClose.size()-1 ; i >= 0 ; --i )
			{	
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

			if (!vClExtremityClose.empty()) 
			{
				ign::geometry::Point newPt = mp.asMultiPoint().getCentroid();
				newPt.setZ(0);
				
				newClGeom.setPointN(newPt, ending == START ? 0 : newClGeom.numPoints()-1);
				fCL.setGeometry(newClGeom);

				ENDING ending_ = ending;
				std::map<std::string, std::pair<CFeatGenerationOp::ENDING, ign::feature::Feature>>::iterator mit = mFClModified.find(fCL.getId());
				if ( mit != mFClModified.end() && mit->second.first != ending )
					ending_ = BOTH;

				mFClModified[fCL.getId()] = std::make_pair(ending_, fCL);

				for ( size_t i = 0 ; i < vClExtremityClose.size() ; ++i ) 
				{
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
		std::vector<std::pair<CFeatGenerationOp::ENDING, ign::feature::Feature>> CFeatGenerationOp::_getClExtremityClose(
			double distMaxClClosest,
			CFeatGenerationOp::ENDING ending,
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
			epg::tools::FilterTools::addAndConditions(filterArroundCl, _isClStatement);
			filterArroundCl.setExtent(ptClCurr.getEnvelope().expandBy(distMaxClClosest));

			ign::feature::FeatureIteratorPtr itClArround = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterArroundCl);
			std::vector<std::pair<ENDING, ign::feature::Feature>> vConnectedCl;
			while (itClArround->hasNext())
			{
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
		void CFeatGenerationOp::_mergeIntersectingClWithGraph(
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
			ign::geometry::graph::tools::SnapRoundPlanarizer<GraphType> planarizerCl(graphCl, 1e4);

			//--
			ign::feature::FeatureFilter filterCL("((" + _mIsCountryStatement.begin()->second + ") OR (" + _mIsCountryStatement.rbegin()->second + "))");

			//--
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCL);
			boost::progress_display displayLoad(numFeatures, std::cout, "[ building CL planar graph % complete ]\n");

			ign::feature::FeatureIteratorPtr itCL = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCL);
			while (itCL->hasNext()) 
			{
				++displayLoad;

				ign::feature::Feature fCL = itCL->next();
				planarizerCl.addEdge(fCL.getGeometry().asLineString(), fCL.getId());
			}
			planarizerCl.planarize();

			//fusion des edges adjacents avec les mêmes origines
			_mergingEdgesByOrigin(graphCl);

			GraphType::edge_iterator eit, eitEnd;
			graphCl.edges(eit, eitEnd);
			boost::progress_display display(graphCl.numEdges(), std::cout, "[ merging CL % complete ]\n");
			while (eit != eitEnd) 
			{
				++display;

				std::vector<std::string> vClOrigins = graphCl.origins(*eit);
				ign::geometry::LineString lsCl = graphCl.getGeometry(*eit);

				if (graphCl.origins(*eit).size() == 1) {
					++eit;
					continue;
				}

				std::map<std::string, ign::feature::Feature> mIdClOriginsCountryFront, mIdClOriginsCountryBack;

				//recuperation des edges liés aux CLs
				for (std::vector<std::string>::iterator vit = vClOrigins.begin() ; vit != vClOrigins.end() ; ++vit)
				{
					ign::feature::Feature fCl;
					_fsCL->getFeatureById(*vit, fCl);

					if (fCl.getId().empty())
						continue;

					std::string countryCodeCl = fCl.getAttribute(countryCodeName).toString();
					// if (countryCodeCl == _vCountry[0]) 
					if (countryCodeCl.find(_vCountry.front()) != std::string::npos)
						mIdClOriginsCountryFront[fCl.getId()] = fCl;
					// else if (countryCodeCl == _vCountry[1])
					else if (countryCodeCl.find(_vCountry.back()) != std::string::npos)
						mIdClOriginsCountryBack[fCl.getId()] = fCl;
					else //ne devrait pas arriver
						continue;
				}

				if (mIdClOriginsCountryFront.size() == 0 || mIdClOriginsCountryBack.size() == 0) {
					++eit;
					continue;//pas de fusion, CL de un seul pays
				}

				//recuperation des portions d'edges associées et selection des CLs à fusionner
				std::set<std::string> sEdgesMerged;
				int nbEdgesMerged = -1;

				while (nbEdgesMerged != sEdgesMerged.size() )
				{
					nbEdgesMerged = sEdgesMerged.size();

					double distMin = std::numeric_limits<double>::infinity();
					std::pair<std::string, std::string> cl2merge;

					for (std::map<std::string, ign::feature::Feature>::iterator mit1 = mIdClOriginsCountryFront.begin(); mit1 != mIdClOriginsCountryFront.end(); ++mit1) 
					{
						if (sEdgesMerged.find(mit1->first) != sEdgesMerged.end())//deja utilise pour merge
							continue;

						//--
						std::string idEdgeLinked1 = mit1->second.getAttribute(linkedFeatIdName).toString();

						ign::feature::Feature fEdge1;
						_fsEdge->getFeatureById(idEdgeLinked1, fEdge1);

						ign::geometry::LineString lsEdge1 = fEdge1.getGeometry().asLineString();

						//--
						ign::geometry::GeometryPtr geomClEdge1(_getGeomProjClOnEdge(lsCl, lsEdge1, snapProjCl2edge));

						//--
						for (std::map<std::string, ign::feature::Feature>::iterator mit2 = mIdClOriginsCountryBack.begin(); mit2 != mIdClOriginsCountryBack.end(); ++mit2) 
						{
							if (sEdgesMerged.find(mit2->first) != sEdgesMerged.end())//deja utilise pour merge
								continue;

							//--
							std::string idEdgeLinked2 = mit2->second.getAttribute(linkedFeatIdName).toString();

							ign::feature::Feature fEdge2;
							_fsEdge->getFeatureById(idEdgeLinked2, fEdge2);

							ign::geometry::LineString lsEdge2 = fEdge2.getGeometry().asLineString();

							//--
							ign::geometry::GeometryPtr geomClEdge2(_getGeomProjClOnEdge(lsCl, lsEdge2, snapProjCl2edge));

							double hausdorffDistEdges = ign::geometry::algorithm::OptimizedHausdorffDistanceOp::distance(*geomClEdge1, *geomClEdge2);

							if (hausdorffDistEdges > distMaxEdges) {
								//DEBUG
								ign::feature::Feature feat;
								feat.setGeometry(lsEdge2);
								_shapeLogger->writeFeature("ClDebug", feat, mit1->first + " : " + idEdgeLinked1 + " " + mit2->first + " : " + idEdgeLinked2);

								continue;
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

					ign::feature::Feature fClNew = _attrMerger->merge(
						mIdClOriginsCountryFront[cl2merge.first], 
						mIdClOriginsCountryBack[cl2merge.second]
					);

					lsCl.setFillZ(0);
					fClNew.setGeometry(lsCl);

					std::string idCLNew = _idGeneratorCL->next();
					_fsCL->createFeature(fClNew, idCLNew);
				}
				++eit;
			}

			//suppression des CL sans #
			std::string query = "DELETE  FROM " + _fsCL->getTableName() + " WHERE NOT ("+_isClStatement+")";

			context->getDataBaseManager().getConnection()->update(query);

			//fusion des CL ayant les mêmes linkedFeatIdName et qui se touchent (ou presque)
			// A ce stade il y a des CL avec simple linkedFeatIdName et double linkedFeatIdName ...!?
			ign::feature::FeatureFilter mergeFilter(linkedFeatIdName + " LIKE '%#%'");
			ome2::calcul::detail::ClMerger::mergeAll(_fsCL, mergeFilter, _idGeneratorCL.get());
		}

		///
		///
		///
		void CFeatGenerationOp::_deleteClByAngleAndDistEdges() const 
		{	
			_logger->log(epg::log::TITLE, "[ BEGIN DELETE CL BY ANGLE EDGES FOR : " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
			std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			double angleMax = themeParameters->getValue(CL_EDGE_MAX_ANGLE).toDouble();
			double const distMax = themeParameters->getValue(CL_EDGE_MAX_DIST).toDouble();
			double const snapProjCl2edge = themeParameters->getValue(CL_SNAP_PROJ_CL_2_EDGE_DIST).toDouble();
			angleMax = angleMax * M_PI / 180;

			//--
			ign::feature::FeatureFilter filterCL(_isClStatement);

			int numCl2delete = -1;
			while (numCl2delete != 0)
			{
				//--
				GraphType graphCl;
				_loadGraphCL(graphCl);

				//--
				size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCL);
				boost::progress_display display(numFeatures, std::cout, "[ deleting CL by angle between linked edges % complete ]\n");

				//--
				std::set<std::string> sCl2delete;
				ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCL);
				while (it->hasNext())
				{
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
			_logger->log(epg::log::TITLE, "[ END DELETE CL BY ANGLE EDGES FOR :" + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
		}

		///
		///
		///
		void CFeatGenerationOp::_updateGeomCL(double snapProjCl2edge) const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN UPDATE GEOM CL " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();
			std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			std::string const fictitiousFieldName = themeParameters->getValue(EDGE_FICTITIOUS_NAME).toString();

			std::string countryCode1 = _vCountry.front();
			std::string countryCode2 = _vCountry.back();

			//--
			ign::feature::FeatureFilter filterCL(_isClStatement);
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCL);
			boost::progress_display display(numFeatures, std::cout, "[ computing CL geometry by linked edges interpolation % complete ]\n");

			//--
			ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCL);
			std::set<std::string> sCL2delete;
			while (it->hasNext())
			{
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

			_logger->log(epg::log::TITLE, "[ BEGIN UPDATE GEOM CL " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
		}

		///
		///
		///
		ign::geometry::LineString CFeatGenerationOp::_computeMeanGeom(
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
		ign::geometry::Geometry* CFeatGenerationOp::_getGeomProjClOnEdge(
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

		///
		///
		///
		void CFeatGenerationOp::_deleteCLUnderThreshold() const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN CLEAN CL UNDER THRESHOLD " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());

			//--
			params::ThemeParameters* themeParameters = params::ThemeParametersS::getInstance();
			double minLength = themeParameters->getValue(CL_MIN_LENGTH).toDouble();

			//--
			GraphType graph;
			_loadGraphCL(graph);

			// patience
			boost::progress_display display(graph.numEdges(), std::cout, "[ deleting short CL % complete ]\n");

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

				for (edges_path_const_iterator pit = path.begin() ; pit != path.end() ; ++pit)
				{
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
			_logger->log(epg::log::TITLE, "[ END CLEAN CL UNDER THRESHOLD " + _borderCode + " ] : " + epg::tools::TimeTools::getTime());
		}

		///
		///
		///
		double CFeatGenerationOp::_getLength(
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
		CFeatGenerationOp::edges_path CFeatGenerationOp::_getPath(
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
		void CFeatGenerationOp::_loadGraphCL(GraphType & graphCL) const
		{
			//--
			epg::Context* context = epg::ContextS::getInstance();
			std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

			//--
			ign::feature::FeatureFilter filterCL(_isClStatement);
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsCL, filterCL);
			boost::progress_display display(numFeatures, std::cout, "[ loading CL graph % complete ]\n");

			//--
			ign::geometry::graph::builder::SimpleGraphBuilder<GraphType> graphBuilder(graphCL, 0.01);
			ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsCL, filterCL);
			while (it->hasNext()) {
				++display;
				ign::feature::Feature fCL = it->next();
				graphBuilder.addEdge(fCL.getGeometry().asLineString(), fCL.getId());
			}
		}

		///
		///
		///
		void CFeatGenerationOp::_loadGraphEdges(
			std::string const& country,
			GraphType & graphEdges
		) const {
			ign::feature::FeatureFilter filter(_mIsCountryStatement.at(country));

			//--
			size_t numFeatures = ome2::feature::sql::NotDestroyedTools::NumFeatures(*_fsEdge, filter);
			boost::progress_display display(numFeatures, std::cout, "[ loading " + country + " edges graph % complete ]\n");
			
			ign::feature::FeatureIteratorPtr it = ome2::feature::sql::NotDestroyedTools::GetFeatures(*_fsEdge, filter);

			ign::geometry::graph::builder::SimpleGraphBuilder<GraphType> graphBuilder(graphEdges, 0.01);
			while (it->hasNext()) 
			{
				++display;

				ign::feature::Feature fedge = it->next();
				graphBuilder.addEdge(fedge.getGeometry().asLineString(), fedge.getId());
			}
		}

		///
		///
		///
		bool CFeatGenerationOp::_isConnectedEdges(
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
		std::pair<bool,std::pair<std::string, std::string>> CFeatGenerationOp::_getClLinkedEdges(
			std::string const& linkedFeatIdName,
			GraphType const& graphCL,
			GraphType::edge_descriptor eCl
		) const {
			std::string idCl = graphCL.origins(eCl)[0];

			ign::feature::Feature clFeat;
			_fsCL->getFeatureById(idCl, clFeat);

			if (clFeat.getId().empty())
				return std::make_pair(false, std::make_pair("", ""));

			std::string linkedFeat = clFeat.getAttribute(linkedFeatIdName).toString();

			std::vector<std::string> vEdgeslinkedJ;
			epg::tools::StringTools::Split(linkedFeat, "#", vEdgeslinkedJ);

			return std::make_pair(true,std::make_pair(vEdgeslinkedJ[0], vEdgeslinkedJ[1]));
		}

		///
		///
		///
		bool CFeatGenerationOp::_areParallelEdges(
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
		ign::geometry::Point CFeatGenerationOp::_getLinkedEdgesConnectingPoint(
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
		void CFeatGenerationOp::_setContinuityCl(GraphType const& graphCL) const
		{
			_logger->log(epg::log::TITLE, "[ BEGIN SET CL CONTINUITY ] : " + epg::tools::TimeTools::getTime());

			epg::Context* context = epg::ContextS::getInstance();
			std::string const linkedFeatIdName = context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString();

			GraphType graphEdges1, graphEdges2;
			_loadGraphEdges(_vCountry.front(), graphEdges1);
			_loadGraphEdges(_vCountry.back(), graphEdges2);

			boost::progress_display display(graphCL.numVertices(), std::cout, "[ set CL continuity % complete ]\n");
			std::map<std::string, ign::feature::Feature> mClModified;

			GraphType::vertex_iterator vit, vitEnd;
			graphCL.vertices(vit, vitEnd);
			while (vit != vitEnd)
			{
				++display;

				if (graphCL.degree(*vit) < 2 ) {
					++vit;
					continue;
				}

				std::vector< GraphType::oriented_edge_descriptor > vClsIncidentTemp;
				graphCL.incidentEdges(*vit, vClsIncidentTemp);

				//--
				std::vector<std::vector< GraphType::oriented_edge_descriptor >> vVClsTrueIncident;
				std::set<size_t> sTreated;
				for (size_t i = 0; i < vClsIncidentTemp.size()-1; ++i)
				{
					if ( sTreated.find(i) != sTreated.end() ) continue;

					std::vector< GraphType::oriented_edge_descriptor > vClsIncidentTempConnectI;
					vClsIncidentTempConnectI.push_back(vClsIncidentTemp[i]);

					std::pair<bool, std::pair<std::string, std::string>> pLinkedEdgesI = _getClLinkedEdges(linkedFeatIdName, graphCL, vClsIncidentTemp[i].descriptor);

					//si CL n'existe plus
					if (!pLinkedEdgesI.first)
						continue;

					for (size_t j = i + 1; j < vClsIncidentTemp.size(); ++j)
					{
						if ( sTreated.find(j) != sTreated.end() )
							continue;

						std::pair<bool, std::pair<std::string, std::string>> pLinkedEdgesJ = _getClLinkedEdges(linkedFeatIdName, graphCL, vClsIncidentTemp[j].descriptor);

						if (!pLinkedEdgesJ.first)
							continue;

						bool isConnected1 = pLinkedEdgesI.first == pLinkedEdgesJ.first || _isConnectedEdges(graphEdges1, pLinkedEdgesI.second.first, pLinkedEdgesJ.second.first);
						bool isConnected2 = pLinkedEdgesI.second == pLinkedEdgesJ.second || _isConnectedEdges(graphEdges2, pLinkedEdgesI.second.second, pLinkedEdgesJ.second.second);

						// tester si vClsIncidentTemp[i].descriptor et vClsIncidentTemp[j].descriptor sont paralleles
						if (_areParallelEdges(graphCL, vClsIncidentTemp[i].descriptor, vClsIncidentTemp[j].descriptor))
						{
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
				for (size_t j = 0; j < vVClsTrueIncident.size(); ++j)
				{
					std::vector< GraphType::oriented_edge_descriptor > vClsTrueIncident = vVClsTrueIncident[j];

					ign::geometry::MultiPoint multiPtToConnect;

					for (size_t i = 0; i < vClsTrueIncident.size(); ++i)
					{
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
					for (size_t i = 0; i < vClsTrueIncident.size(); ++i)
					{
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
		void CFeatGenerationOp::_mergingEdgesByOrigin(
			GraphType & graph
		) const {
			std::vector<GraphType::edges_path> vGatheredEdges;
			std::set< GraphType::edge_descriptor > sVisitedEdges;

			boost::progress_display display(graph.numEdges(), std::cout, "[ gathering edges by origins % complete ]\n");

			GraphType::edge_iterator eit, eend;
			for (graph.edges(eit, eend); eit != eend; ++eit)
			{
				++display;

				if (graph.target(*eit) == graph.source(*eit)) continue;
				if (sVisitedEdges.find(*eit) != sVisitedEdges.end()) continue;

				sVisitedEdges.insert(*eit);

				std::vector<std::string> vOriginsRef = graph.origins(*eit);

				GraphType::oriented_edge_descriptor tPivot[] = {
					GraphType::oriented_edge_descriptor(*eit, ign::graph::DIRECT),
					GraphType::oriented_edge_descriptor(*eit, ign::graph::REVERSE)
				};

				//arcs a fusionner
				GraphType::edges_path path;
				path.push_back(tPivot[0]);

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

						if (vOriginsRef != graph.origins(nextEdge.descriptor)) break;

						sVisitedEdges.insert(nextEdge.descriptor);
						if (i == 0) {
							path.push_back(nextEdge);
						} else {
							path.push_front(nextEdge.getReverse());
						}
						
						vTarget = graph.target(nextEdge);

						if (graph.degree(vTarget) != 2) break;
					}
					if (isLoop) break;
				}

				if (path.size() > 1 && !isLoop)
					vGatheredEdges.push_back(path);
			}

			//--
			boost::progress_display display2(vGatheredEdges.size(), std::cout, "[ merging CL by origins % complete ]\n");

			typename std::vector<GraphType::edges_path>::iterator vit = vGatheredEdges.begin();
			for (; vit != vGatheredEdges.end(); ++vit, ++display2)
			{	
				std::vector<ign::geometry::Point> vInter;
				for (GraphType::edges_path_iterator it = vit->begin(); it != vit->end(); ++it)
				{
					ign::geometry::LineString lsTemp = graph.getGeometry(it->descriptor);
					if (it->direction == ign::graph::DIRECT)
					{
						ign::geometry::LineString::iterator itLs = lsTemp.begin();
						for (++itLs ; itLs != lsTemp.end() ; ++itLs)
							vInter.push_back(*itLs);
					} else {
						ign::geometry::LineString::reverse_iterator itLs = lsTemp.rbegin();
						for (++itLs ; itLs != lsTemp.rend() ; ++itLs)
							vInter.push_back(*itLs);
					}
				}
				vInter.pop_back();

				GraphType::vertex_descriptor vSource = graph.source(*vit->begin());
				GraphType::vertex_descriptor vTarget = graph.target(*vit->rbegin());

				GraphType::oriented_edge_descriptor newD = graph.addEdge(vSource, vTarget, vInter, graph[vit->begin()->descriptor]);

				for (GraphType::edges_path_iterator it = vit->begin(); it != vit->end(); )
				{
					typename GraphType::edges_path_iterator itTemp = it++;
					graph.removeEdge(itTemp->descriptor);
				}
			}
		}
	}
}