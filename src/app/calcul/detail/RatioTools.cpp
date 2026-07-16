#include <app/calcul/detail/RatioTools.h>

// EPG
#include <epg/tools/StringTools.h>

namespace app
{
    namespace calcul
    {
        namespace detail
        {

            ///
            ///
            ///
            double RatioTools::GetRatio(
                ign::geometry::LineString const& ls,
                std::string const& country,
                ign::feature::FeatureStore * fs1,
                ign::feature::FeatureStore * fs2
            ) {
                epg::Context *context = epg::ContextS::getInstance();
                epg::params::EpgParameters const& epgParams = context->getEpgParameters();
                std::string const geomName = epgParams.getValue(GEOM).toString();
                std::string const countryCodeName = epgParams.getValue(COUNTRY_CODE).toString();

                std::set<std::string> sCountry;
                epg::tools::StringTools::Split(country, "#", sCountry);

                ign::geometry::GeometryPtr unionPtr(new ign::geometry::Polygon());

                ign::feature::FeatureFilter filter;
                filter.setExtent(ls.getEnvelope());
                ign::feature::FeatureIteratorPtr it1 = fs1->getFeatures(filter);;
                while (it1->hasNext())
                {
                    ign::feature::Feature feat = it1->next();
                    std::string const& featCountry = feat.getAttribute(countryCodeName).toString();
                    ign::geometry::Geometry const& geom = feat.getGeometry();

                    bool hasCommonCountry = false;
                    for ( std::set<std::string>::const_iterator sit = sCountry.begin() ; sit != sCountry.end() ; ++sit )
                        if( featCountry.find(*sit) != std::string::npos )
                            hasCommonCountry = true;

                    if( !hasCommonCountry )
                        continue;

                    if( !geom.intersects(ls) )
                        continue;

                    unionPtr.reset(unionPtr->Union(geom));
                }

                if( fs2 != 0 )
                {
                    ign::feature::FeatureIteratorPtr it2 = fs2->getFeatures(filter);;
                    while (it2->hasNext())
                    {
                        ign::feature::Feature feat = it2->next();
                        std::string const& featCountry = feat.getAttribute(countryCodeName).toString();
                        ign::geometry::Geometry const& geom = feat.getGeometry();

                        bool hasCommonCountry = false;
                        for ( std::set<std::string>::const_iterator sit = sCountry.begin() ; sit != sCountry.end() ; ++sit )
                            if( featCountry.find(*sit) != std::string::npos )
                                hasCommonCountry = true;

                        if( !hasCommonCountry )
                            continue;

                        if( !geom.intersects(ls) )
                            continue;

                        unionPtr.reset(unionPtr->Union(geom));
                    }
                }
                
                if(unionPtr->isEmpty() || unionPtr->isNull()) return 0;

                ign::geometry::GeometryPtr resultPtr(unionPtr->Intersection(ls));

                double lengthInter = _GetLength(*resultPtr);

                return lengthInter / ls.length();
            }

            ///
            ///
            ///
            double RatioTools::_GetLength(
                ign::geometry::Geometry const& geom
            ) {
                double length = 0;

                ign::geometry::Geometry::GeometryType geomType = geom.getGeometryType();
                switch( geomType )
                {
                    case ign::geometry::Geometry::GeometryTypeNull :
                    case ign::geometry::Geometry::GeometryTypePoint :
                    case ign::geometry::Geometry::GeometryTypeMultiPoint :
                    case ign::geometry::Geometry::GeometryTypeTriangle :
                    case ign::geometry::Geometry::GeometryTypeTriangulatedSurface :
                    case ign::geometry::Geometry::GeometryTypePolyhedralSurface :
                    case ign::geometry::Geometry::GeometryTypePolygon :
                    case ign::geometry::Geometry::GeometryTypeMultiPolygon :
                        return 0;
                    case ign::geometry::Geometry::GeometryTypeLineString :
                        {
                            return geom.asLineString().length();
                        }
                        
                    case ign::geometry::Geometry::GeometryTypeMultiLineString : 
                        {
                            ign::geometry::MultiLineString const& mls = geom.asMultiLineString();
                            for( size_t i = 0 ; i < mls.numGeometries() ; ++i ) {
                                length += mls.lineStringN(i).length();
                            }
                            return length;
                        }
                    
                    case ign::geometry::Geometry::GeometryTypeGeometryCollection :
                        {
                            ign::geometry::GeometryCollection const& collection = geom.asGeometryCollection();
                            for( size_t i = 0 ; i < collection.numGeometries() ; ++i ) {
                                length += _GetLength(collection.geometryN(i));
                            }
                            return length;
                        }
                    default :
                        IGN_THROW_EXCEPTION( "Geometry type not allowed" );
                }
            }
        }
    }
}
