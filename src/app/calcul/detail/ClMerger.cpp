//APP
#include <app/calcul/detail/ClMerger.h>

//EPG
#include <epg/Context.h>

//SOCLE
#include <ign/feature/FeatureStore.h>
#include <ign/geometry/algorithm/LineMergerOpGeos.h>

namespace app{
namespace calcul{
namespace detail{

    ///
    ///
    ///
    void ClMerger::mergeAll(
        epg::sql::tools::IdGeneratorInterface* idGenerator,
        ign::feature::FeatureStore* fsCl,
        ign::feature::FeatureFilter const& filter
    ) {
        ClMerger merger(idGenerator, fsCl);
        merger._mergeAll(filter);
    }

    ///
    ///
    ///
    ign::geometry::LineString ClMerger::merge(
        ign::feature::FeatureStore* fsCl,
        ign::feature::Feature const& refClFeat,
        std::string const& linkedFeatureId,
        std::set<std::string> & sMergedCl
    ) {
        ClMerger merger(0, fsCl);
        return merger._mergecl(refClFeat, linkedFeatureId, sMergedCl);
    }

    ///
    ///
    ///
    ClMerger::ClMerger(
        epg::sql::tools::IdGeneratorInterface* idGenerator,
        ign::feature::FeatureStore* fsCl
    ):
        _fsCl(fsCl),
        _logger(0),
        _idGenerator(idGenerator)
    {
        _logger = epg::log::EpgLoggerS::getInstance();
    };

    ///
    ///
    ///
    ClMerger::~ClMerger() {
    };

    ///
    ///
    ///
    void ClMerger::_mergeAll(ign::feature::FeatureFilter const& filter) const {
        epg::Context* context = epg::ContextS::getInstance();
        epg::params::EpgParameters const& epgParams = context->getEpgParameters();
        std::string const linkedFeatureIdName = epgParams.getValue(LINKED_FEATURE_ID).toString();

        std::set<std::string> sMergedCl;
        ign::feature::FeatureIteratorPtr itCl = _fsCl->getFeatures(filter);
        while (itCl->hasNext())
        {
            ign::feature::Feature fCl = itCl->next();
            ign::geometry::LineString const& clGeom = fCl.getGeometry().asLineString();
            std::string linkedFeatureId = fCl.getAttribute(linkedFeatureIdName).toString();

            if ( sMergedCl.find(fCl.getId()) != sMergedCl.end() ) continue;

            std::set<std::string> sMergedCl2;
            ign::geometry::LineString mergedClGeom = _mergecl(fCl, linkedFeatureId, sMergedCl2);

            if (sMergedCl2.size() < 2) continue;

            sMergedCl.insert(sMergedCl2.begin(), sMergedCl2.end());

            fCl.setGeometry(mergedClGeom);
            _fsCl->createFeature(fCl, _idGenerator->next());

            _logger->log(epg::log::INFO, "New merged CL [id] " + fCl.getId());
        }

        for(std::set<std::string>::const_iterator sit = sMergedCl.begin() ; sit != sMergedCl.end() ; ++sit) {
            _fsCl->deleteFeature(*sit);
        }
    }

    ///
    ///
    ///
    ign::geometry::LineString ClMerger::_mergecl(
        ign::feature::Feature const& refClFeat,
        std::string const& linkedFeatureId,
        std::set<std::string> & sMergedCl
    ) const {
        epg::Context* context = epg::ContextS::getInstance();
        epg::params::EpgParameters const& epgParams = context->getEpgParameters();
        std::string const linkedFeatureIdName = epgParams.getValue(LINKED_FEATURE_ID).toString();

        ign::geometry::LineString const& refClGeom = refClFeat.getGeometry().asLineString();
        ign::geometry::Point startPoint = refClGeom.startPoint();
        ign::geometry::Point endPoint = refClGeom.endPoint();
        bool bNewStart = true;
        bool bNewEnd = true;

        std::vector<ign::feature::Feature> vCandidates;
        ign::feature::FeatureIteratorPtr itCl = _fsCl->getFeatures(linkedFeatureIdName + " LIKE '%" + linkedFeatureId + "%'");
        while (itCl->hasNext())
        {
            ign::feature::Feature const& fCl = itCl->next();
            vCandidates.push_back(fCl);
        }

        std::vector<ign::geometry::LineString> vGeom2Merge;
        vGeom2Merge.push_back(refClGeom);
        sMergedCl.insert(refClFeat.getId());

        do {
            std::vector<ign::feature::Feature>::const_iterator vit;
            size_t before = sMergedCl.size();
            for ( vit = vCandidates.begin() ; vit != vCandidates.end() ; ++vit ) {
                if (sMergedCl.find(vit->getId()) != sMergedCl.end()) continue;
                
                ign::geometry::LineString const& clGeom = vit->getGeometry().asLineString();

                //DEBUG
                double bTouchingEnd_d = startPoint.distance(clGeom.endPoint());
                double bTouchingStart_d = startPoint.distance(clGeom.startPoint());

                bool bTouchingEnd = startPoint.distance(clGeom.endPoint()) < 1e-1;
                bool bTouchingStart = startPoint.distance(clGeom.startPoint()) < 1e-1;
                if ( bTouchingEnd || bTouchingStart ) {
                    sMergedCl.insert(vit->getId());
                    vGeom2Merge.push_back(clGeom);
                    if (bTouchingEnd) vGeom2Merge.back().endPoint() = startPoint;
                    if (bTouchingStart) vGeom2Merge.back().startPoint() = startPoint;
                    startPoint = bTouchingEnd ? clGeom.startPoint() : clGeom.endPoint();
                    break;
                }
            }
            if (before == sMergedCl.size()) bNewStart = false;
        }while (bNewStart);

        do {
            std::vector<ign::feature::Feature>::const_iterator vit;
            size_t before = sMergedCl.size();
            for ( vit = vCandidates.begin() ; vit != vCandidates.end() ; ++vit ) {
                if (sMergedCl.find(vit->getId()) != sMergedCl.end()) continue;

                ign::geometry::LineString const& clGeom = vit->getGeometry().asLineString();
                bool bTouchingEnd = endPoint.distance(clGeom.endPoint()) < 1e-5;
                bool bTouchingStart = endPoint.distance(clGeom.startPoint()) < 1e-5;
                if ( bTouchingEnd || bTouchingStart ) {
                    sMergedCl.insert(vit->getId());
                    vGeom2Merge.push_back(clGeom);
                    if (bTouchingEnd) vGeom2Merge.back().endPoint() = endPoint;
                    if (bTouchingStart) vGeom2Merge.back().startPoint() = endPoint;
                    endPoint = bTouchingEnd ? clGeom.startPoint() : clGeom.endPoint();
                    break;
                }
            }
            if (before == sMergedCl.size()) bNewEnd = false;
        } while (bNewEnd);

        if (vGeom2Merge.size() == 1) return vGeom2Merge.front();

        std::vector<ign::geometry::LineString> vMergedGeom = ign::geometry::algorithm::LineMergerOpGeos::MergeLineStrings(vGeom2Merge);
        if ( vMergedGeom.size() > 1 ) _logger->log(epg::log::WARN, "Merging adjacent CL gives a MultilineString [ref CL id] " + refClFeat.getId());

        return vMergedGeom.front();
    };

}
}
}