#ifndef _APP_CALCUL_DETAIL_CLMERGER_H_
#define _APP_CALCUL_DETAIL_CLMERGER_H_

//EPG
#include <epg/log/EpgLogger.h>
#include <epg/sql/tools/IdGeneratorInterface.h>

//SOCLE
#include <ign/feature/FeatureStore.h>


namespace app{
namespace calcul{
namespace detail{

    class ClMerger {
        public:

        //--
        static void mergeAll(
            epg::sql::tools::IdGeneratorInterface* idGenerator,
            ign::feature::FeatureStore* fsCl,
            ign::feature::FeatureFilter const& filter
        );

        //--
        static ign::geometry::LineString merge(
            ign::feature::FeatureStore* fsCl,
            ign::feature::Feature const& refClFeat,
            std::string const& linkedFeatureId,
            std::set<std::string> & sMergedCl
        );

        private:
            //--
            ign::feature::FeatureStore*              _fsCl;
            //--
		    epg::log::EpgLogger*                     _logger;
            //--
            epg::sql::tools::IdGeneratorInterface*   _idGenerator;

        private:

        //--
        ClMerger(
            epg::sql::tools::IdGeneratorInterface* idGenerator,
            ign::feature::FeatureStore* fsCl
        );

        //--
        ~ClMerger();

        //--
        void _mergeAll(ign::feature::FeatureFilter const& filter) const;

        //--
        ign::geometry::LineString _mergecl(
            ign::feature::Feature const& refClFeat,
            std::string const& linkedFeatureId,
            std::set<std::string> & sMergedCl
        ) const;

    };

}
}
}

#endif