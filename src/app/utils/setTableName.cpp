#include <app/utils/setTableName.h>

//EPG
#include <epg/Context.h>
#include <epg/tools/StringTools.h>

//APP
#include <app/params/ThemeParameters.h>

namespace app{
namespace utils{

    //--
    void setTableName(TN_PARAMETERS TABLE_NAME_PARAM) {
        epg::Context* context = epg::ContextS::getInstance();
        app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();

        std::string const tableName = themeParameters->getValue(TABLE_NAME_PARAM).toString();

        std::vector<std::string> vNameParts;
        epg::tools::StringTools::Split(tableName, ".", vNameParts);

        if( vNameParts.size() != 2 ) return;

        context->getDataBaseManager().addSchemaToSearchPath(vNameParts.front());
        
        themeParameters->setParameter(TABLE_NAME_PARAM, ign::data::String(vNameParts.back()));
    }

    //--
    void setTableName(EPG_PARAMETERS TABLE_NAME_PARAM) {
        epg::Context* context = epg::ContextS::getInstance();
        epg::params::EpgParameters & epgParams = context->getEpgParameters();

        std::string const tableName = epgParams.getValue(TABLE_NAME_PARAM).toString();

        std::vector<std::string> vNameParts;
        epg::tools::StringTools::Split(tableName, ".", vNameParts);

        if( vNameParts.size() != 2 ) return;

        context->getDataBaseManager().addSchemaToSearchPath(vNameParts.front());
        
        epgParams.setParameter(TABLE_NAME_PARAM, ign::data::String(vNameParts.back()));
    }
}
}