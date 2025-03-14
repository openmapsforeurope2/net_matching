//APP
#include <app/utils/createCpClTables.h>
#include <app/params/ThemeParameters.h>

//EPG
#include <epg/Context.h>


namespace app{
namespace utils{

    //--
    void createCpClTables(std::string const& edgeTableName) {
        epg::Context* context = epg::ContextS::getInstance();
        app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();

        std::string const idName = context->getEpgParameters().getValue(ID).toString();
	    std::string const geomName = context->getEpgParameters().getValue(GEOM).toString();
        std::string const countryCodeName = context->getEpgParameters().getValue(COUNTRY_CODE).toString();

        std::string const cpTableName = themeParameters->getValue(CP_TABLE).toString();
        if (!context->getDataBaseManager().tableExists(cpTableName)) {
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
                << "ALTER TABLE " << cpTableName << " ADD COLUMN " << context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString() << " character varying(255);";

            context->getDataBaseManager().getConnection()->update(ss.str());
        }
        std::string const clTableName = themeParameters->getValue(CL_TABLE).toString();
        if (!context->getDataBaseManager().tableExists(clTableName)) {
            std::ostringstream ss;
            ss << "DROP TABLE IF EXISTS " << clTableName << " ;";
            ss << "CREATE TABLE " << clTableName
                << " AS TABLE " << edgeTableName
                << " WITH NO DATA;"
                << "ALTER TABLE " << clTableName << " ALTER COLUMN "
                << geomName << " type geometry(LineStringZ, 0);"
                << "ALTER TABLE " << clTableName << " ALTER COLUMN "
                << idName << " type varchar(255);"
                << "ALTER TABLE " << clTableName << " ALTER COLUMN "
                << countryCodeName << " type varchar(255);"
                << "ALTER TABLE " << clTableName << " ALTER COLUMN "
                << "w_national_identifier" << " type varchar(255);"
                << "ALTER TABLE " << clTableName << " ADD COLUMN " << context->getEpgParameters().getValue(LINKED_FEATURE_ID).toString() << " character varying(255);";
            context->getDataBaseManager().getConnection()->update(ss.str());
        }
    }
}
}