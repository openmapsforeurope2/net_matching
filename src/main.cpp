
// BOOST
#include <boost/program_options.hpp>

// EPG
#include <epg/Context.h>
#include <epg/log/EpgLogger.h>
#include <epg/log/ShapeLogger.h>
#include <epg/tools/TimeTools.h>
#include <epg/params/tools/loadParameters.h>

// OME2
#include <ome2/utils/setTableName.h>

// APP
#include <app/params/ThemeParameters.h>
#include <app/step/tools/initSteps.h>
#include <app/utils/createCpClTables.h>


namespace po = boost::program_options;

int main(int argc, char *argv[])
{
   // ign::geometry::PrecisionModel::SetDefaultPrecisionModel(ign::geometry::PrecisionModel(ign::geometry::PrecisionModel::FIXED, 1.0e5, 1.0e7) );

    epg::Context* context = epg::ContextS::getInstance();

    std::string     logDirectory = "";
    std::string     epgParametersFile = "";
    std::string     themeParametersFile = "";
    std::string     dbName = "";
    std::string     suffix = "";
    std::string     areaSuffix = "";
    std::string     stepCode = "";
    std::string     borderCode = "";
    std::string     table = "";
    bool            verbose = true;

    epg::step::StepSuite< app::params::ThemeParametersS > stepSuiteWaterCourseLink, stepSuiteRoadLink, stepSuiteRailwayLink;
	app::step::tools::initStepsWatercourseLink(stepSuiteWaterCourseLink);
	app::step::tools::initStepsRoadLink(stepSuiteRoadLink);
	app::step::tools::initStepsRailwayLink(stepSuiteRailwayLink);

	std::ostringstream OperatorDetail;
	OperatorDetail << "set step :" << std::endl
        << "road_link:" << std::endl
		<< stepSuiteRoadLink.toString() 
		<< "watercourse_link:" << std::endl
		<< stepSuiteWaterCourseLink.toString()
        << "railway_link:" << std::endl
		<< stepSuiteRailwayLink.toString();

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("c" , po::value< std::string >(&epgParametersFile)     , "conf file" )
        ("d" , po::value< std::string >(&dbName)                , "data base name" )
        ("t" , po::value< std::string >(&table)                 , "table" )
        ("s" , po::value< std::string >(&suffix)                , "working table suffix" )
        ("as", po::value< std::string >(&areaSuffix)            , "area working tables suffix" )
        ("sp", po::value< std::string >(&stepCode), OperatorDetail.str().c_str())
    ;

    //main log
    std::string     logFileName = "log.txt";
    std::ofstream   logFile( logFileName.c_str() ) ;

    logFile << "[START] " << epg::tools::TimeTools::getTime() << std::endl;

    int returnValue = 0;
    try{

        po::parsed_options parsed = po::command_line_parser(argc, argv)
                                    .options(desc)
                                    .allow_unregistered()
                                    .run();

        po::variables_map vm;
        po::store( parsed, vm );
        po::notify( vm );

        if ( vm.count( "help" ) ) {
            std::cout << desc << std::endl;
            return 1;
        }

        // Récupérer les arguments libres (non reconnus)
        std::vector<std::string> countries = po::collect_unrecognized(parsed.options, po::include_positional);

        if ( countries.size() != 2 ) {
            std::string mError = "spécifier au moins deux et seulement deux pays en argument";
            IGN_THROW_EXCEPTION(mError);
        }
        if( countries.front() > countries.back() )
            std::swap(countries.front(), countries.back());
        borderCode = countries.front()+"#"+countries.back();

        epg::step::StepSuite< app::params::ThemeParametersS >* stepSuitePtr = 0;
        if (table == "watercourse_link")
			stepSuitePtr = &stepSuiteWaterCourseLink;
		else if (table == "road_link")
			stepSuitePtr = &stepSuiteRoadLink;
        else if (table == "railway_link")
			stepSuitePtr = &stepSuiteRailwayLink;
        else {
            std::string mError = "unknown table " + table;
            IGN_THROW_EXCEPTION(mError);
        }

        if (stepCode.empty()) stepCode = stepSuitePtr->getStepsRange();

        //parametres EPG
		context->loadEpgParameters( epgParametersFile, table );

        //Initialisation du log de prod
        logDirectory = context->getConfigParameters().getValue( LOG_DIRECTORY ).toString();

        //test si le dossier de log existe sinon le creer
        boost::filesystem::path logDir(logDirectory);
        if (!boost::filesystem::is_directory(logDir))
        {
            if (!boost::filesystem::create_directory(logDir))
            {
                std::string mError = "impossible to create directory " + logDirectory;
                IGN_THROW_EXCEPTION(mError);
            }
        }

        //repertoire de travail
        context->setLogDirectory( logDirectory );

        //epg logger
        epg::log::EpgLogger* logger = epg::log::EpgLoggerS::getInstance();
        // logger->setProdOfstream( logDirectory+"/net_matching.log" );
        logger->setDevOfstream( context->getLogDirectory()+"/net_matching.log" );
        
        //theme parameters
        themeParametersFile = context->getConfigParameters().getValue( THEME_PARAMETER_FILE ).toString();
		app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
        epg::params::tools::loadParams( *themeParameters, themeParametersFile, borderCode );
        if (themeParameters->getParameter(COUNTRY_CODE_W).getValue().toString() == "")
            IGN_THROW_EXCEPTION("country code " + borderCode + " unknown in theme parameter file");

        //info de connection db
        context->loadEpgParameters( themeParameters->getValue(DB_CONF_FILE).toString() );
        if( dbName != "" )
            context->getConfigParameters().setParameter(DATABASE, ign::data::String(dbName));

        //tables des surfaces
        if ( !areaSuffix.empty() && table == "watercourse_link" ) {
            std::string watercourseTableBaseName = themeParameters->getValue(WATERCOURSE_AREA_TABLE_BASE).toString();
            std::string watercourseTableName = watercourseTableBaseName + "_" + countries.front() + "_" + countries.back() + "_" + areaSuffix;
            themeParameters->setParameter(WATERCOURSE_AREA_TABLE, ign::data::String(watercourseTableName));

            std::string standingWaterTableBaseName = themeParameters->getValue(STANDING_WATER_TABLE_BASE).toString();
            std::string standingWaterTableName = standingWaterTableBaseName + "_" + countries.front() + "_" + countries.back() + "_" + areaSuffix;
            themeParameters->setParameter(STANDING_WATER_TABLE, ign::data::String(standingWaterTableName));

            std::string matchedWatercourseTableBaseName = themeParameters->getValue(MATCHED_WATERCOURSE_AREA_TABLE_BASE).toString();
            std::string matchedWatercourseTableName = matchedWatercourseTableBaseName + "_" + countries.front() + "_" + countries.back() + "_" + areaSuffix;
            themeParameters->setParameter(MATCHED_WATERCOURSE_AREA_TABLE, ign::data::String(matchedWatercourseTableName));

            std::string matchedStandingWaterTableBaseName = themeParameters->getValue(MATCHED_STANDING_WATER_TABLE_BASE).toString();
            std::string matchedStandingWaterTableName = matchedStandingWaterTableBaseName + "_" + countries.front() + "_" + countries.back() + "_" + areaSuffix;
            themeParameters->setParameter(MATCHED_STANDING_WATER_TABLE, ign::data::String(matchedStandingWaterTableName));
        }
        
        //tables de travail
        if ( !suffix.empty() ) {
            std::string tableBaseName = themeParameters->getValue(EDGE_TABLE_INIT_BASE).toString();
            std::string tableName = tableBaseName + "_" + countries.front() + "_" + countries.back() + "_" + suffix;
            themeParameters->setParameter(EDGE_TABLE_INIT, ign::data::String(tableName));
        }
        if ( themeParameters->getValue(CL_TABLE).toString() == "" )
            themeParameters->setParameter(CL_TABLE, ign::data::String(themeParameters->getValue(EDGE_TABLE_INIT).toString() + themeParameters->getValue(CL_TABLE_SUFFIX).toString()));
        if ( themeParameters->getValue(CP_TABLE).toString() == "" ) 
            themeParameters->setParameter(CP_TABLE, ign::data::String(themeParameters->getValue(EDGE_TABLE_INIT).toString() + themeParameters->getValue(CP_TABLE_SUFFIX).toString()));

        //set BDD search path
        context->getDataBaseManager().setSearchPath(themeParameters->getValue(WORKING_SCHEMA).toString());
        ome2::utils::setTableName<app::params::ThemeParametersS>(LANDMASK_TABLE);
        ome2::utils::setTableName<app::params::ThemeParametersS>(BOUNDARY_SMOOTHED_TABLE);
        ome2::utils::setTableName<app::params::ThemeParametersS>(ALL_WATERCOURSE_AREA_TABLE);
        ome2::utils::setTableName<app::params::ThemeParametersS>(ALL_STANDING_WATER_TABLE);
        ome2::utils::setTableName<app::params::ThemeParametersS>(ALL_EDGE_TABLE);
        ome2::utils::setTableName<epg::params::EpgParametersS>(TARGET_BOUNDARY_TABLE);

        //créer les tables CP et CL vides si elles n'existent pas
        app::utils::createCpClTables(themeParameters->getValue(EDGE_TABLE_INIT).toString());

        logger->log(epg::log::INFO, "[START EDGE MATCHING PROCESS ] " + epg::tools::TimeTools::getTime());

        //lancement du traitement
        stepSuitePtr->run(stepCode, verbose);

		logger->log(epg::log::INFO, "[END EDGE MATCHING PROCESS ] " + epg::tools::TimeTools::getTime());
    }
    catch( ign::Exception &e )
    {
        std::cerr << e.diagnostic() << std::endl;
        epg::log::EpgLoggerS::getInstance()->log( epg::log::ERROR, std::string(e.diagnostic()));
        logFile << e.diagnostic() << std::endl;
        returnValue = 1;
    }
    catch( std::exception &e )
    {
        std::cerr << e.what() << std::endl;
        epg::log::EpgLoggerS::getInstance()->log( epg::log::ERROR, std::string(e.what()));
        logFile << e.what() << std::endl;
        returnValue = 1;
    }

    logFile << "[END] " << epg::tools::TimeTools::getTime() << std::endl;

    epg::ContextS::kill();
    epg::log::EpgLoggerS::kill();
    epg::log::ShapeLoggerS::kill();
    epg::params::EpgParametersS::kill();
    app::params::ThemeParametersS::kill();
    
    logFile.close();

    return returnValue ;
}