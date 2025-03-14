
//BOOST
#include <boost/program_options.hpp>

//EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/params/tools/loadParameters.h>

//OME2
#include <ome2/utils/setTableName.h>

//APP
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
    std::string     stepCode = "";
    std::string     countryCode = "";
    std::string     theme = "";
    bool            verbose = true;


    epg::step::StepSuite< app::params::ThemeParametersS > stepSuiteTn, stepSuiteHy, stepSuiteRa;
	app::step::tools::initStepsHy(stepSuiteHy);
	app::step::tools::initStepsTn(stepSuiteTn);
	app::step::tools::initStepsRa(stepSuiteRa);


	std::ostringstream OperatorDetail;
	OperatorDetail << "set step :" << std::endl
        << "TN:" << std::endl
		<< stepSuiteTn.toString() 
		<< "HY:" << std::endl
		<< stepSuiteHy.toString()
        << "RA:" << std::endl
		<< stepSuiteRa.toString();


    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("c" , po::value< std::string >(&epgParametersFile)     , "conf file" )
        ("T" , po::value< std::string >(&theme)                 , "theme" )
        ("cc" , po::value< std::string >(&countryCode)          , "country code" )
        ("sp", po::value< std::string >(&stepCode), OperatorDetail.str().c_str())
    ;

    //main log
    std::string     logFileName = "log.txt";
    std::ofstream   logFile( logFileName.c_str() ) ;

    logFile << "[START] " << epg::tools::TimeTools::getTime() << std::endl;

    int returnValue = 0;
    try{

        po::variables_map vm;
        int style = po::command_line_style::default_style & ~po::command_line_style::allow_guessing;
        po::store( po::parse_command_line( argc, argv, desc, style ), vm );
        po::notify( vm );    

        if ( vm.count( "help" ) ) {
            std::cout << desc << std::endl;
            return 1;
        }

        epg::step::StepSuite< app::params::ThemeParametersS >* stepSuitePtr = 0;
        if (theme == "hy")
			stepSuitePtr = &stepSuiteHy;
		else if (theme == "tn")
			stepSuitePtr = &stepSuiteTn;
        else if (theme == "ra")
			stepSuitePtr = &stepSuiteRa;
        else {
            std::string mError = "unknown theme " + theme;
            IGN_THROW_EXCEPTION(mError);
        }

        if (stepCode.empty()) stepCode = stepSuitePtr->getStepsRange();

        //parametres EPG
		context->loadEpgParameters( epgParametersFile, theme );

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
        // logger->setProdOfstream( logDirectory+"/au_merging.log" );
        logger->setDevOfstream( context->getLogDirectory()+"/net_matching.log" );

        //shape logger
        epg::log::ShapeLogger* shapeLogger = epg::log::ShapeLoggerS::getInstance();
	    shapeLogger->setDataDirectory( context->getLogDirectory()+"/shape" );
        
        //theme parameters
        themeParametersFile = context->getConfigParameters().getValue( THEME_PARAMETER_FILE ).toString();
		app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
        epg::params::tools::loadParams( *themeParameters, themeParametersFile, countryCode );
        if ( themeParameters->getValue(CL_TABLE).toString() == "" )
            themeParameters->setParameter(CL_TABLE, ign::data::String(themeParameters->getValue(EDGE_TABLE_INIT).toString() + themeParameters->getValue(CL_TABLE_SUFFIX).toString()));
        if ( themeParameters->getValue(CP_TABLE).toString() == "" ) 
            themeParameters->setParameter(CP_TABLE, ign::data::String(themeParameters->getValue(EDGE_TABLE_INIT).toString() + themeParameters->getValue(CP_TABLE_SUFFIX).toString()));

        //info de connection db
        context->loadEpgParameters( themeParameters->getValue(DB_CONF_FILE).toString() );

        //set BDD search path
        context->getDataBaseManager().setSearchPath(themeParameters->getValue(WORKING_SCHEMA).toString());
        ome2::utils::setTableName<app::params::ThemeParametersS>(LANDMASK_TABLE);
        ome2::utils::setTableName<epg::params::EpgParametersS>(TARGET_BOUNDARY_TABLE);

        //créer les tables CP et CL vides si elles n'existent pas
        app::utils::createCpClTables(themeParameters->getValue(EDGE_TABLE_INIT).toString());

        logger->log(epg::log::INFO, "[START EDGE-MATCHING PROCESS ] " + epg::tools::TimeTools::getTime());

        //lancement du traitement
        stepSuitePtr->run(stepCode, verbose);

		logger->log(epg::log::INFO, "[END EDGE-MATCHING PROCESS ] " + epg::tools::TimeTools::getTime());

    }
    catch( ign::Exception &e )
    {
        std::cerr<< e.diagnostic() << std::endl;
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
    app::params::ThemeParametersS::kill();
    
    logFile.close();

    return returnValue ;
}