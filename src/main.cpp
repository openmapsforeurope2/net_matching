
//BOOST
#include <boost/program_options.hpp>

//EPG
#include <epg/Context.h>
#include <epg/tools/TimeTools.h>
#include <epg/params/tools/loadParameters.h>

//APP
#include <app/params/ThemeParameters.h>
#include <app/calcul/CFeatConnectionOp.h>
#include <app/calcul/CFeatGenerationOp.h>
#include <app/calcul/EdgeCleaningOp.h>
#include <app/calcul/EdgeConnectorOp.h>


namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    // ign::geometry::PrecisionModel::SetDefaultPrecisionModel(ign::geometry::PrecisionModel(ign::geometry::PrecisionModel::FIXED, 1.0e5, 1.0e7) );
    epg::Context* context = epg::ContextS::getInstance();

    std::string     logDirectory = "";
    std::string     epgParametersFile = "";
    std::string     themeParametersFile = "";
    std::string     countryCode = "";
    bool            verbose = true;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("c" , po::value< std::string >(&epgParametersFile)     , "conf file" )
        ("cc" , po::value< std::string >(&countryCode)          , "country code" )
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

        //parametres EPG
		context->loadEpgParameters( epgParametersFile );

        //Initialisation du log de prod
        logDirectory = context->getConfigParameters().getValue( LOG_DIRECTORY ).toString();

        //test si le dossier de log existe sinon le creer
        boost::filesystem::path logDir(logDirectory);
        if (!boost::filesystem::is_directory(logDir))
        {
            if (!boost::filesystem::create_directory(logDir))
            {
                std::string mError = "le dossier " + logDirectory + " ne peut être cree";
                IGN_THROW_EXCEPTION(mError);
            }
        }

        //repertoire de travail
        context->setLogDirectory( logDirectory );

        //epg logger
        epg::log::EpgLogger* logger = epg::log::EpgLoggerS::getInstance();
        // logger->setProdOfstream( logDirectory+"/au_merging.log" );
        logger->setDevOfstream( context->getLogDirectory()+"/tn_matching.log" );

        //shape logger
        epg::log::ShapeLogger* shapeLogger = epg::log::ShapeLoggerS::getInstance();
	    shapeLogger->setDataDirectory( context->getLogDirectory()+"/shape" );
        
		logger->log(epg::log::INFO, "[START EDGE-MATCHING PROCESS ] " + epg::tools::TimeTools::getTime());
        
        //theme parameters
        themeParametersFile = context->getConfigParameters().getValue( THEME_PARAMETER_FILE ).toString();
		app::params::ThemeParameters* themeParameters = app::params::ThemeParametersS::getInstance();
        epg::params::tools::loadParams( *themeParameters, themeParametersFile );

        //lancement du traitement
		bool clCompute = themeParameters->getValue(CL_COMPUTE).toBoolean();
		
		app::calcul::CFeatGenerationOp cFeatGenerationOp(countryCode);
        app::calcul::CFeatConnectionOp cFeatConnectionOp(countryCode, verbose);
		
		if (clCompute)
		{
			cFeatGenerationOp.computeCL();
            cFeatConnectionOp.computeCl();
            cFeatConnectionOp.computeClImport();
		}

		cFeatGenerationOp.computeCP();
        cFeatConnectionOp.computeCp();

        // nettoyage
		app::calcul::EdgeCleaningOp edgeCleaningOp(countryCode, verbose);
        std::set<std::string> sTreatedFeatures;
        edgeCleaningOp.cleanFaces();
        // edgeCleaningOp.cleanPathsOutOfCountry();
        edgeCleaningOp.cleanParalelleEdges();
        edgeCleaningOp.cleanFacesAndAntennaByCountry(sTreatedFeatures);

        app::calcul::EdgeConnectorOp::compute(countryCode, verbose);

        // nettoyage
        edgeCleaningOp.cleanParalelleEdges();
        edgeCleaningOp.cleanFaces2(ign::feature::FeatureFilter("form_of_way <> 'bicycle_road'"));

        // app::calcul::EdgeConnectorOp::compute(countryCode, verbose);

        // edgeCleaningOp.cleanFaces2(ign::feature::FeatureFilter("form_of_way <> 'bicycle_road'"));

        // edgeCleaningOp.cleanAntennas(sTreatedFeatures);

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