// APP
#include <app/geometry/tools/getBuffer.h>

// EPG
#include <epg/Context.h>
#include <epg/sql/DataBaseManager.h>
#include <epg/log/EpgLogger.h>

// SOCLE
#include <ign/geometry/io/WkbReader.h>


namespace app{
namespace geometry{
namespace tools{

    ///
    ///
    ///
    ign::geometry::Geometry* getBuffer(
        ign::geometry::Geometry const& geom,
        double distBuffer,
        size_t quadrantSegments,
		ign::geometry::algorithm::BufferOpGeos::EndCapStyle endCapStyle
    ) {
        ign::geometry::Geometry* result = 0;
        ign::geometry::algorithm::BufferOpGeos buffOp;
        try {
            result = buffOp.buffer(geom, distBuffer, quadrantSegments, endCapStyle);
        } catch (...) {
            epg::log::EpgLoggerS::getInstance()->log(epg::log::WARN, "ign::geometry::algorithm::BufferOpGeos failed");

            // context
            epg::Context* context = epg::ContextS::getInstance();

            //--
            std::ostringstream ss;
            ss  << "SELECT ST_Buffer(ST_GeomFromText('" + geom.toString() + "'), "
                << distBuffer << ", ";

            if( quadrantSegments != 0 )
                ss << quadrantSegments << ")";
            else {
                std::string style = "'endcap=";
                switch ( endCapStyle )
                {
		            case ign::geometry::algorithm::BufferOpGeos::CAP_ROUND :{
			            style += "round";
                        break;
                    }
                    case ign::geometry::algorithm::BufferOpGeos::CAP_FLAT :{
			            style += "flat";
                        break;
                    }
                    case ign::geometry::algorithm::BufferOpGeos::CAP_SQUARE :{
			            style += "square";
                        break;
                    }
			    }
                ss << style << " join=round')";
            }

            //--
            ign::sql::SqlResultSetPtr resultPtr = context->getDataBaseManager().getConnection()->query( ss.str() );

            //--
            result = ign::geometry::io::ReadHexaWkb( resultPtr->getFieldValue(0,0).toString() );

        }
        return result;
    }

}
}
}