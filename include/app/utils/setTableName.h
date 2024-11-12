#ifndef _APP_UTILS_SETTABLENAME_H_
#define _APP_UTILS_SETTABLENAME_H_

#include <string>
#include <app/params/ThemeParameters.h>
#include <epg/params/EpgParameters.h>

namespace app{
namespace utils{

    //--
    void setTableName(TN_PARAMETERS TABLE_NAME_PARAM);

    //--
    void setTableName(EPG_PARAMETERS TABLE_NAME_PARAM);
}
}

#endif