#ifndef _APP_STEP_CONNECTIONCONNECTINGPOINT_H_
#define _APP_STEP_CONNECTIONCONNECTINGPOINT_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class ConnectionConnectingPoint : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 250; };

		/// \brief
		std::string getName() { return "ConnectionConnectingPoint"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif