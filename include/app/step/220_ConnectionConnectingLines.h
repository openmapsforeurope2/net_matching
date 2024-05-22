#ifndef _APP_STEP_CONNECTIONCONNECTINGLINES_H_
#define _APP_STEP_CONNECTIONCONNECTINGLINES_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class ConnectionConnectingLines : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 220; };

		/// \brief
		std::string getName() { return "ConnectionConnectingLines"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif