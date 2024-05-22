#ifndef _APP_STEP_SNAPCONNECTINGLINES_H_
#define _APP_STEP_SNAPCONNECTINGLINES_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class SnapConnectingLines : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 212; };

		/// \brief
		std::string getName() { return "SnapConnectingLines"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif