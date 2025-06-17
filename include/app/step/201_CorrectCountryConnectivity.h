#ifndef _APP_STEP_CORRECTCOUNTRYCONNECTIVITY_H_
#define _APP_STEP_CORRECTCOUNTRYCONNECTIVITY_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app {
	namespace step {

		class CorrectCountryConnectivity : public epg::step::StepBase< app::params::ThemeParametersS > {

		public:

			/// \brief
			int getCode() { return 201; };

			/// \brief
			std::string getName() { return "CorrectCountryConnectivity"; };

			/// \brief
			void onCompute(bool);

			/// \brief
			void init();

		};

	}
}

#endif