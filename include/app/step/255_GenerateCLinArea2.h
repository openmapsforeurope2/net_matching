#ifndef _APP_STEP_GENERATECLINSINAREA2_H_
#define _APP_STEP_GENERATECLINSINAREA2_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app {
	namespace step {

		class GenerateCLinArea2 : public epg::step::StepBase< app::params::ThemeParametersS > {

		public:

			/// \brief
			int getCode() { return 22554; };

			/// \brief
			std::string getName() { return "GenerateCLinArea2"; };

			/// \brief
			void onCompute(bool);

			/// \brief
			void init();

		};

	}
}

#endif