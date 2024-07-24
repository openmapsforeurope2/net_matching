#ifndef _APP_STEP_GENERATECLINSINAREA_H_
#define _APP_STEP_GENERATECLINSINAREA_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app {
	namespace step {

		class GenerateCLinArea : public epg::step::StepBase< app::params::ThemeParametersS > {

		public:

			/// \brief
			int getCode() { return 204; };

			/// \brief
			std::string getName() { return "GenerateCLinArea"; };

			/// \brief
			void onCompute(bool);

			/// \brief
			void init();

		};

	}
}

#endif