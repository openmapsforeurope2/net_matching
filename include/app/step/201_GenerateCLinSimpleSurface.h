#ifndef _APP_STEP_GENERATECLINSIMPLESURFACE_H_
#define _APP_STEP_GENERATECLINSIMPLESURFACE_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app {
	namespace step {

		class GenerateCLinSimpleSurface : public epg::step::StepBase< app::params::ThemeParametersS > {

		public:

			/// \brief
			int getCode() { return 201; };

			/// \brief
			std::string getName() { return "GenerateCLinSimpleSurface"; };

			/// \brief
			void onCompute(bool);

			/// \brief
			void init();

		};

	}
}

#endif