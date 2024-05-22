#ifndef _APP_STEP_GENERATECONNECTINGPOINT_H_
#define _APP_STEP_GENERATECONNECTINGPOINT_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class GenerateConnectingPoint : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 240; };

		/// \brief
		std::string getName() { return "GenerateConnectingPoint"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif