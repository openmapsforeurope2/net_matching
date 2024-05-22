#ifndef _APP_STEP_EDGECLEANING1_H_
#define _APP_STEP_EDGECLEANING1_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class EdgeCleaning1 : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 260; };

		/// \brief
		std::string getName() { return "EdgeCleaning1"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif