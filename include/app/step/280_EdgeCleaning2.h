#ifndef _APP_STEP_EDGECLEANING2_H_
#define _APP_STEP_EDGECLEANING2_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class EdgeCleaning2 : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 280; };

		/// \brief
		std::string getName() { return "EdgeCleaning2"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif