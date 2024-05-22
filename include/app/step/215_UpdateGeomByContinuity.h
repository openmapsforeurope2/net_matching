#ifndef _APP_STEP_UPDATEGEOMBYCONTINUITY_H_
#define _APP_STEP_UPDATEGEOMBYCONTINUITY_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class UpdateGeomByContinuity : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 215; };

		/// \brief
		std::string getName() { return "UpdateGeomByContinuity"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif