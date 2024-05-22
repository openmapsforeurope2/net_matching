#ifndef _APP_STEP_UPDATEGEOMCONNECTINLINES_H_
#define _APP_STEP_UPDATEGEOMCONNECTINLINES_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class UpdateGeomConnectinLines : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 214; };

		/// \brief
		std::string getName() { return "UpdateGeomConnectinLines"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif