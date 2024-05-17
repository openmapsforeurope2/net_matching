#ifndef _APP_STEP_DELETECONNECTINGLINES_H_
#define _APP_STEP_DELETECONNECTINGLINES_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class DeleteConnectingLines : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 213; };

		/// \brief
		std::string getName() { return "DeleteConnectingLines"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif