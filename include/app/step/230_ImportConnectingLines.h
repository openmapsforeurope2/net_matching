#ifndef _APP_STEP_IMPORTCONNECTINGLINES_H_
#define _APP_STEP_IMPORTCONNECTINGLINES_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class ImportConnectingLines : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 230; };

		/// \brief
		std::string getName() { return "ImportConnectingLines"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif