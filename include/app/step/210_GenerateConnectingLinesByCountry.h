#ifndef _APP_STEP_GENERATECONNECTINGLINESBYCOUNTRY_H_
#define _APP_STEP_GENERATECONNECTINGLINESBYCOUNTRY_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class GenerateConnectingLinesByCountry : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 210; };

		/// \brief
		std::string getName() { return "GenerateConnectingLinesByCountry"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif