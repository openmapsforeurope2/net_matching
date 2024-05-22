#ifndef _APP_STEP_EDGECONNECTION_H_
#define _APP_STEP_EDGECONNECTION_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class EdgeConnection : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 270; };

		/// \brief
		std::string getName() { return "EdgeConnection"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif