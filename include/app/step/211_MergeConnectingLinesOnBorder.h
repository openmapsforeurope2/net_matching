#ifndef _APP_STEP_MERGECONNECTINGLINESONBORDER_H_
#define _APP_STEP_MERGECONNECTINGLINESONBORDER_H_

#include <epg/step/StepBase.h>
#include <app/params/ThemeParameters.h>

namespace app{
namespace step{

	class MergeConnectingLinesOnBorder : public epg::step::StepBase< app::params::ThemeParametersS > {

	public:

		/// \brief
		int getCode() { return 211; };

		/// \brief
		std::string getName() { return "MergeConnectingLinesOnBorder"; };

		/// \brief
		void onCompute( bool );

		/// \brief
		void init();

	};

}
}

#endif