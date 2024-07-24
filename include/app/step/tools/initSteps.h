#ifndef _APP_STEP_TOOLS_INITSTEPS_H_
#define _APP_STEP_TOOLS_INITSTEPS_H_

//EPG
#include <epg/step/StepSuite.h>
#include <epg/step/factoryNew.h>

//APP
#include <app/step/202_JunctionMatching.h>
#include <app/step/204_GenerateCLinArea.h>
#include <app/step/210_GenerateConnectingLinesByCountry.h>
#include <app/step/211_MergeConnectingLinesOnBorder.h>
#include <app/step/212_SnapConnectingLines.h>
#include <app/step/213_DeleteConnectingLines.h>
#include <app/step/214_UpdateGeomConnectingLines.h>
#include <app/step/215_UpdateGeomByContinuity.h>
#include <app/step/220_ConnectionConnectingLines.h>
#include <app/step/230_ImportConnectingLines.h>
#include <app/step/240_GenerateConnectingPoint.h>
#include <app/step/250_ConnectionConnectingPoint.h>
#include <app/step/260_EdgeCleaning1.h>
#include <app/step/270_EdgeConnection.h>
#include <app/step/280_EdgeCleaning2.h>


namespace app{
namespace step{
namespace tools{

	template<  typename StepSuiteType >
	void initStepsHy( StepSuiteType& stepSuite )
	{
		// stepSuite.addStep( epg::step::factoryNew< JunctionMatching >() );
		stepSuite.addStep( epg::step::factoryNew< GenerateCLinArea >() );
		stepSuite.addStep( epg::step::factoryNew< GenerateConnectingLinesByCountry >() );
		stepSuite.addStep( epg::step::factoryNew< MergeConnectingLinesOnBorder >() );
		stepSuite.addStep( epg::step::factoryNew< SnapConnectingLines >() );
		stepSuite.addStep( epg::step::factoryNew< DeleteConnectingLines >() );
		stepSuite.addStep( epg::step::factoryNew< UpdateGeomConnectingLines >() );
		stepSuite.addStep( epg::step::factoryNew< UpdateGeomByContinuity >() );
		stepSuite.addStep( epg::step::factoryNew< ConnectionConnectingLines >() );
		stepSuite.addStep( epg::step::factoryNew< ImportConnectingLines >() );
		stepSuite.addStep( epg::step::factoryNew< GenerateConnectingPoint >() );
		stepSuite.addStep( epg::step::factoryNew< ConnectionConnectingPoint >() );
		stepSuite.addStep( epg::step::factoryNew< EdgeCleaning1 >() );
		stepSuite.addStep( epg::step::factoryNew< EdgeConnection >() );
		stepSuite.addStep( epg::step::factoryNew< EdgeCleaning2 >() );
	}

	template<  typename StepSuiteType >
	void initStepsTn(StepSuiteType& stepSuite)
	{
		stepSuite.addStep(epg::step::factoryNew< GenerateConnectingLinesByCountry >());
		stepSuite.addStep(epg::step::factoryNew< MergeConnectingLinesOnBorder >());
		stepSuite.addStep(epg::step::factoryNew< SnapConnectingLines >());
		stepSuite.addStep(epg::step::factoryNew< DeleteConnectingLines >());
		stepSuite.addStep(epg::step::factoryNew< UpdateGeomConnectingLines >());
		stepSuite.addStep(epg::step::factoryNew< UpdateGeomByContinuity >());
		stepSuite.addStep(epg::step::factoryNew< ConnectionConnectingLines >());
		stepSuite.addStep(epg::step::factoryNew< ImportConnectingLines >());
		stepSuite.addStep(epg::step::factoryNew< GenerateConnectingPoint >());
		stepSuite.addStep(epg::step::factoryNew< ConnectionConnectingPoint >());
		stepSuite.addStep(epg::step::factoryNew< EdgeCleaning1 >());
		stepSuite.addStep(epg::step::factoryNew< EdgeConnection >());
		stepSuite.addStep(epg::step::factoryNew< EdgeCleaning2 >());
	}

}
}
}

#endif