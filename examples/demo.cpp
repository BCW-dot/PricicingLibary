#include "PricingLibrary/utils/matlib.h"
#include "PricingLibrary/utils/geometry.h"
#include "PricingLibrary/utils/textfunctions.h"
#include "PricingLibrary/options/CallOption.h"
#include "PricingLibrary/options/PutOption.h"
#include "PricingLibrary/charts/PieChart.h"
#include "PricingLibrary/charts/LineChart.h"
#include "PricingLibrary/models/BlackScholesModel.h"
#include "PricingLibrary/charts/Histogram.h"
#include "PricingLibrary/pricing/MonteCarloPricer.h"
#include "PricingLibrary/options/UpAndOutOption.h"
#include "PricingLibrary/options/UpAndInOption.h"
#include "PricingLibrary/options/DownAndOutOption.h"
#include "PricingLibrary/portfolio/Portfolio.h"
#include "PricingLibrary/portfolio/HedgingSimulator.h"
#include "PricingLibrary/options/AsianOption.h"
#include "PricingLibrary/models/TwoLevelModel.h"

#include "PricingLibrary/pde/pde.h"
#include "PricingLibrary/pde/fdm.h"

//using namespace std;

int main() {
    /*
    testMatlib();
    testGeometry();
    testPieChart();
    testCallOption();
    testPutOption();
    testBlackScholesModel();
    testLineChart();
    testTextFunctions();
    testHistogram();
    testMonteCarloPricer();
    testUpAndOutOption();
	testUpAndInOption();
    testDownAndOutOption();
	testAsianOption();
    testContinuousTimeOptionBase();
    testPortfolio();
	testHedgingSimulator();
	*/
    //testMonteCarloPricer();
    //testHedgingSimulator();
    //testTwoLevelModel();
    testBlackScholesPDE();
    return 0;
}
