#include "matlib.h"
#include "geometry.h"
#include "textfunctions.h"
#include "CallOption.h"
#include "PutOption.h"
#include "PieChart.h"
#include "LineChart.h"
#include "BlackScholesModel.h"
#include "Histogram.h"
#include "MonteCarloPricer.h"
#include "UpAndOutOption.h"
#include "UpAndInOption.h"
#include "DownAndOutOption.h"
#include "Portfolio.h"
#include "HedgingSimulator.h"
#include "AsianOption.h"
#include "TwoLevelModel.h"

#include "pde.h"
#include "fdm.h"
#include "HestonModel.h"
#include "ScottChesney.h"

#include <ctime>

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
    //testBlackScholesModel();
    //testMonteCarloPricer();
    //testHedgingSimulator();
    //testTwoLevelModel();
    //testBlackScholesPDE();
    //testMatlib();
    //testHestonModel();
    testScottChesneyModel();
    return 0;
}
