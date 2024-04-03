#include "HedgingSimulator.h"
#include "testing.h"
#include "matlib.h"

using namespace std;

HedgingSimulator::HedgingSimulator() {
	// Choose default models and options
	shared_ptr<BlackScholesModel> model(new BlackScholesModel());
    model->setStockPrice(1);
    model->setDate(0);
    model->setRiskFreeRate(0.05);
    model->setVolatility(0.2);
    model->setDrift(0.1);

	shared_ptr<CallOption> option = make_shared<CallOption>();
	option->setStrike(model->getStockPrice());
	option->setMaturity(1);

	setToHedge(option);
    setPricingModel(model);
	setSimulationModel(model);
	nSteps = 10;
}
/*  Runs a number of simulations and returns
a vector of the profit and loss */
std::vector<double> HedgingSimulator::runSimulations(int nSimulations) const {
	std::vector<double> ret(nSimulations);
	for (int i = 0; i < nSimulations; i++) {
		ret[i] = runSimulation();
	}
	return ret;
}
/*  Run a simulation and compute the profit and loss */
double HedgingSimulator::runSimulation() const {
    //i = 0 in developed theory/equations
	double T = toHedge->getMaturity();
	double S0 = simulationModel->getStockPrice();
    vector<double> pricePath = simulationModel->generatePricePath(T,nSteps);//generateRiskNeutralPricePath(T, nSteps);
	double dt = T / nSteps;
	double charge = chooseCharge(S0); //toHedge->price(*pricingModel);
    double stockQuantity = selectStockQuantity(0,S0);
	double bankBalance = charge - stockQuantity*S0; //charge is the BSM price for the call
    
    //i = N-1
	for (int i = 0; i< nSteps-1; i++) {
		double balanceWithInterest = bankBalance * exp(simulationModel->getRiskFreeRate()*dt);
		double S = pricePath[i];
		double date = dt*(i + 1);
		double newStockQuantity = selectStockQuantity(date, S);
		double costs = (newStockQuantity - stockQuantity)*S;
        bankBalance = balanceWithInterest - costs;
		stockQuantity = newStockQuantity;
	}
    
    //i = N
	double balanceWithInterest = bankBalance * exp(simulationModel->getRiskFreeRate()*dt);
	double S = pricePath[nSteps - 1];
	double stockValue = stockQuantity*S;
	double payout = toHedge->payoff(pricePath);
	return balanceWithInterest + stockValue - payout;
}

/*  How much should we charge the customer */
double HedgingSimulator::chooseCharge(double stockPrice) const {
	// create a copy of the pricing model
	BlackScholesModel pm = *pricingModel;
    pm.setStockPrice(stockPrice);
	return toHedge->price(pm);
}
/*  How much stock should we hold */
double HedgingSimulator::selectStockQuantity(double date, double stockPrice) const {
	// create a copy of the pricing model
	BlackScholesModel pm = *pricingModel;
    pm.setStockPrice(stockPrice);
    pm.setDate(date);
	return toHedge->delta(pm);
}


//
//
//   Tests
//
//


static void testDeltaHedgingMeanPayoff() {
  rng("default");
  HedgingSimulator simulator;
  simulator.setNSteps(1000);
  vector<double> result = simulator.runSimulations(10000);
  ASSERT_APPROX_EQUAL(result[0], 0.0, 0.01);
  hist("deltaHedgingFrequentlyPNL.html", result, 20);
}

static void testPlotDeltaHedgingHistogram() {
	rng("default");
	HedgingSimulator simulator;
	simulator.setNSteps(100);
	vector<double> result = simulator.runSimulations(10000);
	hist("deltaHedgingInfrequentlyPNL.html", result, 20);
}
void testHedgingSimulator() {
	TEST(testDeltaHedgingMeanPayoff);
	TEST(testPlotDeltaHedgingHistogram);
}
