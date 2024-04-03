#include "AsianOption.h"

#include "matlib.h"
#include "MonteCarloPricer.h"
#include "CallOption.h"

double AsianOption::payoff(const std::vector<double>& stockPrices) const {
	double m = mean(stockPrices);
	if (m > getStrike()) {
		return m - getStrike();
	}
	else {
		return 0;
	}
}

static void testMonteCarloPriceAsian() {
	rng("default");
	AsianOption option;
	option.setMaturity(1.0);
	option.setStrike(100.0);
    
	BlackScholesModel model;
	model.setStockPrice(95);
    model.setVolatility(0.1);
    model.setDrift(0);
    model.setRiskFreeRate(0);
    
	MonteCarloPricer pricer;
	pricer.nSteps = 100;
	pricer.nScenarios = 3;
	double price = pricer.price(option, model);

	CallOption callOption;
	callOption.setMaturity(1.0);
	callOption.setStrike(100.0);
	double europeanPrice = callOption.price(model);
	// without studying any theory on asian options we
	// can't guess much about the price, but it is presumably
	// not too wildly different from the european option price
	assert(price > 0);
	assert(price < 2*europeanPrice);
}

void testAsianOption() {
	TEST(testMonteCarloPriceAsian);
}
