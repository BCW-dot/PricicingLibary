#include "CallOption.h"
#include "matlib.h"


double CallOption::payoff( double stockAtMaturity ) const {
    if (stockAtMaturity>getStrike()) {
        return stockAtMaturity-getStrike();
    } else {
        return 0.0;
    }
}

double CallOption::price( const BlackScholesModel& bsm ) const {
    double S = bsm.getStockPrice();
    double K = getStrike();
    double sigma = bsm.getVolatility();
    double r = bsm.getRiskFreeRate();
    double T = getMaturity() - bsm.getDate();

    double numerator = log( S/K ) + ( r + sigma*sigma*0.5)*T;
    double denominator = sigma * sqrt(T );
    double d1 = numerator/denominator;
    double d2 = d1 - denominator;
    return S*normcdf(d1) - exp(-r*T)*K*normcdf(d2);
}

double CallOption::delta( const BlackScholesModel& bsm ) const {
	double S = bsm.getStockPrice();
	double K = getStrike();
	double sigma = bsm.getVolatility();
	double r = bsm.getRiskFreeRate();
	double T = getMaturity() - bsm.getDate();

	double numerator=log(S/K)+(r+sigma*sigma*0.5)*T;
	double denominator = sigma * sqrt(T);
	double d1 = numerator / denominator;
	return normcdf(d1);
}





//////////////////////////
//
//  Test the call option class
//  
//
//////////////////////////

static void testCallOptionPrice() {
    CallOption callOption;
    callOption.setStrike( 105.0 );
    callOption.setMaturity( 2.0 );
    
    BlackScholesModel bsm;
    bsm.setRiskFreeRate(0.05);
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100.0);
    bsm.setDate(1.0);

    double price = callOption.price( bsm );
    ASSERT_APPROX_EQUAL( price, 4.046, 0.01);
}

static void testCallOptionDelta() {
	CallOption callOption;
	callOption.setStrike(105.0);
	callOption.setMaturity(2.0);

    BlackScholesModel bsm;
    bsm.setRiskFreeRate(0.05);
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100.0);
    bsm.setDate(1.0);

	double d = callOption.delta(bsm);
	ASSERT_APPROX_EQUAL(d, 0.524, 0.001);
}

void testCallOption() {
    TEST( testCallOptionPrice );
	TEST(testCallOptionDelta);
}
