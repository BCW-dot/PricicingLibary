#include "MonteCarloPricer.h"
#include "BlackScholesModel.h"
#include "StockPriceModel.h"

#include "matlib.h"
#include "CallOption.h"
#include "PutOption.h"
#include "HestonModel.h"
#include "ScottChesney.h"

using namespace std;

MonteCarloPricer::MonteCarloPricer() :
    nScenarios(10000),
    nSteps(100) {
}

double MonteCarloPricer::price(const ContinuousTimeOption& option, const StockPriceModel& model ) {
    int nSteps = this->nSteps;
    /*
    if (!option.isPathDependent()) {
        nSteps = 1;
    }
     */
    double total = 0.0;
    for (int i=0; i<nScenarios; i++) {
        vector<double> path = model.generateRiskNeutralPricePath(option.getMaturity(), nSteps );
        double payoff = option.payoff( path );
        total+= payoff;
    }
    double mean = total/nScenarios;
    double r = model.getRiskFreeRate();
    double T = option.getMaturity() - model.getDate();
    return exp(-r*T)*mean;
}

double MonteCarloPricer::delta(const ContinuousTimeOption& option, const BlackScholesModel &model){
    int nSteps = this->nSteps;
    if (!option.isPathDependent()) {
        nSteps = 1;
    }
    
    double h = model.getStockPrice() * pow(10,-6);
    double stock_price_right = model.getStockPrice() + h;
    double stock_price_left = model.getStockPrice() - h;
    
    double total_right = 0.0;
    double total_left = 0.0;
    
    for(int i = 0; i<nScenarios; i++){
        vector<double> epsilon = randn(nSteps);
        
        /* S_0+h payoff*/
        vector<double> path = model.generateRiskNeutralDeltaPath(option.getMaturity(), nSteps, stock_price_right, epsilon);
        double payoff_right = option.payoff(path);
        total_right += payoff_right;
        
        /* S_0-h payoff*/
        path = model.generateRiskNeutralDeltaPath(option.getMaturity(), nSteps, stock_price_left, epsilon);
        double payoff_left = option.payoff(path);
        total_left += payoff_left;
    }
    double mean_right = total_right/nScenarios;
    double mean_left = total_left/nScenarios;
    
    double r = model.getRiskFreeRate();
    double T = option.getMaturity() - model.getDate();
    
    return exp(-r*T)*(mean_right - mean_left)/(2*h);
}

//////////////////////////////////////
//
//   Tests
//
//////////////////////////////////////

static void testPriceCallOption() {
    my_rng();

    CallOption c;
    c.setStrike( 110 );
    c.setMaturity( 2 );

    BlackScholesModel m;
    m.setRiskFreeRate(0.05);
    m.setVolatility(0.1);
    m.setStockPrice(100.0);
    m.setDate(1.0);
    m.setDrift(0.1);
    
    MonteCarloPricer pricer;
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    double price = pricer.price( c, m );
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << " Option were priced in "<< pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
    double expected = c.price( m );
    ASSERT_APPROX_EQUAL( price, expected, 0.1 );
}

static void testHestonCallOption(){
    my_rng();
    
    double S_0 = 100.0;    // Initial spot price
    double r = 0.03;     // Risk-free rate
    double v_0 = 0.01; // Initial volatility
    double T = 1.00;       // One year until expiry
    double K = 100;
    double date_0 = 0.0;

    double rho = 0.8;     // Correlation of asset and volatility
    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.019;  // Long run average volatility
    double xi = 0.61;      // "Vol of vol"

    HestonModel hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(date_0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;      // "Vol of vol"

    // Create the PayOff, Option and Heston objects
    CallOption c;
    c.setStrike(K);
    c.setMaturity(T);
    
    BlackScholesModel bsm;
    bsm.setRiskFreeRate(r);
    bsm.setVolatility(v_0);
    bsm.setStockPrice(S_0);
    bsm.setDate(date_0);
    
    MonteCarloPricer pricer;
    pricer.nSteps = 1000;
    std::cout << "Monte Carlo pricer Heston: " << pricer.price(c, hest_euler) << std::endl;
    std::cout << "Closed form: " << c.price(bsm) << std::endl;
    
    //ASSERT_APPROX_EQUAL(pricer.price(c, hest_euler), c.price(bsm), 0.1);
}

static void testHestonVsScott(){
    my_rng();
    
    double S_0 = 100.0;    // Initial spot price
    double r = 0.03;     // Risk-free rate
    double v_0 = 0.01; // Initial volatility
    double T = 1.00;       // One year until expiry
    double date_0 = 0.0;
    double K = 100;

    double rho = 0.8;     // Correlation of asset and volatility
    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.019;  // Long run average volatility
    double xi = 0.61;      // "Vol of vol"

    HestonModel hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(date_0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;
    
    ScottChesney scott_euler;
    scott_euler.setRiskFreeRate(r);
    scott_euler.setVolatility(v_0);
    scott_euler.setStockPrice(S_0);
    scott_euler.setDate(date_0);
    scott_euler.rho = rho;
    scott_euler.kappa = kappa;
    scott_euler.theta = theta;
    scott_euler.xi = xi;

    // Create the PayOff, Option and Heston objects
    CallOption c;
    c.setStrike(K);
    c.setMaturity(T);
    
    MonteCarloPricer pricer;
    pricer.nSteps = 100;
    std::cout << "Monte Carlo pricer Heston: " << pricer.price(c, hest_euler) << std::endl;
    std::cout << "Monte Carlo pricer Scott: " << pricer.price(c, scott_euler) << std::endl;
}

static void testDeltaCall(){
    CallOption c;
    c.setStrike( 110 );
    c.setMaturity( 2 );

    BlackScholesModel m;
    m.setRiskFreeRate(0.05);
    m.setVolatility(0.1);
    m.setStockPrice(100.0);
    m.setDate(1.0);
    m.setDrift(0.1);
    /*
    m.volatility = 0.1;
    m.riskFreeRate = 0.05;
    m.stockPrice = 100.0;
    m.drift = 0.1;
    m.date = 1;
    */
    
    MonteCarloPricer pricer;
    double answer = pricer.delta(c, m);
    double expected = c.delta(m);
    
    ASSERT_APPROX_EQUAL(answer, expected, 0.1);
}

static void testDeltaPut(){
    PutOption p;
    p.setStrike( 90 );
    p.setMaturity( 2 );

    BlackScholesModel m;
    m.setRiskFreeRate(0.05);
    m.setVolatility(0.1);
    m.setStockPrice(100.0);
    m.setDate(1.0);
    m.setDrift(0.1);
    
    MonteCarloPricer pricer;
    double answer = pricer.delta(p, m);
    double expected = p.delta(m);
    
    ASSERT_APPROX_EQUAL(answer, expected, 0.1);
}

void testMonteCarloPricer() {
    TEST( testPriceCallOption );
    TEST( testDeltaCall );
    TEST( testDeltaPut );
    TEST( testHestonCallOption );
    TEST( testHestonVsScott );
}
