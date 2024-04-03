#include "MonteCarloPricer.h"
#include "BlackScholesModel.h"

#include "matlib.h"
#include "CallOption.h"
#include "PutOption.h"

using namespace std;

MonteCarloPricer::MonteCarloPricer() :
    nScenarios(10000),
    nSteps(100) {
}

double MonteCarloPricer::price(const ContinuousTimeOption& option, const BlackScholesModel& model ) {
    int nSteps = this->nSteps;
    if (!option.isPathDependent()) {
        nSteps = 1;
    }
    double total = 0.0;
    for (int i=0; i<nScenarios; i++) {
        vector<double> path= model.generateRiskNeutralPricePath(option.getMaturity(), nSteps );
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
    rng("default");

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
    double price = pricer.price( c, m );
    double expected = c.price( m );
    ASSERT_APPROX_EQUAL( price, expected, 0.1 );
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
    /*
    m.volatility = 0.1;
    m.riskFreeRate = 0.05;
    m.stockPrice = 100.0;
    m.drift = 0.1;
    m.date = 1;
     */
    
    MonteCarloPricer pricer;
    double answer = pricer.delta(p, m);
    double expected = p.delta(m);
    
    ASSERT_APPROX_EQUAL(answer, expected, 0.1);
}

void testMonteCarloPricer() {
    TEST( testPriceCallOption );
    TEST( testDeltaCall );
    TEST( testDeltaPut );
}
