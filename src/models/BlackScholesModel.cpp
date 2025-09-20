#include "PricingLibrary/models/BlackScholesModel.h"

using namespace std;

#include "PricingLibrary/utils/matlib.h"

/*
BlackScholesModel::BlackScholesModel() :
    drift(0.0),
    stockPrice(0.0),
    volatility(0.0),
    riskFreeRate(0.0),
    date(0.0) {
}
 */

/**
 *  Creates a price path according to the model parameters
 */
vector<double> BlackScholesModel::generateRiskNeutralPricePath(double toDate, int nSteps ) const {
    return generatePricePath( toDate, nSteps, getRiskFreeRate() );
}

/**
 *  Creates a price path according to the model parameters
 */
vector<double> BlackScholesModel::generatePricePath(double toDate, int nSteps ) const {
    return generatePricePath( toDate, nSteps, getDrift() );
}


/**
 *  Creates a price path according to the model parameters
 */
vector<double> BlackScholesModel::generatePricePath(double toDate, int nSteps, double drift ) const {
    vector<double> path(nSteps,0.0);
    vector<double> epsilon = randn( nSteps );
    double dt = (toDate-getDate())/nSteps;
    double a = (drift - getVolatility()*getVolatility()*0.5)*dt;
    double b = getVolatility()*sqrt(dt);
    double currentLogS = log( getStockPrice() );
    for (int i=0; i<nSteps; i++) {
        double dLogS = a + b*epsilon[i];
        double logS = currentLogS + dLogS;
        path[i] = exp( logS );
        currentLogS = logS;
    }
    return path;
}

std::vector<double> BlackScholesModel::generateRiskNeutralDeltaPath(double toDate, int nSteps, double price, vector<double> epsilon) const{
    return generateDeltaPath(toDate, nSteps, getRiskFreeRate(), price, epsilon);
}

std::vector<double> BlackScholesModel::generateDeltaPath(double toDate, int nSteps, double drift, double price, vector<double> epsilon) const{
    vector<double> path(nSteps,0.0);
    
    double dt = (toDate-getDate())/nSteps;
    double a = (drift - 0.5*getVolatility()*getVolatility())*dt;
    double b = getVolatility()*sqrt(dt);
    double s = log(price);
    
    for(int i = 0; i<nSteps; i++){
    s = s + a + b*epsilon[i];
    path[i] = exp(s);
    }
    return path;
}

////////////////////////////////
//
//   TESTS
//
////////////////////////////////

void testRiskNeutralPricePath() {
    rng("default");

    BlackScholesModel bsm;
    bsm.setRiskFreeRate(0.05);
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100.0);
    bsm.setDate(2.0);

    int nPaths = 10000;
    int nsteps = 5;
    double maturity = 4.0;
    vector<double> finalPrices(nPaths,0.0);
    for (int i=0; i<nPaths; i++) {
        vector<double> path = bsm.generateRiskNeutralPricePath( maturity, nsteps );
        finalPrices[i] = path.back();
    }
    ASSERT_APPROX_EQUAL( mean( finalPrices ), exp( bsm.getRiskFreeRate()*2.0)*bsm.getStockPrice(), 0.5);
}

void testVisually() {
    BlackScholesModel bsm;
    bsm.setRiskFreeRate(0.05);
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100.0);
    bsm.setDate(2.0);

    int nSteps = 1000;
    double maturity = 4.0;

    vector<double> path = bsm.generatePricePath( maturity,
                                                 nSteps );
    double dt = (maturity-bsm.getDate())/nSteps;
    vector<double> times = linspace(dt, maturity, nSteps );
    plot("examplePricePath.html",times, path );
}


void testBlackScholesModel() {
    TEST( testRiskNeutralPricePath );
    TEST( testVisually );
}

