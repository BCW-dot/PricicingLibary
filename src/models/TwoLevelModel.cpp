#include "PricingLibrary/models/TwoLevelModel.h"
#include "PricingLibrary/utils/matlib.h"
#include "PricingLibrary/models/BlackScholesModel.h"

using namespace std;

TwoLevelModel::TwoLevelModel():
    volatility2(0.1)
{}
 
std::vector<double> TwoLevelModel::generatePricePath(double toDate, int nSteps) const {
    vector<double> path(nSteps,0.0);
    vector<double> epsilon = randn( nSteps );
    double dt = (toDate-getDate())/nSteps;
    double a = (getDrift() - getVolatility()*getVolatility()*0.5)*dt;
    double b = getVolatility()*sqrt(dt);
    double currentLogS = log( getStockPrice() );
    
    int n_halfe = int(nSteps/2);
    for (int i=0; i<n_halfe; i++) {
        double dLogS = a + b*epsilon[i];
        double logS = currentLogS + dLogS;
        path[i] = exp( logS );
        currentLogS = logS;
    }
    
    a = (getDrift() - volatility2*volatility2*0.5)*dt;
    b = volatility2*sqrt(dt);
    for (int i=n_halfe; i<nSteps; i++) {
        double dLogS = a + b*epsilon[i];
        double logS = currentLogS + dLogS;
        path[i] = exp( logS );
        currentLogS = logS;
    }
    
    return path;
}

static void testTwoVolatilities(){
    rng("default");
    TwoLevelModel model;
    
    model.setRiskFreeRate(0.05);
    model.setStockPrice(100.0);
    model.setDate(2.0);
    model.setVolatility(0.1);
    model.volatility2 = model.getVolatility();
    
    model.setRiskFreeRate(0.3);
    model.setDrift(model.getRiskFreeRate());

    int nPaths = 10000;
    int nsteps = 5;
    double maturity = 4.0;
    vector<double> finalPrices(nPaths,0.0);
    for (int i=0; i<nPaths; i++) {
        vector<double> path = model.generatePricePath( maturity, nsteps );
        finalPrices[i] = path.back();
    }
    ASSERT_APPROX_EQUAL( mean( finalPrices ), exp( model.getRiskFreeRate()*2.0)*model.getStockPrice(), 0.5);
}

static void testVisually(){
    TwoLevelModel model;
    
    model.setRiskFreeRate(0.05);
    model.setVolatility(0.1);
    model.setStockPrice(100.0);
    model.setDate(1.0);
    model.volatility2 = 0.8;

    int nSteps = 1000;
    double maturity = 4.0;

    vector<double> path = model.generatePricePath( maturity,
                                                 nSteps );
    double dt = (maturity-model.getDate())/nSteps;
    vector<double> times = linspace(dt, maturity, nSteps );
    plot("examplePricePathTwoVolatilities.html", times, path );
}

void testTwoLevelModel(){
    TEST( testTwoVolatilities );
    TEST( testVisually );
}

