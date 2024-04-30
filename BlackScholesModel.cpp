#include "BlackScholesModel.h"

using namespace std;

#include "matlib.h"

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
/*
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
 */

vector<double> BlackScholesModel::generatePricePath(double toDate, int nSteps, double drift ) const {
    vector<double> path(nSteps,0.0);
    path[0] = getStockPrice();
    vector<double> epsilon = randn( nSteps );
    double dt = (toDate-getDate())/nSteps;
    double sqrt_dt = sqrt(dt);
    for (int i=1; i<nSteps; i++) {
        path[i] = path[i-1] + drift * path[i-1] * dt + getVolatility() * path[i-1] * sqrt_dt * epsilon[i];
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
    //rng("default");
    my_rng();

    BlackScholesModel bsm;
    bsm.setRiskFreeRate(0.05);
    bsm.setVolatility(0.1);
    bsm.setStockPrice(10.0);
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
    //std::cout << mean( finalPrices ) << std::endl;
    //std::cout << exp( bsm.getRiskFreeRate()*2.0)*bsm.getStockPrice() << std::endl;
}

void testVisually() {
    my_rng();
    
    BlackScholesModel bsm;
    bsm.setRiskFreeRate(0.05);
    bsm.setVolatility(0.1);
    bsm.setDrift(0.05);
    bsm.setStockPrice(100.0);
    bsm.setDate(0.0);

    int nSteps = 1000;
    double maturity = 1.0;

    vector<double> path = bsm.generateRiskNeutralPricePath( maturity, nSteps );
    
    double dt = (maturity-bsm.getDate())/nSteps;
    vector<double> times = linspace(dt, maturity, nSteps );
    plot("examplePricePath.html",times, path );
    //open_plot("examplePricePath.html");
}

void testRiskNeutralPricing(){
    /* This test checks that GBM under the risk neutral measure is a martinagle and normally distributed */
    my_rng();
    // First we create the parameter list
    // Note that you could easily modify this code to input the parameters
    // either from the command line or via a file
    unsigned num_sims = 1000;   // Number of simulated asset paths
    unsigned num_intervals = 1000;  // Number of intervals for the asset path to be sampled

    double S_0 = 100.0;    // Initial spot price
    double r = 0.03;     // Risk-free rate
    double v_0 = 0.01; // Initial volatility
    double date_0 = 0.0;    //Initial date
    double T = 1.00;       // One year until expiry

    BlackScholesModel bsm;
    bsm.setRiskFreeRate(r);
    bsm.setVolatility(v_0);
    bsm.setStockPrice(S_0);
    bsm.setDate(date_0);

    std::vector<double> spot_prices;//(num_intervals, S_0);  // Vector of initial spot prices
    std::vector<double> hist_of_last_prices(num_intervals, 0.0);
    
    for (int i=0; i<num_sims; i++) {
      spot_prices = bsm.generateRiskNeutralPricePath(T, num_intervals);
      hist_of_last_prices[i] = exp(-bsm.getRiskFreeRate()*T)*spot_prices.back();
    }
    hist("GBMRiskNeutralPrices.html",hist_of_last_prices,100);
    open_hist("GBMRiskNeutralPrices.html");
}

/*
void testRandomNumberGeneration(){
     int n = 10000;
     int M = 3;
     n = M;
     //my_rng();
     //vector<double> nor = my_randn(n);
     
     for(int i =0; i<M; i++){
         vector<double> nor = randn(n);
         print_window(nor, 0, 10);
         std::cout << "Mean: " << mean(nor) << ", " << std::endl;
     }
     //hist("NormalTestMYrandn.html", nor, 100);
     //open_hist("NormalTestMYrandn.html");
 }

void testRandomNumberGeneration1(){
    int n = 10000;
    vector<double> nor = randn(n);
    print_window(nor, 0, 10);
    std::cout << mean(nor) << std::endl;
    hist("NormalTestRandn.html", nor, 100);
    //open_hist("NormalTestRandn.html");
}

static void testGBM(){
    BlackScholesModel bsm;
    bsm.setRiskFreeRate(0.05);
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100.0);
    bsm.setDate(0.0);

    int nSteps = 1000;
    double maturity = 1.0;

    vector<double> path = bsm.generatePricePath( maturity, nSteps );
    vector<double> path1 = bsm.StockPriceModel::generateRiskNeutralPricePath(maturity, nSteps);
    
    double dt = (maturity-bsm.getDate())/nSteps;
    vector<double> times = linspace(dt, maturity, nSteps );
    
    plot("examplePricePathmy.html",times, path );
    plot("examplePricePathmy1.html",times, path1 );
    //open_plot("examplePricePathmy.html");
    //open_plot("examplePricePathmy1.html");
}
 
 void testRandomNumberGeneration(){
     int n = 10000;
     int M = 3;
     n = M;
     //my_rng();
     //vector<double> nor = my_randn(n);
     for(int i =0; i<M; i++){
         vector<double> nor = randn(n);
         print_window(nor, 0, 10);
         std::cout << "Mean: " << mean(nor) << ", " << std::endl;
     }
    
     //hist("NormalTestMYrandn.html", nor, 100);
     //open_hist("NormalTestMYrandn.html");
 }
*/

void testBlackScholesModel() {
    TEST( testRiskNeutralPricePath );
    TEST( testVisually );
    TEST( testRiskNeutralPricing );
    /*
    TEST( testRandomNumberGeneration );
    TEST( testRandomNumberGeneration1 );
    TEST( testGBM );
     */
}

