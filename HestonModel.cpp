#include "HestonModel.h"
#include "CallOption.h"
#include "matlib.h"
#include "LineChart.h"

using namespace std;

std::vector<double> HestonModel::generatePricePath(double toDate, int nSteps) const{
    std::vector<double> vol_draws(nSteps, 0.0);
    std::vector<double> spot_draws(nSteps, 0.0);
    fill_correlated_normals(vol_draws, spot_draws, rho);
    
    std::vector<double> vol_path(nSteps, 0.0);
    vol_path = calc_vol_path(toDate, nSteps, vol_draws);
    return generatePricePath(toDate, nSteps, vol_path, spot_draws);
}

std::vector<double> HestonModel::generateRiskNeutralPricePath(double toDate, int nSteps) const{
    std::vector<double> vol_draws(nSteps, 0.0);
    std::vector<double> spot_draws(nSteps, 0.0);
    fill_correlated_normals(vol_draws, spot_draws, rho);
    
    std::vector<double> vol_path(nSteps, 0.0);
    vol_path = calc_vol_path(toDate, nSteps, vol_draws);
    return generateRiskNeutralPricePath(toDate, nSteps, vol_path, spot_draws);
}

vector<double> HestonModel::generatePricePath(double toDate, int nSteps, std::vector<double>& vol_path, std::vector<double>& spot_draws) const{
    vector<double> spot_path(nSteps,0.0);
    spot_path[0] = getStockPrice();
    double dt = (toDate-getDate())/nSteps;
    double sqrt_dt = sqrt(dt);

    // Create the spot price path making use of the volatility
    // path. Uses a similar Euler Truncation method to the vol path.
    for (int i=1; i<nSteps; i++) {
        double v_max = std::max(vol_path[i-1], 0.0);
        spot_path[i] = spot_path[i-1] + getDrift() * spot_path[i-1] * dt + sqrt(v_max) *sqrt_dt * spot_path[i-1] * spot_draws[i-1];
    }
    
    return spot_path;
}

/* This need correction in the Q-measure term, sqrt_dt placement is unsure right now
  * should be working now, there was a mistake in the equation of x[i+1] and in Q_measure_brownian
 * chnaged from: //double Q_measure_brownian = sqrt_dt * spot_draws[i-1] + ( (getDrift() - getRiskFreeRate())/sqrt_v_max ) * dt;
                 //spot_path[i] = spot_path[i-1] + getRiskFreeRate() * spot_path[i-1] * dt + sqrt_v_max * spot_path[i-1] * Q_measure_brownian;
 * to current code, i.e. a - in Q_bm and getRiskfreerate() to getDrift()
 */
vector<double> HestonModel::generateRiskNeutralPricePath(double toDate, int nSteps, std::vector<double>& vol_path, std::vector<double>& spot_draws) const{
    vector<double> spot_path(nSteps,0.0);
    spot_path[0] = getStockPrice();
    double dt = (toDate-getDate())/nSteps;
    double sqrt_dt = sqrt(dt);
    double epsilon = 1e-5; //need this since we would devite by zero else in Q_measure_BM

    for (int i=1; i<nSteps; i++) {
        double sqrt_v_max = sqrt(std::max(vol_path[i-1], epsilon));
        //double Q_measure_brownian = sqrt_dt * spot_draws[i-1] + ( (getDrift() - getRiskFreeRate())/sqrt_v_max ) * dt;
        //changing + to - from page 207
        double Q_measure_brownian = sqrt_dt * spot_draws[i-1] - ( (getDrift() - getRiskFreeRate())/sqrt_v_max ) * dt;
        spot_path[i] = spot_path[i-1] + getDrift() * spot_path[i-1] * dt + sqrt_v_max * spot_path[i-1] * Q_measure_brownian;

        //spot_path[i] = spot_path[i-1] + getRiskFreeRate() * spot_path[i-1] * dt + sqrt_v_max * spot_path[i-1] * Q_measure_brownian;
    }
    return spot_path;
}

vector<double> HestonModel::calc_vol_path(double toDate, int nSteps, const std::vector<double>& vol_draws) const{
    vector<double> vol_path(nSteps,0.0);
    
    //size_t vec_size = vol_draws.size();
    double dt = (toDate-getDate())/nSteps;

      // Iterate through the correlated random draws vector and
      // use the 'Full Truncation' scheme to create the volatility path
      for (int i=1; i<nSteps; i++) {
        double v_max = std::max(vol_path[i-1], 0.0);
        vol_path[i] = vol_path[i-1] + kappa * dt * (theta - v_max) + xi * sqrt(v_max*dt) * vol_draws[i-1];
      }
    
    return vol_path;
}

/* this method was written in order to test the risk neutral pricing. Here we first do a Cholesky decomposition to transfrom the SDE with correlates
 Brownain motions to an SDE with uncorrelated Brownians. This means we can generate normals while in the loop and not before hand. Secondly on this transformed SDE a measure chance is applied in order for the price process to be a local martinglae. */
vector<double> HestonModel::independent_price_path(double toDate, int nSteps) const{
    vector<double> spot_path(nSteps,0.0);
    spot_path[0] = getStockPrice();
    double v = getVolatility();
    
    double dt = (toDate-getDate())/nSteps;
    double sqrt_dt = sqrt(dt);
    
    vector<double> b1 = randn( nSteps );
    vector<double> b2 = randn( nSteps );
    
    double sqrt_rho = sqrt(1-rho*rho);

    // Create the spot price path making use of the volatility
    // path. Uses a similar Euler Truncation method to the vol path.
    for (int i=1; i<nSteps; i++) {
        //vol "path"
        v = std::max(v + kappa * dt * (theta - v) + xi * sqrt(v) * (rho * sqrt_dt * b1[i] + sqrt_dt * sqrt_rho * b2[i]), 0.0);
        //price path
        spot_path[i] = spot_path[i-1] + getDrift() * spot_path[i-1] * dt + sqrt(v) * spot_path[i-1] * sqrt_dt * b1[i];
    }
    
    return spot_path;
}

vector<double> HestonModel::Q_measure_price_path(double toDate, int nSteps) const{
    vector<double> spot_path(nSteps,0.0);
    spot_path[0] = getStockPrice();
    double v = getVolatility();
    
    double dt = (toDate-getDate())/nSteps;
    double sqrt_dt = sqrt(dt);
    
    vector<double> b1 = randn( nSteps );
    vector<double> b2 = randn( nSteps );
    
    double sqrt_rho = sqrt(1-rho*rho);

    // Create the spot price path making use of the volatility
    // path. Uses a similar Euler Truncation method to the vol path.
    for (int i=1; i<nSteps; i++) {
        //vol "path"
        v = std::max(v + kappa * dt * (theta - v) + xi * sqrt(v) * (rho * sqrt_dt * b1[i] + sqrt_dt * sqrt_rho * b2[i]), 0.001);
        //Q measure
        double W_Q = sqrt_dt * b1[i] + (getDrift() - getRiskFreeRate())/sqrt(v) * dt;
        //price path
        spot_path[i] = spot_path[i-1] + getRiskFreeRate() * spot_path[i-1] * dt + sqrt(v) * spot_path[i-1] * W_Q;
    }
    
    return spot_path;
}

/*
 *
 *
 * TEST
 *
 *
 */

static void testHestonVisually(){
    my_rng();
    // First we create the parameter list
    // Note that you could easily modify this code to input the parameters
    // either from the command line or via a file
    unsigned num_intervals = 1000;   // Number of simulated asset paths

    double S_0 = 100.0;    // Initial spot price
    double r = 0.03;     // Risk-free rate
    double v_0 = 0.01; // Initial volatility
    double T = 1.00;       // One year until expiry

    double rho = 0.8;     // Correlation of asset and volatility
    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.019;  // Long run average volatility
    double xi = 0.61;      // "Vol of vol"

    HestonModel hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(0.0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;

    // Create the spot and vol initial normal and price paths
    std::vector<double> spot_draws(num_intervals, 0.0);  // Vector of initial spot normal draws
    std::vector<double> vol_draws(num_intervals, 0.0);   // Vector of initial correlated vol normal draws
    std::vector<double> spot_prices(num_intervals, 0);  // Vector of spot path
    std::vector<double> vol_prices(num_intervals, 0);   // Vector of vol path

    //generate_normal_correlation_paths(rho, spot_draws, vol_draws);
    fill_correlated_normals(vol_draws, spot_draws, rho);
    vol_prices = hest_euler.getVolPath(T, num_intervals, vol_draws);// hest_euler.calc_vol_path(T, num_intervals, vol_draws);
    spot_prices = hest_euler.getPath(T, num_intervals, vol_prices, spot_draws);
    
    double dt = (T-hest_euler.getDate())/num_intervals;
    vector<double> times = linspace(dt, T, num_intervals );
    
    //std::string filename = "/Users/benewilkens/Documents/ArmstrongFin/CSV/Heston.csv";
    //writeVectorsToCsv(spot_prices, times, vol_prices, filename);

    plot("exampleHestonPricePath.html",times, spot_prices );
    plot("exampleHestonVolPath.html",times, vol_prices );
}

static void testCallPrice(){
    my_rng();
    // First we create the parameter list
    // Note that you could easily modify this code to input the parameters
    // either from the command line or via a file
    unsigned num_sims = 1000;   // Number of simulated asset paths
    unsigned num_intervals = 1000;  // Number of intervals for the asset path to be sampled

    double S_0 = 100.0;    // Initial spot price
    double K = 100.0;      // Strike price
    double r = 0.03;     // Risk-free rate
    //double drift = 0.04;
    double v_0 = 0.01; // Initial volatility
    double T = 1.00;       // One year until expiry
    double date_0 = 0.0;

    double rho = 0.0;     // Correlation of asset and volatility
    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.019;  // Long run average volatility
    double xi = 0.61;      // "Vol of vol"

    // Create the PayOff, Option and Heston objects
    CallOption c;
    c.setStrike(K);
    c.setMaturity(T);

    HestonModel hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setDrift(hest_euler.getRiskFreeRate());
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(date_0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;

    // Create the spot and vol initial normal and price paths
    std::vector<double> spot_draws(num_intervals, 0.0);  // Vector of initial spot normal draws
    std::vector<double> vol_draws(num_intervals, 0.0);   // Vector of initial correlated vol normal draws
    std::vector<double> spot_prices(num_intervals, S_0);  // Vector of initial spot prices
    std::vector<double> vol_prices(num_intervals, pow(v_0,2));   // Vector of initial vol prices

    // Monte Carlo options pricing
    double payoff_sum = 0.0;
    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (unsigned i=0; i<num_sims; i++) {
      fill_correlated_normals(vol_draws, spot_draws, rho);
      //vol_prices = hest_euler.calc_vol_path(T, num_intervals, vol_draws);
      spot_prices = hest_euler.getPath(T, num_intervals, vol_prices, spot_draws);
      payoff_sum += c.payoff(spot_prices);
    }
    /*
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << " Option was priced in "<< pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
    double option_price = (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
    double option_price1 = (payoff_sum1 / static_cast<double>(num_sims)) * exp(-r*T);
    std::cout << "Option Price Heston: " << option_price << std::endl;
    std::cout << "Option Price Heston 1: " << option_price1 << std::endl;

    BlackScholesModel bsm;
    bsm.setRiskFreeRate(r);
    bsm.setVolatility(v_0);
    bsm.setStockPrice(S_0);
    bsm.setDate(date_0);
    
    std::cout << "Option Price Black Scholes: " << c.price(bsm) << std::endl;
    */
    double option_price = (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
    
    BlackScholesModel bsm;
    bsm.setRiskFreeRate(r);
    bsm.setVolatility(v_0);
    bsm.setStockPrice(S_0);
    bsm.setDate(date_0);
    
    ASSERT_APPROX_EQUAL(option_price, c.price(bsm), 0.1);
}

static void testHestonPayoff(){
    my_rng();
    // First we create the parameter list
    // Note that you could easily modify this code to input the parameters
    // either from the command line or via a file
    unsigned num_sims = 1000;   // Number of simulated asset paths
    unsigned num_intervals = 1000;  // Number of intervals for the asset path to be sampled

    double S_0 = 100.0;    // Initial spot price
    double K = 100.0;      // Strike price
    double r = 0.03;     // Risk-free rate
    double v_0 = 0.01; // Initial volatility
    double date_0 = 0.0;    //Initial date
    double T = 1.00;       // One year until expiry

    double rho = 0.8;     // Correlation of asset and volatility
    double kappa = 3.0;   // Mean-reversion rate
    double theta = 0.02;  // Long run average volatility
    double xi = 0.01;      // "Vol of vol"

    // Create the PayOff, Option and Heston objects
    CallOption c;
    c.setStrike(K);
    c.setMaturity(T);

    HestonModel hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(date_0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;

    std::vector<double> spot_prices;//(num_intervals, S_0);  // Vector of initial spot prices
    std::vector<double> hist_of_last_prices(num_intervals, 0.0);
    
    // Monte Carlo options pricing
    double payoff_sum = 0.0;
    for (int i=0; i<num_sims; i++) {
      spot_prices = hest_euler.generateRiskNeutralPricePath(T, num_intervals);//hest_euler.generatePricePath(T, num_intervals, vol_prices, spot_draws);
      payoff_sum += c.payoff(spot_prices);
      hist_of_last_prices[i] = spot_prices.back();
    }
    double option_price = (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
    std::cout << "Option Price Heston: " << option_price << std::endl;
    //hist("HestonPrices.html",hist_of_last_prices,100);
    //open_hist("HestonPrices.html");
    
    BlackScholesModel bsm;
    bsm.setRiskFreeRate(r);
    bsm.setVolatility(v_0);
    bsm.setStockPrice(S_0);
    bsm.setDate(date_0);
    
    std::cout << "Black Shcoles price: " << c.price(bsm) << std::endl;
}

static void testGenPricePath(){
    unsigned num_intervals = 1000;   // Number of simulated asset paths

    double S_0 = 100.0;    // Initial spot price
    double r = 0.03;     // Risk-free rate
    double v_0 = 0.01; // Initial volatility
    double T = 1.00;       // One year until expiry

    double rho = 0.8;     // Correlation of asset and volatility
    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.019;  // Long run average volatility
    double xi = 0.61;      // "Vol of vol"

    HestonModel hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(0.0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;

    std::vector<double> spot_prices(num_intervals,0.0);
    std::vector<double> spot_prices_RiskNeutral(num_intervals,0.0);
    
    spot_prices = hest_euler.generatePricePath(T, num_intervals);
    spot_prices_RiskNeutral = hest_euler.generateRiskNeutralPricePath(T, num_intervals);
    
    double dt = (T-hest_euler.getDate())/num_intervals;
    vector<double> times = linspace(dt, T, num_intervals );

    plot("exampleHestonGenPath.html",times, spot_prices );
    //open_plot("exampleHestonGenPath.html");
    
    plot("exampleHestonGenPathRiskNeutral.html",times, spot_prices_RiskNeutral );
    //open_plot("exampleHestonGenPathRiskNeutral.html");
}

//Here is a mistake, but it got "fixed" in RiskNeutralPricingIndependet(), how to fix it here:
/*
 1. check formulas for spot price adn vol path again
 2. check measure change
 3. check correlated vectors are correctly produced
 4. give up
 */
static void testRiskNeutralPricing(){
    /* This tests Risk Neutral pricing of the heston model */

    /* DOESNT WORK YET, but should work now, still doesnt tho*/
    my_rng();
    unsigned num_sims = 1000;   // Number of simulated asset paths
    unsigned num_intervals = 1000;  // Number of intervals for the asset path to be sampled

    double S_0 = 100.0;    // Initial spot price
    double r = 0.03;     // Risk-free rate
    double v_0 = 0.09; // Initial volatility
    double date_0 = 0.0;    //Initial date
    double T = 1.00;       // One year until expiry

    double rho = 0.8;     // Correlation of asset and volatility
    double kappa = 3.0;   // Mean-reversion rate
    double theta = 0.02;  // Long run average volatility
    double xi = 0.09;      // "Vol of vol"

    HestonModel hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(date_0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;

    std::vector<double> spot_prices;//(num_intervals, S_0);  // Vector of initial spot prices
    std::vector<double> hist_of_last_prices(num_intervals, 0.0);
    
    for (int i=0; i<num_sims; i++) {
      spot_prices = hest_euler.generateRiskNeutralPricePath(T, num_intervals);
      hist_of_last_prices[i] = spot_prices.back(); //exp(-hest_euler.getRiskFreeRate()*T)*
    }
    std::cout << mean(hist_of_last_prices) << std::endl;
    std::cout << exp(hest_euler.getRiskFreeRate()*T)*hest_euler.getStockPrice() << std::endl;
    hist("HestonRiskNeutralPrices.html",hist_of_last_prices,100);
    //open_hist("HestonRiskNeutralPrices.html");
}

static void testIndependentPath(){
    my_rng();
    unsigned num_intervals = 1000;   // Number of simulated asset paths

    double S_0 = 100.0;    // Initial spot price
    double r = 0.03;     // Risk-free rate
    double v_0 = 0.01; // Initial volatility
    double T = 1.00;       // One year until expiry

    double rho = 0.8;     // Correlation of asset and volatility
    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.019;  // Long run average volatility
    double xi = 0.61;      // "Vol of vol"

    HestonModel hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(0.0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;
    
    vector<double> spot_prices(num_intervals,0.0);
    vector<double> spot_prices_underQ(num_intervals,0.0);

    spot_prices = hest_euler.independent_price_path(T, num_intervals);
    //my_rng();
    spot_prices_underQ = hest_euler.Q_measure_price_path(T, num_intervals);
    
    double dt = (T-hest_euler.getDate())/num_intervals;
    vector<double> times = linspace(dt, T, num_intervals );
    plot("IndependentHestonPricePath.html",times, spot_prices );
    plot("QMeasureHestonPricePath.html",times, spot_prices_underQ );
}

static void testRiskNeutralQPricingIndpendent(){
    /* This tests Risk Neutral pricing of the heston model */
    my_rng();
    unsigned num_sims = 1000;   // Number of simulated asset paths
    unsigned num_intervals = 1000;  // Number of intervals for the asset path to be sampled

    double S_0 = 100.0;    // Initial spot price
    double r = 0.03;     // Risk-free rate
    double v_0 = 0.09; // Initial volatility
    double date_0 = 0.0;    //Initial date
    double T = 1.00;       // One year until expiry

    double rho = 0.8;     // Correlation of asset and volatility
    double kappa = 3.0;   // Mean-reversion rate
    double theta = 0.02;  // Long run average volatility
    double xi = 0.09;      // "Vol of vol"

    HestonModel hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(date_0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;

    std::vector<double> spot_prices;//(num_intervals, S_0);  // Vector of initial spot prices
    std::vector<double> hist_of_last_prices(num_intervals, 0.0);
    
    for (int i=0; i<num_sims; i++) {
      spot_prices = hest_euler.Q_measure_price_path(T, num_intervals);
      hist_of_last_prices[i] = spot_prices.back(); //exp(-hest_euler.getRiskFreeRate()*T)*
    }
    std::cout << mean(hist_of_last_prices) << std::endl;
    std::cout << exp(hest_euler.getRiskFreeRate()*T)*hest_euler.getStockPrice() << std::endl;
    hist("HestonRiskNeutralindependentPrices.html",hist_of_last_prices,100);
    //open_hist("HestonRiskNeutralindependentPrices.html");
}

void testHestonModel(){
    TEST( testHestonPayoff );
    TEST( testHestonVisually );
    TEST( testCallPrice );
    TEST( testGenPricePath );
    TEST( testRiskNeutralPricing );
    TEST( testIndependentPath );
    TEST( testRiskNeutralQPricingIndpendent );
}
