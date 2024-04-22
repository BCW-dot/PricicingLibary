#include "ScottChesney.h"
#include "CallOption.h"
#include "matlib.h"
#include "LineChart.h"
//#include "HestonModel.h"

using namespace std;

std::vector<double> ScottChesney::generatePricePath(double toDate, int nSteps) const{
    std::vector<double> vol_draws(nSteps, 0.0);
    std::vector<double> spot_draws(nSteps, 0.0);
    fill_correlated_normals(vol_draws, spot_draws, rho);
    
    std::vector<double> vol_path(nSteps, 0.0);
    vol_path = calc_vol_path(toDate, nSteps, vol_draws);

    return generatePricePath(toDate, nSteps, vol_path, spot_draws);
}

/*
std::vector<double> ScottChesney::generateRiskNeutralPricePath(double toDate, int nSteps) const{
    std::vector<double> vol_draws(nSteps, 0.0);
    std::vector<double> spot_draws(nSteps, 0.0);
    fill_correlated_normals(vol_draws, spot_draws, rho);
    
    std::vector<double> vol_path(nSteps, 0.0);
    vol_path = calc_vol_path(toDate, nSteps, vol_draws);
    //this doesnt work! need to add correlate_(vol,spot)
    return generatePricePath(toDate, nSteps, vol_path, spot_draws);
}
*/

vector<double> ScottChesney::generatePricePath(double toDate, int nSteps, std::vector<double>& vol_path, std::vector<double>& spot_draws) const{
    vector<double> spot_path(nSteps,0.0);
    spot_path[0] = getStockPrice();
    double dt = (toDate-getDate())/nSteps;
    double sqrt_dt = sqrt(dt);

    // Create the spot price path making use of the volatility
    for (int i=1; i<nSteps; i++) {
      double v_e = std::exp(vol_path[i-1]);
      //std::cout<< "V_max: " << v_max << ", ";
      //spot_path[i] = spot_path[i-1] * exp( (getRiskFreeRate() - 0.5*v_e)*dt + v_e * sqrt(dt)*spot_draws[i-1]);
      spot_path[i] = spot_path[i-1] + getDrift() * spot_path[i-1] * dt + v_e * spot_path[i-1] * sqrt_dt * spot_draws[i-1];
    }
    
    return spot_path;
}

vector<double> ScottChesney::calc_vol_path(double toDate, int nSteps, const std::vector<double> &vol_draws) const{
    vector<double> vol_path(nSteps,0.0);
    double dt = (toDate-getDate())/nSteps;
    
    for (int i=1; i<nSteps; i++) {
      //double v_max = std::max(vol_path[i-1], 0.0);
      //vol_path[i] = vol_path[i-1] + kappa * dt * (theta - v_max) + xi * sqrt(dt) * vol_draws[i-1];
      vol_path[i] = std::max(vol_path[i-1] + kappa * dt * (theta - vol_path[i-1]) + xi * sqrt(dt) * vol_draws[i-1], 0.0);
    }
    
    return vol_path;
}

static void testScottChesneyVisually(){
    my_rng();
    // First we create the parameter list
    unsigned num_intervals = 1000;  // Number of intervals for the asset path to be sampled

    double S_0 = 100.0;    // Initial spot price
    //double K = 100.0;      // Strike price
    double r = 0.03;     // Risk-free rate
    double drift = 0.04;
    double v_0 = 0.01; // Initial volatility
    double T = 1.00;       // One year until expiry
    double date_0 = 0.0;

    double rho = 0.9;     // Correlation of asset and volatility
    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.02;  // Long run average volatility
    double xi = 0.02;      // "Vol of vol"

    ScottChesney hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setDrift(drift);
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
    std::vector<double> vol_prices(num_intervals, v_0);   // Vector of initial vol prices

    fill_correlated_normals(vol_draws, spot_draws, rho);

    vol_prices = hest_euler.calc_vol_path(T, num_intervals, vol_draws);
    spot_prices = hest_euler.generatePricePath(T, num_intervals, vol_prices, spot_draws);
    
    double dt = (T-hest_euler.getDate())/num_intervals;
    vector<double> times = linspace(dt, T, num_intervals );
    
    std::string filename = "/Users/benewilkens/Documents/ArmstrongFin/CSV/ScottChesney.csv";
    writeVectorsToCsv(spot_prices, times, vol_prices, filename);

    plot("exampleScottChesneyPricePath.html",times, spot_prices );
    plot("exampleScottChesneyVolPath.html",times, vol_prices );
    //open_plot("exampleScottChesneyPricePath.html");
    //open_plot("exampleScottChesneyVolPath.html");
}

static void testScottChesneyPayoff(){
    my_rng();
    // First we create the parameter list
    // Note that you could easily modify this code to input the parameters
    // either from the command line or via a file
    unsigned num_sims = 1000;   // Number of simulated asset paths
    unsigned num_intervals = 1000;  // Number of intervals for the asset path to be sampled

    double S_0 = 100.0;    // Initial spot price
    double K = 100.0;      // Strike price
    double r = 0.03;     // Risk-free rate
    double drift = 0.04;
    double v_0 = 0.010201; // Initial volatility
    double T = 1.00;       // One year until expiry
    double date_0 = 0.0;

    double rho = 0.9;     // Correlation of asset and volatility
    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.02;  // Long run average volatility
    double xi = 0.1;      // "Vol of vol"

    ScottChesney hest_euler;
    hest_euler.setRiskFreeRate(r);
    hest_euler.setDrift(drift);
    hest_euler.setVolatility(v_0);
    hest_euler.setStockPrice(S_0);
    hest_euler.setDate(date_0);
    hest_euler.rho = rho;
    hest_euler.kappa = kappa;
    hest_euler.theta = theta;
    hest_euler.xi = xi;
    
    CallOption c;
    c.setStrike(K);
    c.setMaturity(T);

    // Create the spot and vol initial normal and price paths
    std::vector<double> spot_draws(num_intervals, 0.0);  // Vector of initial spot normal draws
    std::vector<double> vol_draws(num_intervals, 0.0);   // Vector of initial correlated vol normal draws
    std::vector<double> spot_prices(num_intervals, S_0);  // Vector of initial spot prices
    std::vector<double> vol_prices(num_intervals, v_0);   // Vector of initial vol prices

    // Monte Carlo options pricing
    double payoff_sum = 0.0;
    //double payoff_sum1 = 0.0;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (unsigned i=0; i<num_sims; i++) {
      //std::cout << "Calculating path " << i+1 << " of " << num_sims << std::endl;
      fill_correlated_normals(vol_draws, spot_draws, rho);
      vol_prices = hest_euler.calc_vol_path(T, num_intervals, vol_draws);
      spot_prices = hest_euler.generatePricePath(T, num_intervals, vol_prices, spot_draws);
      payoff_sum += c.payoff(spot_prices);
        //payoff_sum1 += c.payoff(hest_euler.generatePricePath1(T, num_intervals, vol_prices, spot_draws));
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << " Option was priced in "<< pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
    double option_price = (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
    std::cout << " Option Price: " << option_price << std::endl;
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

    ScottChesney hest_euler;
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
    std::vector<double> vol_prices(num_intervals, log(v_0));   // Vector of initial vol prices

    // Monte Carlo options pricing
    double payoff_sum = 0.0;
    for (unsigned i=0; i<num_sims; i++) {
      //std::cout << "Calculating path " << i+1 << " of " << num_sims << std::endl;
      fill_correlated_normals(vol_draws, spot_draws, rho);
      //vol_prices = hest_euler.calc_vol_path(T, num_intervals, vol_draws);
      spot_prices = hest_euler.generatePricePath(T, num_intervals, vol_prices, spot_draws);
      payoff_sum += c.payoff(spot_prices);
    }
    double option_price = (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
    
    BlackScholesModel bsm;
    bsm.setRiskFreeRate(r);
    bsm.setVolatility(v_0);
    bsm.setStockPrice(S_0);
    bsm.setDate(date_0);
    
    ASSERT_APPROX_EQUAL(option_price, c.price(bsm), 0.1);
}

void testScottChesneyModel(){
    TEST( testScottChesneyPayoff );
    TEST( testScottChesneyVisually );
    TEST( testCallPrice );
}
