#include "fdm.h"
#include "stdafx.h"
#include "CallOption.h"
#include "BlackScholesModel.h"
#include "LineChart.h"
#include "matlib.h"

#include <fstream>

using namespace std;

FDMBase::FDMBase(double _x_dom, unsigned long _J,
                 double _t_dom, unsigned long _N,
                 std::shared_ptr<ConvectionDiffusionPDE> _pde)
  : x_dom(_x_dom), J(_J), t_dom(_t_dom), N(_N), pde(_pde) {}

FDMEulerExplicit::FDMEulerExplicit(double _x_dom, unsigned long _J,
                                   double _t_dom, unsigned long _N,
                                   std::shared_ptr<ConvectionDiffusionPDE> _pde)
  : FDMBase(_x_dom, _J, _t_dom, _N, _pde) {
  calculate_step_sizes();
  set_initial_conditions();
}

void FDMEulerExplicit::calculate_step_sizes() {
  dx = x_dom/static_cast<double>(J-1);
  dt = t_dom/static_cast<double>(N-1);
}

void FDMEulerExplicit::set_initial_conditions() {
  // Spatial settings
  double cur_spot = 0.0;

  old_result.resize(J, 0.0);
  new_result.resize(J, 0.0);
  x_values.resize(J, 0.0);

  for (unsigned long j=0; j<J; j++) {
    cur_spot = static_cast<double>(j)*dx;
    old_result[j] = pde->init_cond(cur_spot);
    x_values[j] = cur_spot;
  }

  // Temporal settings
  prev_t = 0.0;
  cur_t = 0.0;
}

void FDMEulerExplicit::calculate_boundary_conditions() {
  new_result[0] = pde->boundary_left(prev_t, x_values[0]);
  new_result[J-1] = pde->boundary_right(prev_t, x_values[J-1]);
}

void FDMEulerExplicit::calculate_inner_domain() {
  // Only use inner result indices (1 to J-2)
  for (unsigned long j=1; j<J-1; j++) {
    // Temporary variables used throughout
    double dt_sig = dt * (pde->diff_coeff(prev_t, x_values[j]));
    double dt_sig_2 = dt * dx * 0.5 * (pde->conv_coeff(prev_t, x_values[j]));

    // Differencing coefficients (see \alpha, \beta and \gamma in text)
    alpha = dt_sig - dt_sig_2;
    beta = dx * dx - (2.0 * dt_sig) + (dt * dx * dx * (pde->zero_coeff(prev_t, x_values[j])));
    gamma = dt_sig + dt_sig_2;

    // Update inner values of spatial discretisation grid (Explicit Euler)
      new_result[j] = ( (alpha * old_result[j-1]) +
                      (beta * old_result[j]) +
                      (gamma * old_result[j+1]) )/(dx*dx) -
      (dt*(pde->source_coeff(prev_t, x_values[j])));
  }
}
void FDMEulerExplicit::step_march() {
  std::ofstream fdm_out("/Users/benewilkens/Documents/ArmstrongFin/CSV/fdm.csv");

  while(cur_t < t_dom) {
    cur_t = prev_t + dt;
    calculate_boundary_conditions();
    calculate_inner_domain();
    for (int j=0; j<J; j++) {
      fdm_out << x_values[j] << " " << prev_t << " " << new_result[j] << std::endl;
    }
    
    old_result = new_result;
    prev_t = cur_t;
  }

  fdm_out.close();
}

//DOESNT work wrll cince the threads are all accecsing the same functions, which are not parallel somehwoe
void FDMEulerExplicit::calculate_inner_domain_parallel() {
    #pragma omp parallel for
    for (unsigned long j = 1; j < J - 1; j++) {
        // Temporary variables used throughout
        double dt_sig = dt * (pde->diff_coeff(prev_t, x_values[j]));
        double dt_sig_2 = dt * dx * 0.5 * (pde->conv_coeff(prev_t, x_values[j]));

        // Differencing coefficients (see \alpha, \beta and \gamma in text)
        double alpha = dt_sig - dt_sig_2;
        double beta = dx * dx - (2.0 * dt_sig) + (dt * dx * dx * (pde->zero_coeff(prev_t, x_values[j])));
        double gamma = dt_sig + dt_sig_2;

        // Update inner values of spatial discretisation grid (Explicit Euler)
        double new_value = ((alpha * old_result[j - 1]) +
                            (beta * old_result[j]) +
                            (gamma * old_result[j + 1])) / (dx * dx) -
                           (dt * (pde->source_coeff(prev_t, x_values[j])));

        // Protect the update of new_result[j] using critical section
        #pragma omp critical
        {
            new_result[j] = new_value;
        }
    }
}
void FDMEulerExplicit::step_march_parallel() {
    std::ofstream fdm_out("/Users/benewilkens/Documents/ArmstrongFin/CSV/fdm.csv");

  while(cur_t < t_dom) {
    cur_t = prev_t + dt;
    calculate_boundary_conditions();
    calculate_inner_domain_parallel();
    for (int j=0; j<J; j++) {
      fdm_out << x_values[j] << " " << prev_t << " " << new_result[j] << std::endl;
    }
    
    old_result = new_result;
    prev_t = cur_t;
  }

  fdm_out.close();
}

void FDMEulerExplicit::copy_result(std::vector<double>& v){
    for(int i = 0; i<this->new_result.size(); i++){
        v[i] = this->new_result[i];
    }
}

static void testPDEvisually(){
    std::shared_ptr<CallOption> c = std::make_shared<CallOption>();
    c->setStrike(50.5);
    c->setMaturity(1.0);
    
    std::shared_ptr<BlackScholesModel> bsm = std::make_shared<BlackScholesModel>();
    bsm->setRiskFreeRate(0.05);
    bsm->setVolatility(0.2);

      // FDM discretisation parameters
      double x_dom = 100.0;       // Spot goes from [0.0, x_dom]
      unsigned long J = 40;       //1/J space discretization
      double t_dom = c->getMaturity();         // Time period as for the option
      unsigned long N = 100;       //1/N time discretization
        
      // Create the PDE and FDM objects
      std::shared_ptr<BlackScholesPDE> bs_pde = std::make_shared<BlackScholesPDE>(c, bsm);
      FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, bs_pde);

      // Run the FDM solver
    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    fdm_euler.step_march();
    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //std::cout << option_value.size() << " Options were priced in "<< pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    //std::cout << "Time difference = " << pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
    
    std::vector<double> option_value(J,0.0);
    fdm_euler.copy_result(option_value);
    double dx = x_dom/J;
    //std::cout << "Option price at t_0 for different S_0: ";
    //print(option_value);
    
    vector<double> spots = linspace(dx, x_dom, J );
    LineChart linechart;
    linechart.setSeries(spots, option_value);
    linechart.setTitle("x-Axis = S_0 spots; y-Axis = option value at t_0");
    linechart.writeAsHTML("PDE_call_prices.html");
    
    //plot("S0_call_prices.html",spots, option_value );
}

static void testPDEcallPricesVsBSMprices(){
    std::shared_ptr<CallOption> c = std::make_shared<CallOption>();
    c->setStrike(0.5);
    c->setMaturity(1.0);
    
    std::shared_ptr<BlackScholesModel> bsm = std::make_shared<BlackScholesModel>();
    bsm->setRiskFreeRate(0.05);
    bsm->setVolatility(0.2);

      // FDM discretisation parameters
      double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
      unsigned long J = 20;
      double t_dom = c->getMaturity();         // Time period as for the option
      unsigned long N = 20;
        
      // Create the PDE and FDM objects
      std::shared_ptr<BlackScholesPDE> bs_pde = std::make_shared<BlackScholesPDE>(c, bsm);
      FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, bs_pde);

    // Run the FDM solver
    fdm_euler.step_march();
    
    std::vector<double> option_value(J,0.0);
    fdm_euler.copy_result(option_value);
    double dx = x_dom/J;
    //std::cout << "Option price at t_0 for different S_0: ";
    //print(option_value);
    
    vector<double> spots = linspace(dx, x_dom, J );
    vector<double> call_prices_closed_forsolution(J,0.0);
    
    CallOption option;
    option.setStrike(c->getStrike());
    option.setMaturity(c->getMaturity());
    
    BlackScholesModel m;
    m.setRiskFreeRate(bsm->getRiskFreeRate());
    m.setVolatility(bsm->getVolatility());
    m.setStockPrice(0.0);
    m.setDate(0.0);
    
    for(int i = 0; i<J; i++){
        double answer = option.price(m);
        call_prices_closed_forsolution[i] = answer;
        ASSERT_APPROX_EQUAL(call_prices_closed_forsolution[i], option_value[i], 0.05);
        m.setStockPrice(dx+i*dx);
    }
    
    plot("testCallprices.html", spots, call_prices_closed_forsolution );
    plot("CallPDEprices.html", spots, option_value);
}

/*
static void testPDEparallel(int M){
    std::shared_ptr<CallOption> c = std::make_shared<CallOption>();
    c->setStrike(50.5);
    c->setMaturity(1.0);
    
    std::shared_ptr<BlackScholesModel> bsm = std::make_shared<BlackScholesModel>();
    bsm->setRiskFreeRate(0.05);
    bsm->setVolatility(0.2);

      // FDM discretisation parameters
      double x_dom = 100.0;       // Spot goes from [0.0, 1.0]
      unsigned long J = M;
      double t_dom = c->getMaturity();         // Time period as for the option
      unsigned long N = 100;
        
      // Create the PDE and FDM objects
      std::shared_ptr<BlackScholesPDE> bs_pde = std::make_shared<BlackScholesPDE>(c, bsm);
      FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, bs_pde);
    
    // Run the parallel FDM solver
  std::cout << "\n";
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  fdm_euler.step_march_parallel();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Parallel PDE solver time: " << pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
    
}
*/

static void testPDEspeed(){
    std::shared_ptr<CallOption> c = std::make_shared<CallOption>();
    c->setStrike(50.5);
    c->setMaturity(1.0);
    
    std::shared_ptr<BlackScholesModel> bsm = std::make_shared<BlackScholesModel>();
    bsm->setRiskFreeRate(0.05);
    bsm->setVolatility(0.2);

      // FDM discretisation parameters
      double x_dom = 100.0;       // Spot goes from [0.0, x_dom]
      unsigned long J = 40;
      double t_dom = c->getMaturity();         // Time period as for the option
      unsigned long N = 100;
        
      // Create the PDE and FDM objects
      std::shared_ptr<BlackScholesPDE> bs_pde = std::make_shared<BlackScholesPDE>(c, bsm);
      FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, bs_pde);
    
    // Run the FDM solver
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    fdm_euler.step_march();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << J << " Options were priced in "<< pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
}

/*
static void testSimpleVsParallel(){
    int J = 40;
    //testPDEparallel(J);
    testPDEsimple(J);
}
*/

void testBlackScholesPDE(){
    TEST( testPDEvisually );
    TEST( testPDEspeed );
    TEST( testPDEcallPricesVsBSMprices );
}
