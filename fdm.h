#pragma once

#include "pde.h"
#include <vector>
#include "FdmBase.h"

class FDMEulerExplicit : public FDMBase {
 protected:
  void calculate_step_sizes();
  void set_initial_conditions();
  void calculate_boundary_conditions();
  void calculate_inner_domain();
    void calculate_inner_domain_parallel();

 public:
  FDMEulerExplicit(double _x_dom, unsigned long _J,
                   double _t_dom, unsigned long _N,
                   std::shared_ptr<ConvectionDiffusionPDE> _pde);

  void step_march();
  void step_march_parallel();
    
  void copy_result(std::vector<double>& v);
};

void testBlackScholesPDE();

/*
std::ofstream fdm_out("/Users/benewilkens/Documents/ArmstrongFin/CSV/fdm.csv");
void FDMEulerExplicit::copy_result(std::vector<double>& v){
    for(int i = 0; i<this->new_result.size(); i++){
        v[i] = this->new_result[i];
    }
}


static void testBSMpde_visually(){
   // Create the PayOff and Option objects
    shared_ptr<CallOption> c = make_shared<CallOption>();
    c->setStrike(0.5);
    c->setMaturity(1.0);

    shared_ptr<BlackScholesModel> bsm = make_shared<BlackScholesModel>();
    bsm->setVolatility(0.2);
    bsm->setRiskFreeRate(0.05);
    
    // FDM discretisation parameters
    double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
    unsigned long J = 20;//40;
    double t_dom = c->getMaturity();    // Time period as for the option
    unsigned long N = 20;//100;
      
    // Create the PDE and FDM objects
    shared_ptr<BlackScholesPDE> bs_pde = make_shared<BlackScholesPDE>(bsm, c);
    FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, bs_pde);
    // Run the FDM solver
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    fdm_euler.step_march();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    std::cout << "Time difference = " << pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
    
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

static void testBSMpde_callPrices(){
    shared_ptr<CallOption> c = make_shared<CallOption>();
    c->setStrike(0.5);
    c->setMaturity(1.0);

    shared_ptr<BlackScholesModel> bsm = make_shared<BlackScholesModel>();
    bsm->setVolatility(0.2);
    //bsm->setStockPrice(1.0);
    bsm->setRiskFreeRate(0.05);
    
    // FDM discretisation parameters
    double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
    unsigned long J = 20;//40;
    double t_dom = c->getMaturity();    // Time period as for the option
    unsigned long N = 20;//100;
      
    // Create the PDE and FDM objects
    shared_ptr<BlackScholesPDE> bs_pde = make_shared<BlackScholesPDE>(bsm, c);
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
    plot("testCallPDE.html", spots, call_prices_closed_forsolution );
}


void testBlackScholesPDE(){
    TEST( testBSMpde_visually );
    TEST( testBSMpde_callPrices );
}
 */


