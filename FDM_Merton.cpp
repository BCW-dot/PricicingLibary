#include "fdm.h"
#include "stdafx.h"
#include "CallOption.h"
#include "BlackScholesModel.h"
#include "LineChart.h"
#include "matlib.h"

#include <fstream>

using namespace std;
#include "FDM_Merton.h"

FDM_Merton::FDM_Merton(double _x_dom, unsigned long _J,
                       double _t_dom, unsigned long _N, double extraP,
                       std::shared_ptr<ConvectionDiffusionPDE> _pde)
  : FDMBase(_x_dom, _J, _t_dom, _N, _pde) {
  calculate_step_sizes();
  set_initial_conditions();
      extraP = extraP;
}

void FDM_Merton::calculate_step_sizes() {
  dx = x_dom/static_cast<double>(J-1);
  dt = t_dom/static_cast<double>(N-1);
}

void FDM_Merton::set_initial_conditions() {
  // Spatial settings
  double cur_spot = 0.0;

  old_result.resize(J, 0.0);
  new_result.resize(J, 0.0);
  x_values.resize(J, 0.0);

  for (unsigned long j=0; j<J; j++) {
    cur_spot = j*dx;
    old_result[j] = pde->init_cond(cur_spot);
    x_values[j] = cur_spot;
  }

  // Temporal settings
  prev_t = 0.0;
  cur_t = 0.0;
}

//!!!!!!!CHECK THIS!!!!!!!
void FDM_Merton::calculate_boundary_conditions() {
  new_result[0] = pde->boundary_left(prev_t, x_values[0]);
    
    for (int i = J - extraP - 1; i < J; ++i) {
        new_result[i] = pde->boundary_right(prev_t, x_values[i]);// exp(x[i]) - K * exp(-r * j * dt);
    }
    for (int i = 0; i <= extraP; ++i) {
        new_result[i] = pde->boundary_left(prev_t, x_values[0]);
    }
    
  new_result[J-(extraP-1)] = pde->boundary_right(prev_t, x_values[J-(extraP-1)]);
}

void FDM_Merton::calculate_inner_domain() {
  // Only use inner result indices (1 to J-2)
  for (unsigned long j=1; j<J-1; j++) {
    // Temporary variables used throughout
    double dt_sig = dt * (pde->diff_coeff(prev_t, x_values[j]));
    double dt_sig_2 = dt * dx * 0.5 * (pde->conv_coeff(prev_t, x_values[j]));

    // Differencing coefficients (see \alpha, \beta and \gamma in text)
      alpha = 0;
    beta = dx * dx - (2.0 * dt_sig) + (dt * dx * dx * (pde->zero_coeff(prev_t, x_values[j])));
    gamma = dt_sig + dt_sig_2;

    // Update inner values of spatial discretisation grid (Explicit Euler)
      new_result[j] = ( (alpha * old_result[j-1]) +
                      (beta * old_result[j]) +
                      (gamma * old_result[j+1]) )/(dx*dx) -
      (dt*(pde->source_coeff(prev_t, x_values[j])));
  }
}
void FDM_Merton::step_march() {
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

void testFDM_MertonPDE(){
}
