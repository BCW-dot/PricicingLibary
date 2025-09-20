#include "PricingLibrary/pde/pde.h"
#include <math.h>

BlackScholesPDE::BlackScholesPDE(std::shared_ptr<PathIndependentOption> option, std::shared_ptr<BlackScholesModel> model) : option(option),
        model(model)
{}

// Diffusion coefficient
double BlackScholesPDE::diff_coeff(double t, double x) const {
  double vol = model->getVolatility();
  return 0.5*vol*vol*x*x;  // \frac{1}{2} \sigma^2 S^2
}

// Convection coefficient
double BlackScholesPDE::conv_coeff(double t, double x) const {
  return (model->getRiskFreeRate())*x;  // rS
}

// Zero-term coefficient
double BlackScholesPDE::zero_coeff(double t, double x) const {
  return -(model->getRiskFreeRate());  // -r
}

// Source coefficient
double BlackScholesPDE::source_coeff(double t, double x) const {
  return 0.0;
}

// Left boundary-condition (vanilla call option)
double BlackScholesPDE::boundary_left(double t, double x) const {
  return 0.0;  // Specifically for a CALL option
}

// Right boundary-condition (vanilla call option)
double BlackScholesPDE::boundary_right(double t, double x) const {
  // This is via Put-Call Parity and works for a call option
  return (x-(option->getStrike())*exp(-(model->getRiskFreeRate())*((option->getMaturity())-t)));
}

// Initial condition (vanilla call option)
double BlackScholesPDE::init_cond(double x) const {
  return option->payoff(x);
}


