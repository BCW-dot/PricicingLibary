#pragma once
#include "../stdafx.h"
#include "../options/PathIndependentOption.h"
#include "../models/BlackScholesModel.h"

// Convection Diffusion Equation - Second-order PDE
class ConvectionDiffusionPDE {
 public:
  // PDE Coefficients
  virtual double diff_coeff(double t, double x) const = 0;
  virtual double conv_coeff(double t, double x) const = 0;
  virtual double zero_coeff(double t, double x) const = 0;
  virtual double source_coeff(double t, double x) const = 0;

  // Boundary and initial conditions
  virtual double boundary_left(double t, double x) const = 0;
  virtual double boundary_right(double t, double x) const = 0;
  virtual double init_cond(double x) const = 0;
};

// Black-Scholes PDE
class BlackScholesPDE : public ConvectionDiffusionPDE {
 public:
  std::shared_ptr<PathIndependentOption> option;
  std::shared_ptr<BlackScholesModel> model;
    
  BlackScholesPDE(std::shared_ptr<PathIndependentOption> option, std::shared_ptr<BlackScholesModel> model);
  
  double diff_coeff(double t, double x) const;
  double conv_coeff(double t, double x) const;
  double zero_coeff(double t, double x) const;
  double source_coeff(double t, double x) const;

  double boundary_left(double t, double x) const;
  double boundary_right(double t, double x) const;
  double init_cond(double x) const;
};
