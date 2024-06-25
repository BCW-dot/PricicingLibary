#pragma once
#include "stdafx.h"
#include "PathIndependentOption.h"
#include "BlackScholesModel.h"

#include "DiffsuionPDE.h"

class MertonPDE : public ConvectionDiffusionPDE {
 public:
  std::shared_ptr<PathIndependentOption> option;
  std::shared_ptr<BlackScholesModel> model;
    
    MertonPDE(std::shared_ptr<PathIndependentOption> option, std::shared_ptr<BlackScholesModel> model,
              double poisson_parameter, double mean_jump_size, double variance_jump_size, double m_hat, double lamda_hat, double nu);//, double extraP);
  
  double diff_coeff(double t, double x) const;
  double conv_coeff(double t, double x) const;
  double zero_coeff(double t, double x) const;
  double source_coeff(double t, double x) const;
  double integral_item(double t, double x) const;

  double boundary_left(double t, double x) const;
  double boundary_right(double t, double x) const;
  double init_cond(double x) const;
    
    double poisson_parameter;
    double mean_jump_size;
    double variance_jump_size;
    double m_hat;
    double lamda_hat;
    std::vector<double> nu;
    //double extraP;
};

void testMertonPDE();
