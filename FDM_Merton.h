#pragma once

#include "MertonPDE.h"
#include <vector>
//#include "FdmBase.h"

class FDM_Merton : public FDMBase {
 protected:
  void calculate_step_sizes();
  void set_initial_conditions();
  void calculate_boundary_conditions();
  void calculate_inner_domain();
  void calculate_inner_domain_parallel();

 public:
  FDM_Merton(double _x_dom, unsigned long _J,
             double _t_dom, unsigned long _N, double extraP,
             std::shared_ptr<ConvectionDiffusionPDE> _pde);

  void step_march();
    
  void copy_result(std::vector<double>& v);
  double extraP;
};

void testFDM_MertonPDE();
