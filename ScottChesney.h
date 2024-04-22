#pragma once

#include "stdafx.h"
#include "StockPriceModel.h"


class ScottChesney : public StockPriceModel{
public:
    std::vector<double> generatePricePath(double toDate, int nSteps) const;
    //std::vector<double> generateRiskNeutralPricePath(double toDate, int nSteps) const;
    
    std::vector<double> generatePricePath(double toDate, int nSteps, std::vector<double>& vol_path, std::vector<double>& spot_path) const;
    //std::vector<double> generateRiskNeutralpricePath
    std::vector<double> calc_vol_path(double toDate, int nSteps, const std::vector<double>& vol_draws) const;
    
    double kappa;
    double theta;
    double xi;
    double rho;
};

void testScottChesneyModel();
