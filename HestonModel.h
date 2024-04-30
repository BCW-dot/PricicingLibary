#pragma once

#include "stdafx.h"
#include "StockPriceModel.h"


class HestonModel : public StockPriceModel{
public:
    std::vector<double> generatePricePath(double toDate, int nSteps) const;
    std::vector<double> generatePricePath(double toDate, int nSteps, std::vector<double>& vol_path, std::vector<double>& spot_path) const;
    
    std::vector<double> generateRiskNeutralPricePath(double toDate, int nSteps) const;
    std::vector<double> generateRiskNeutralPricePath(double toDate, int nSteps, std::vector<double>& vol_path, std::vector<double>& spot_path) const;
    
    //std::vector<double> calc_vol_path(double toDate, int nSteps, const std::vector<double>& vol_draws) const;
    std::vector<double> getVolPath(double toDate, int nSteps, const std::vector<double>& vol_draws){
        return calc_vol_path(toDate, nSteps, vol_draws);
    }
    
    double kappa;
    double theta;
    double xi;
    double rho;
private:
    std::vector<double> calc_vol_path(double toDate, int nSteps, const std::vector<double>& vol_draws) const;
};

void testHestonModel();
