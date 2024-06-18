#pragma once

#include "stdafx.h"
#include "StockPriceModel.h"


class HestonModel : public StockPriceModel{
public:
    std::vector<double> generatePricePath(double toDate, int nSteps) const;
    std::vector<double> generateRiskNeutralPricePath(double toDate, int nSteps) const;
    
   /* This shoudl only be used for plotting */
    std::vector<double> getVolPath(double toDate, int nSteps, const std::vector<double>& vol_draws){
        return calc_vol_path(toDate, nSteps, vol_draws);
    }
    
    std::vector<double> getPath(double toDate, int nSteps, std::vector<double>& vol_path, std::vector<double>& spot_path){
        return generatePricePath(toDate, nSteps, vol_path, spot_path);
    }
    
    double kappa;
    double theta;
    double xi;
    double rho;
    
    std::vector<double> independent_price_path(double toDate, int nSteps) const;
    std::vector<double> Q_measure_price_path(double toDate, int nSteps) const;
    
private:
    std::vector<double> calc_vol_path(double toDate, int nSteps, const std::vector<double>& vol_draws) const;
    
    std::vector<double> generatePricePath(double toDate, int nSteps, std::vector<double>& vol_path, std::vector<double>& spot_path) const;
    
    std::vector<double> generateRiskNeutralPricePath(double toDate, int nSteps, std::vector<double>& vol_path, std::vector<double>& spot_path) const;
};

void testHestonModel();
