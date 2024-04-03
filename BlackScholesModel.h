#pragma once

#include "stdafx.h"
#include "StockPriceModel.h"

class BlackScholesModel : public StockPriceModel {
public:
    std::vector<double> generatePricePath(double toDate, int nSteps) const;

    std::vector<double> generateRiskNeutralPricePath(double toDate, int nSteps) const;
    
    std::vector<double> generateRiskNeutralDeltaPath(double toDate, int nSteps, double price, std::vector<double> epsilon) const;

private:
    std::vector<double> generatePricePath(double toDate, int nSteps, double drift) const;

    std::vector<double> generateDeltaPath(double toDate, int nSteps, double drift, double price, std::vector<double> epsilon) const;
};

void testBlackScholesModel();
