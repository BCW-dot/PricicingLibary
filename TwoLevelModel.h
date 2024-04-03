#pragma once

#include "stdafx.h"
#include "StockPriceModel.h"

/*
class TwoLevelModel : StockPriceModel{
public:
    TwoLevelModel();
    double volatility2;
    
    std::vector<double> generatePricePath(double toDate, int nSteps);
    
};
 */

class TwoLevelModel : public StockPriceModel {
public:
    TwoLevelModel();
    double volatility2;
    std::vector<double> generatePricePath(double toDate, int nSteps) const;
    std::vector<double> generateRiskNeutralPricePath(double toDate, int nSteps) const;
};

void testTwoLevelModel();
