#include "StockPriceModel.h"

std::vector<double> StockPriceModel::generateRiskNeutralPricePath(double toDate, int nSteps) const{
    return generatePricePath(toDate, nSteps);
}

