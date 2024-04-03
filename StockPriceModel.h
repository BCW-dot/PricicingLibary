#pragma once
#include "stdafx.h"

class StockPriceModel{
public:
    virtual ~StockPriceModel() {}
    /* generate P-measure price paths */
    virtual std::vector<double> generatePricePath(double toDate, int nSteps) const =0;
    
    double getDrift() const {
        return drift;
    }

    void setDrift( double drift ) {
        this->drift = drift;
    }

    double getStockPrice() const {
        return stockPrice;
    }
    
    void setStockPrice( double stockPrice ) {
        this->stockPrice = stockPrice;
    }
    
    double getVolatility()  const {
        return volatility;
    }

    void setVolatility( double volatility ) {
        this->volatility = volatility;
    }

    double getRiskFreeRate()  const {
        return riskFreeRate;
    }
    
    void setRiskFreeRate( double riskFreeRate) {
        this->riskFreeRate = riskFreeRate;
    }
    
    double getDate() const {
        return date;
    }
    
    void setDate( double date ) {
        this->date = date;
    }
    
private:
    double drift;
    double stockPrice;
    double volatility;
    double riskFreeRate;
    double date;
};
