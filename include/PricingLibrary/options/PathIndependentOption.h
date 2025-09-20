#pragma once

#include "../stdafx.h"
#include "ContinuousTimeOptionBase.h"

/**
 *   This states that all path independent options
 *   have a payoff determined by the final stock price
 */
class PathIndependentOption : public ContinuousTimeOptionBase {
public:
    /*  A virtual destructor */
    virtual ~PathIndependentOption() {}
    /*  Returns the payoff at maturity */
    virtual double payoff( double finalStockPrice) const
        = 0;
    /*  Compute the payoff from a price path */
    double payoff(
            const std::vector<double>& stockPrices ) const {
        return payoff(stockPrices.back());
    }
    /*  Is the option path dependent? */
    bool isPathDependent() const {
        return false;
    };
};
