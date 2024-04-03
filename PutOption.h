#pragma once

#include "stdafx.h"
#include "BlackScholesModel.h"
#include "PathIndependentOption.h"

class PutOption : public ContinuousTimeOptionBase {
public:

    /*  Calculate the payoff of the option given
        a history of prices */
    double payoff(
        const std::vector<double>& stockPrices
        ) const;

    double price( const BlackScholesModel& bsm )
        const;
    
    double delta( const BlackScholesModel& bsm ) const;

    bool isPathDependent() const {
        return false;
    };
};

void testPutOption();
