#pragma once

#include "BarrierOption.h"

class UpAndOutOption : public BarrierOption {
public:
    double payoff(
        const std::vector<double>& prices ) const;
};


void testUpAndOutOption();