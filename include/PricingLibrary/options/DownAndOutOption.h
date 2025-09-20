#pragma once

#include "BarrierOption.h"

class DownAndOutOption : public BarrierOption {
public:
    double payoff(
        const std::vector<double>& prices ) const;
};


void testDownAndOutOption();