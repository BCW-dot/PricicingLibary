#pragma once

#include "BarrierOption.h"

class UpAndInOption : public BarrierOption {
public:
	double payoff(
		const std::vector<double>& prices) const;
};


void testUpAndInOption();