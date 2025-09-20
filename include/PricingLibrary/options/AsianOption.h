#pragma once

#include "ContinuousTimeOptionBase.h"

/**
 *  This class represents an Asian call option
 */
class AsianOption : public ContinuousTimeOptionBase {
public:
	virtual double payoff(const std::vector<double>& stockPrices) const;

	bool isPathDependent() const {
		return true;
	}
};

void testAsianOption();
