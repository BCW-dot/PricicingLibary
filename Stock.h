#pragma once

#include "Priceable.h"

class Stock : public Priceable {
public:
	/*  Compute the price of the security in the
	Black Scholes world */
	virtual double price( const BlackScholesModel& model ) const {
		return model.getStockPrice();
	}
};
