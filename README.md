# PricingLibrary

A C++ learning project focused on understanding virtual dispatch and inheritance through financial derivatives pricing.

## Purpose

This project was created to gain hands-on experience with:
- C++ virtual functions and polymorphism
- Object-oriented design patterns
- Basic quantitative finance concepts
- CMake build systems

## Architecture

The library demonstrates inheritance through a simple pricing framework:

```
Priceable (abstract base)
    ↓
PathIndependentOption → CallOption, PutOption
    ↓
ContinuousTimeOptionBase → AsianOption, BarrierOptions

StockPriceModel (abstract base)
    ↓
BlackScholesModel, TwoLevelModel
```

## Features

- **Options**: Call, Put, Asian, and Barrier options with virtual dispatch
- **Models**: Black-Scholes and two-level stock price models
- **Pricing**: Monte Carlo simulation and finite difference methods
- **Utilities**: Basic mathematical functions and charting

## Building

```bash
mkdir build && cd build
cmake ..
make
```

## Example Usage

```cpp
#include "PricingLibrary/options/CallOption.h"
#include "PricingLibrary/models/BlackScholesModel.h"

BlackScholesModel model;
model.setStockPrice(100.0);
model.setVolatility(0.2);
model.setRiskFreeRate(0.05);

CallOption call;
call.setStrike(105.0);
call.setMaturity(1.0);

double price = call.price(model); // Virtual dispatch in action
```

## Project Structure

```
├── include/PricingLibrary/    # Headers organized by component
│   ├── models/               # Stock price models
│   ├── options/              # Option types (inheritance hierarchy)
│   ├── pricing/              # Pricing engines
│   └── utils/                # Math utilities
├── src/                      # Implementation files
├── examples/                 # Demo application
└── tests/                    # Test cases
```

## Note

This is a learning project created to explore C++ virtual dispatch patterns in a financial context. The implementations prioritize educational value over production-ready numerical accuracy.