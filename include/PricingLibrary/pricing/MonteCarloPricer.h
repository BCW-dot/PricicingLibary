#pragma once

#include "../stdafx.h"
#include "../options/ContinuousTimeOption.h"
#include "../models/BlackScholesModel.h"

class MonteCarloPricer {
public:
    /*  Constructor */
    MonteCarloPricer();
    /*  Number of scenarios */
    int nScenarios;
    /*  The number of steps in the calculation */
    int nSteps;
    /*  Price the option */
    double price( const ContinuousTimeOption& option,
                  const BlackScholesModel& model );
    /* Delta of an Option */
    double delta(const ContinuousTimeOption& option, const BlackScholesModel &model);
};

void testMonteCarloPricer();

