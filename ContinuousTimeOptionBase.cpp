#include "ContinuousTimeOptionBase.h"
#include "MonteCarloPricer.h"
#include "UpAndOutOption.h"
#include "CallOption.h"
#include "matlib.h"
#include "testing.h"

double ContinuousTimeOptionBase::price( const BlackScholesModel& model ) const {
    MonteCarloPricer pricer;
    return pricer.price( *this, model );
}

double ContinuousTimeOptionBase::delta( const BlackScholesModel& model ) const {
    MonteCarloPricer pricer;
    return pricer.delta( *this, model );
}

//////////////////////////////
//
//  Tests
//
//////////////////////////////

static void testPrice() {
    rng("default");

    BlackScholesModel bsm;
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100);
    bsm.setRiskFreeRate(0);

    CallOption callOption;
    callOption.setStrike(100);
    callOption.setMaturity(1.0);

    UpAndOutOption knockoutOption;
    knockoutOption.setStrike(100);
    knockoutOption.setMaturity(1.0);
    knockoutOption.setBarrier(1000);

    ContinuousTimeOptionBase& o1=callOption;
    ContinuousTimeOptionBase& o2=knockoutOption;

    double p1 = o1.price( bsm );
    double p2 = o2.price( bsm );

    ASSERT_APPROX_EQUAL( p1, p2, 0.1);
}

void testContinuousTimeOptionBase() {
    TEST( testPrice );
}
