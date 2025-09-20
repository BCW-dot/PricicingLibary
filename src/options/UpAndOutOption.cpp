#include "PricingLibrary/options/UpAndOutOption.h"
#include "PricingLibrary/options/BarrierOption.h"

using namespace std;

double UpAndOutOption::payoff(
        const std::vector<double>& prices ) const {
    int n = prices.size();
    for (int i=0; i<n; i++) {
        if (prices[i]>getBarrier()) {
            return 0.0;
        }
    }
    double stockAtMaturity = prices.back();
    if (stockAtMaturity>getStrike()) {
        return stockAtMaturity-getStrike();
    } else {
        return 0.0;
    }
}

/////////////////////////////////////
//
//   TESTS
//
/////////////////////////////////////


void testUpAndOutOption() {
    UpAndOutOption o;
    o.setBarrier(100);
    o.setStrike(70);
    vector<double> prices;
    prices.push_back(120);
    prices.push_back(80);
    ASSERT_APPROX_EQUAL( o.payoff( prices ), 0.0, 0.001);
    prices[0] = 90;
    ASSERT_APPROX_EQUAL( o.payoff( prices ), 10.0, 0.001);
    prices[1] = 60;
    ASSERT_APPROX_EQUAL( o.payoff( prices ), 0.0, 0.001);
}