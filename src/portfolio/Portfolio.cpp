#include "PricingLibrary/portfolio/Portfolio.h"
#include "PricingLibrary/options/CallOption.h"
#include "PricingLibrary/options/PutOption.h"
#include "PricingLibrary/options/UpAndInOption.h"
#include "PricingLibrary/options/UpAndOutOption.h"
#include "PricingLibrary/utils/ZeroCouponBond.h"
#include "PricingLibrary/utils/Stock.h"

using namespace std;

/*  
 *  By using an abstract interface class with a factory constructor
 *  and only having the implementation in the C++ file we increase information
 *  hiding. Nobody knows about the PortfolioImpl class outside of the C++ file,
 *  so we can change it without any impact on anything else.
 */
class PortfolioImpl : public Portfolio {
public:
    /*  Returns the number of items in the portflio */
    int size() const;    
    /*  Add a new security to the portfolio, returns the index
        at which it was added */
    int add( double quantity,
             shared_ptr<Priceable> security );
    /*  Update the quantity at a given index */
    void setQuantity( int index, double quantity );
    /*  Compute the current price */
    double price( const BlackScholesModel& model ) const;    

//private:
    vector<double> quantities;
    vector< shared_ptr<Priceable> > securities;
};

int PortfolioImpl::size() const {
    return quantities.size();
}

int PortfolioImpl::add( double quantity,
            shared_ptr<Priceable> security ) {
    quantities.push_back( quantity );
    securities.push_back( security );
    return quantities.size();
}

double PortfolioImpl::price(
        const BlackScholesModel& model ) const {
    double ret = 0;
    int n = size();
    for (int i=0; i<n; i++) {
        ret += quantities[i] * securities[i]->price( model );
    }
    return ret;
}

void PortfolioImpl::setQuantity( int index,
        double quantity ) {
    quantities[index] = quantity;
}

/**
 *   Create a Portfolio
 */
shared_ptr<Portfolio> Portfolio::newInstance() {
    shared_ptr<Portfolio> ret=make_shared<PortfolioImpl>();
    return ret;
}


/////////////////////////////
//  Tests
/////////////////////////////

static void testSingleSecurity() {
    shared_ptr<Portfolio> portfolio = Portfolio::newInstance();
    
    shared_ptr<CallOption> c=make_shared<CallOption>();
    c->setStrike(110);
    c->setMaturity(1.0);

    portfolio->add( 100, c );

    BlackScholesModel bsm;
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100);
    bsm.setRiskFreeRate(0);
    
    double unitPrice = c->price( bsm );
    double portfolioPrice = portfolio->price( bsm );
    ASSERT_APPROX_EQUAL( 100*unitPrice, portfolioPrice, 0.0001);
}

static void testPutCallParity() {
    shared_ptr<Portfolio> portfolio
        = Portfolio::newInstance();
    
    shared_ptr<CallOption> c=make_shared<CallOption>();
    c->setStrike(110);
    c->setMaturity(1.0);

    shared_ptr<PutOption> p=make_shared<PutOption>();
    p->setStrike(110);
    p->setMaturity(1.0);


    portfolio->add( 100, c );
    portfolio->add( -100, p );

    BlackScholesModel bsm;
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100);
    bsm.setRiskFreeRate(0);


    double expected = bsm.getStockPrice() - c->getStrike();
    double portfolioPrice = portfolio->price( bsm );
    
    ASSERT_APPROX_EQUAL(100*expected,portfolioPrice,0.0001);

}

static void testStocksAndBonds() {

	shared_ptr<Stock> stock = make_shared<Stock>();
	shared_ptr<ZeroCouponBond> bond = make_shared<ZeroCouponBond>();
	bond->maturity = 1.0;

	shared_ptr<Portfolio> portfolio
		= Portfolio::newInstance();
	portfolio->add(1, stock);
	portfolio->add(1, bond);
	
    BlackScholesModel bsm;
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100);
    bsm.setRiskFreeRate(0);

	double portfolioPrice = portfolio->price(bsm);

	ASSERT_APPROX_EQUAL(101.0, portfolioPrice, 0.0001);

}

static void testBarrierOptions() {

	double strike = 100;
	double barrier = 105;
	double maturity = 1.0;
	shared_ptr<UpAndOutOption> upAndOut = make_shared<UpAndOutOption>();
	shared_ptr<UpAndInOption> upAndIn = make_shared<UpAndInOption>();
	shared_ptr<CallOption> option = make_shared<CallOption>();
    
	option->setStrike(strike);
    
	upAndOut->setStrike(strike);
	upAndOut->setBarrier(barrier);
    
	upAndIn->setStrike(strike);
	upAndIn->setBarrier(barrier);
    
	option->setMaturity(maturity);
    
	upAndOut->setMaturity(maturity);
	upAndIn->setMaturity(maturity);

	shared_ptr<Portfolio> portfolio = Portfolio::newInstance();
	portfolio->add(1, upAndOut);
	portfolio->add(1, upAndIn);
	portfolio->add(-1, option);

    BlackScholesModel bsm;
    bsm.setVolatility(0.1);
    bsm.setStockPrice(100);
    bsm.setRiskFreeRate(0);

	double portfolioPrice = portfolio->price(bsm);

	ASSERT_APPROX_EQUAL(0.0, portfolioPrice, 0.1);

}


void testPortfolio() {
    TEST( testSingleSecurity );
    TEST( testPutCallParity );
	TEST( testStocksAndBonds );
	TEST(testBarrierOptions);
}

