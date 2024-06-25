#include "matlib.h"

#include "geometry.h"
#include "LineChart.h"
#include "Histogram.h"
 

using namespace std;

/*  Create a linearly spaced vector */
vector<double> linspace( double from, double to, int numPoints ) {
    ASSERT( numPoints>=2 );
    vector<double> ret(numPoints,0.0);
    double step = (to-from)/(numPoints-1);
    double current = from;
    for (int i=0; i<numPoints; i++) {
        ret[i]=current;
        current+=step;
    }
    return ret;
}

/**
 *  Find the sum of the elements in an array
 */
double sum( const std::vector<double>& v ) {
	double total = 0.0;
	int n = v.size();
	for (int i=0; i<n; i++) {
		total+= v[i];
	}
	return total;
}


/*  Compute the mean of a vector */
double mean( const vector<double>& v ) {
    int n = v.size();
    ASSERT( n > 0);
    return sum(v)/n;
}

/*  Compute the standard deviation of a vector */
double standardDeviation( const vector<double>& v, bool population ) {
    int n = v.size();
    double total = 0.0;
    double totalSq = 0.0;
    for (int i=0; i<n; i++) {
        total+=v[i];
        totalSq+=v[i]*v[i];
    }
    if (population)  {
        ASSERT( n > 0 );
        return sqrt( (totalSq - total*total/n)/n );
    } else {
        ASSERT( n > 1 );
        return sqrt( (totalSq - total*total/n)/(n-1) );
    }
}

/*  Find the minimum of a vector */
double min( const vector<double>& v ) {
    int n = v.size();
    ASSERT( n > 0);
    double min = v[0];
    for (int i=1; i<n; i++) {
        if (v[i]<min) {
            min=v[i];
        }
    }
    return min;
}

/*  Find the maximum of a vector */
double max( const vector<double>& v ) {
    int n = v.size();
    ASSERT( n > 0);
    double max = v[0];
    for (int i=1; i<n; i++) {
        if (v[i]>max) {
            max=v[i];
        }
    }
    return max;
}

/*  Create uniformly distributed random numbers using the C random number API*/
vector<double> randuniformOld( int n ) {
    vector<double> ret(n, 0.0);
    for (int i=0; i<n; i++) {
        int randInt = rand();
        ret[i] = (randInt + 0.5)/(RAND_MAX+1.0);
    }
    return ret;
}

/*  MersenneTwister random number generator */
static mt19937 mersenneTwister;

/*  Reset the random number generator. We've borrowed the library call
    from MATLAB, though we're ignoring the description string */
void rng( const string& description ) {
    ASSERT( description=="default" );
    mersenneTwister.seed(mt19937::default_seed);
}

void my_rng() {
    mersenneTwister.seed(time( NULL));
}

/*  Create uniformly distributed random numbers using 
    the Mersenne Twister algorithm.*/
vector<double> randuniform( int n ) {
    vector<double> ret(n, 0.0);
    for (int i=0; i<n; i++) {
        ret[i] = (mersenneTwister()+0.5)/(mersenneTwister.max()+1.0);
    }
    return ret;
}

vector<double> my_randuniform( int n ) {
    //mersenneTwister.seed( time( NULL) );
    vector<double> ret(n, 0.0);
    for (int i=0; i<n; i++) {
        ret[i] = (mersenneTwister()+0.5)/(mersenneTwister.max()+1.0);
    }
    return ret;
}


/*  Create normally distributed random numbers */
vector<double> randn( int n ) {
    vector<double> v=randuniform(n);
    for (int i=0; i<n; i++) {
        v[i] = norminv(v[i]);
    }
    return v;
}

vector<double> my_randn( int n ) {
    vector<double> v = my_randuniform(n);
    for (int i=0; i<n; i++) {
        v[i] = norminv(v[i]);
    }
    /*
    vector<double> v(n,0.0);
    std::default_random_engine generator;
    generator.seed( time(NULL) );
    std::normal_distribution<double> distribution(0.0,1.0);
    for (int i=0; i<n; i++) {
        v[i] = distribution(generator);
    }
     */
    return v;
}

/* Create Poisson distrubuted random numbers */
vector<double> randp( int n, double parameter){
    vector<double> v(n,0);
    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::poisson_distribution<> d(parameter);
    for(int i =0; i<n; i++){
        v[i] = d(gen);
    }
    return v;
}

void fill_correlated_normals(vector<double>& n1, vector<double>& n2, double rho){
    int n = n1.size();
    int m = n2.size();
    ASSERT(n = m);
    ASSERT(rho <= 1 && rho >=-1);
    
    vector<double> v1 = randn(n);
    vector<double> v2 = randn(n);
    
    vector<double> v3(n,0.0);
    
    for(int i=0; i<n; i++ ){
        v3[i] = rho*v1[i] + sqrt(1-rho*rho)*v2[i];
    }
    
    n1 = v1;
    n2 = v3;
}

void open_plot(std::string file_name){
    std::string htmlFilePath = "/Users/benewilkens/Documents/ArmstrongFin/Plots/" + file_name;
    //const char* command = "open %s";
    //std::string openCommand[200];
    //sprintf(openCommand, command, htmlFilePath);
    //system(openCommand);
    
    // Command to open the HTML file with the default browser
    //std::string command = "open " + htmlFilePath; // "open" is the command to open files in Mac
    // Execute the command
    //system(command.c_str());
    
    const char* command = "open %s"; // "open" is the command to open files in Mac
        
        // Create a buffer to hold the command with file path
        char openCommand[200]; // Adjust size as needed
        
        // Format the command with the file path
        sprintf(openCommand, command, htmlFilePath.c_str());
        
        // Execute the command
        system(openCommand);
}

void open_hist(std::string file_name){
    std::string htmlFilePath = "/Users/benewilkens/Documents/ArmstrongFin/Histograms/" + file_name;
    //const char* command = "open %s";
    //std::string openCommand[200];
    //sprintf(openCommand, command, htmlFilePath);
    //system(openCommand);
    
    // Command to open the HTML file with the default browser
    //std::string command = "open " + htmlFilePath; // "open" is the command to open files in Mac
    // Execute the command
    //system(command.c_str());
    
    const char* command = "open %s"; // "open" is the command to open files in Mac
        
        // Create a buffer to hold the command with file path
        char openCommand[200]; // Adjust size as needed
        
        // Format the command with the file path
        sprintf(openCommand, command, htmlFilePath.c_str());
        
        // Execute the command
        system(openCommand);
}

void print_window(std::vector<double> v, int start, int end){
    ASSERT(start >= 0);
    ASSERT(end < v.size());
    std::cout << "{";
    for(int i = start; i<end-1; i++){
        std::cout << v[i] << ", ";
    }
    std::cout << v[end] << "}" << std::endl;
}

/**
 *  Sort a vector of doubles
 */
std::vector<double> sort( const std::vector<double>& v ) {
	std::vector<double> copy(v);
	std::sort( copy.begin(), copy.end() );
	return copy;
}

/**
 *  Find the given percentile of a distribution
 */
double prctile( const std::vector<double>& v, double percentage ) {
	// See the text for a precise specification
	// 
	ASSERT( percentage >=0.0 );
	ASSERT( percentage <=100.0 );
	int n = v.size();
	vector<double> sorted = sort( v );

	int indexBelow = (int)(n* percentage/100.0 - 0.5);
	int indexAbove = indexBelow + 1;
	if (indexAbove > n-1 ) {		
		return sorted[n-1];
	} if (indexBelow<0) {
		return sorted[0];
	}
	double valueBelow = sorted[ indexBelow ];
	double valueAbove = sorted[ indexAbove ];
	double percentageBelow = 100.0*(indexBelow+0.5)/n;
	double percentageAbove = 100.0*(indexAbove+0.5)/n;
	if (percentage<=percentageBelow) {
		return valueBelow;
	} 
	if (percentage>=percentageAbove) {
		return valueAbove;
	}
	double correction = (percentage - percentageBelow)*(valueAbove-valueBelow)/(percentageAbove-percentageBelow);
	return valueBelow + correction;
}

/**
 *  Convenience method for generating plots
 */
void plot( const string& file,
           const vector<double>& x, 
           const vector<double>& y ) {
    LineChart lc;
    lc.setSeries(x,y);
    lc.writeAsHTML( file );
}

/**
 *  Convenience method for generating plots
 */
void hist( const string& file,
           const vector<double>& data, 
           int numBuckets ) {
    Histogram h;
    h.setData(data);
	h.setNumBuckets(numBuckets);
    h.writeAsHTML( file );
}

const double ROOT_2_PI = sqrt( 2.0 * PI );


static inline double hornerFunction( double x, double a0, double a1) {
	return a0 + x*a1;
}

static inline double hornerFunction( double x, double a0, double a1, double a2) {
	return a0 + x*hornerFunction( x, a1, a2);
}

static inline double hornerFunction( double x, double a0, double a1, double a2, double a3) {
	return a0 + x*hornerFunction( x, a1, a2, a3);
}

static inline double hornerFunction( double x, double a0, double a1, double a2, double a3, double a4) {
	return a0 + x*hornerFunction( x, a1, a2, a3, a4);
}

static inline double hornerFunction( double x, double a0, double a1, double a2, double a3, double a4,
					   double a5) {
	return a0 + x*hornerFunction( x, a1, a2, a3, a4, a5);
}

static inline double hornerFunction( double x, double a0, double a1, double a2, double a3, double a4,
					   double a5, double a6) {
	return a0 + x*hornerFunction( x, a1, a2, a3, a4, a5, a6);
}

static inline double hornerFunction( double x, double a0, double a1, double a2, double a3, double a4,
					   double a5, double a6, double a7) {
	return a0 + x*hornerFunction( x, a1, a2, a3, a4, a5, a6, a7);
}

static inline double hornerFunction( double x, double a0, double a1, double a2, double a3, double a4,
					   double a5, double a6, double a7, double a8) {
	return a0 + x*hornerFunction( x, a1, a2, a3, a4, a5, a6, a7, a8);
}

/**
 *  Arguably this is a little easier to read than the original normcdf
 *  function as it makes the use of horner's method obvious.
 */
double normcdf( double x ) {
    DEBUG_PRINT( "normcdf("<<x<<")");
	if (x<0) {
		return 1-normcdf(-x);
	}
	double k = 1/(1 + 0.2316419*x);
	double poly = hornerFunction(k,
								 0.0, 0.319381530, -0.356563782,
								 1.781477937,-1.821255978,1.330274429);
	double approx = 1.0 - 1.0/ROOT_2_PI * exp(-0.5*x*x) * poly;
	return approx;
}

/*  Constants required for Moro's algorithm */
static const double a0 = 2.50662823884;
static const double a1 = -18.61500062529;
static const double a2 = 41.39119773534;
static const double a3 = -25.44106049637;
static const double b1 = -8.47351093090;
static const double b2 = 23.08336743743;
static const double b3 = -21.06224101826;
static const double b4 = 3.13082909833;
static const double c0 = 0.3374754822726147;
static const double c1 = 0.9761690190917186;
static const double c2 = 0.1607979714918209;
static const double c3 = 0.0276438810333863;
static const double c4 = 0.0038405729373609;
static const double c5 = 0.0003951896511919;
static const double c6 = 0.0000321767881768;
static const double c7 = 0.0000002888167364;
static const double c8 = 0.0000003960315187;

double norminv( double x ) {
	// We use Moro's algorithm
    DEBUG_PRINT( "norminv(" << x <<")" );
	double y = x - 0.5;
	if (y<0.42 && y>-0.42) {
		double r = y*y;
        DEBUG_PRINT( "Case 1, r=" << r );
		return y*hornerFunction(r,a0,a1,a2,a3)/hornerFunction(r,1.0,b1,b2,b3,b4);
	} else {
		double r;
		if (y<0.0) {
			r = x;
		} else {
			r = 1.0 - x;
		}
        DEBUG_PRINT( "Case 2, r=" << r);
		double s = log( -log( r ));
		double t = hornerFunction(s,c0,c1,c2,c3,c4,c5,c6,c7,c8);
		if (x>0.5) {
			return t;
		} else {
			return -t;
		}
	}
}

void display_matrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << "\t";
        }
        cout << endl;
    }
    cout<< " " << endl;
}

// Function to mimic numpy's linspace
std::vector<double> linspace2(double start, double end, int num) {
    std::vector<double> linspaced;

    if (num == 0) {
        return linspaced;
    }
    if (num == 1) {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        //cout << start + i * delta << ", ";
        linspaced.push_back(start + i * delta);
    }

    // Ensure the last value is exactly 'end' to handle floating-point arithmetic issues
    //linspaced.back() = end;

    return linspaced;
}

double normal_cdf(double x, double mu, double sigma) {
    return 0.5 * erfc(-(x - mu) / (sigma * sqrt(2)));
}

double payoff(double S, double K) {
    return max(S - K, 0.0);
}

/**
 *   Evaluate an integral using the rectangle rule
 */
double integral( RealFunction& f,
                 double a,
                 double b,
                 int nPoints ) {
    double h = (b-a)/nPoints;
    double x = a + 0.5*h;
    double total = 0.0;
    for (int i=0; i<nPoints; i++) {
        double y = f.evaluate(x);
        total+=y;
        x+=h;
    }
    return h*total;
}

void print(std::vector<double> v){
    for(int i = 0; i < v.size(); i++){
        std::cout<< v[i] << ", ";
    }
    std::cout << " " << std::endl;
    //std::cout << v[v.size()] << std::endl;
}

double sum_parallel(std::vector<double>& v){
    double sum = 0;
    int size = int(v.size());
    #pragma omp parallel for reduction(+:totalSum)
    for (int i = 0; i < size; ++i) {
        sum += v[i];
    }
    return sum;
}

double sum_simple(std::vector<double> v){
    double sum = 0;
    int size = int(v.size());
    for(int i = 0; i < size; i++){
        sum += v[i];
    }
    return sum;
}

static void testParallelSum(){
    int n = 1000000;
    std::vector<double> v1(n,1);
    
    double answer_simple;
    double answer_parallel;
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    answer_simple = sum_simple(v1);
    
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    std::cout << "Time difference = " << pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
    std::cout<< "Sum simple = " << answer_simple << std::endl;
    
    /* PARALLEL test*/
    begin = std::chrono::steady_clock::now();
    
    answer_parallel = sum_parallel(v1);
    
    end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    std::cout << "Time difference = " << pow(10,-6)*std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[s]" << std::endl;
    std::cout<< "Sum parallel = " << answer_parallel;
}

void writeVectorsToCsv(const std::vector<double>& v1, const std::vector<double>& t, const std::vector<double>& v2, std::string& filename){
    //file = "/Users/benewilkens/Documents/ArmstrongFin/Plots/" + filename;
    std::ofstream fdm_out(filename);

    // Write data
    for (int j=0; j<t.size(); j++) {
      fdm_out << v1[j] << " " << t[j] << " " << v2[j] << std::endl;
    }
    fdm_out.close();
}

///////////////////////////////////////////////
//
//   TESTS
//
///////////////////////////////////////////////

static vector<double> createTestVector() {
    vector<double> v;
    v.push_back(1);
    v.push_back(5);
    v.push_back(3);
    v.push_back(9);
    v.push_back(7);
    return v;
}

static void testLinspace() {
    vector<double> result = linspace(1.0, 10.0, 4 );
    ASSERT_APPROX_EQUAL( result[0], 1.0, 0.001 );
    ASSERT_APPROX_EQUAL( result[1], 4.0, 0.001 );
    ASSERT_APPROX_EQUAL( result[2], 7.0, 0.001 );
    ASSERT_APPROX_EQUAL( result[3], 10.0, 0.001 );
}

static void testMean() {
    ASSERT_APPROX_EQUAL( mean( createTestVector() ), 5.0, 0.001);
}

static void testStandardDeviation() {
    ASSERT_APPROX_EQUAL( standardDeviation( createTestVector() ), 3.1623, 0.001);
    ASSERT_APPROX_EQUAL( standardDeviation( createTestVector(), true ), 2.8284, 0.001);
}

static void testMin() {
    ASSERT_APPROX_EQUAL( min( createTestVector() ), 1.0, 0.001);
}

static void testMax() {
    ASSERT_APPROX_EQUAL( max( createTestVector() ), 9.0, 0.001);
}

static void testRanduniform() {
    rng("default");
    vector<double> v = randuniform(1000);
    ASSERT( ((int)v.size())==1000 );
    ASSERT_APPROX_EQUAL( mean(v), 0.5, 0.1);
    ASSERT( max(v)<1.0);
    ASSERT( min(v)>0.0);
}

static void testRandn() {
    //rng("default");
    vector<double> v = randn(10000);
    std::cout << "Normal Mean: " << mean(v) << std::endl;
    ASSERT( ((int)v.size())==10000 );
    ASSERT_APPROX_EQUAL( mean(v), 0.0, 0.1);
    ASSERT_APPROX_EQUAL( standardDeviation(v), 1.0, 0.1);
}


static void testNormCdf() {
    ASSERT_APPROX_EQUAL( normcdf( 1.96 ), 0.975, 0.001 );
}

static void testNormInv() {
    ASSERT_APPROX_EQUAL( norminv( 0.975 ), 1.96, 0.01 );
}

static void testPrctile() {
	const vector<double> v = createTestVector();
	ASSERT_APPROX_EQUAL( prctile( v, 100.0 ), 9.0, 0.001 );
	ASSERT_APPROX_EQUAL( prctile( v, 0.0 ), 1.0, 0.001 );
	ASSERT_APPROX_EQUAL( prctile( v, 50.0 ), 5.0, 0.001 );
	ASSERT_APPROX_EQUAL( prctile( v, 17.0 ), 1.7, 0.001 );
	ASSERT_APPROX_EQUAL( prctile( v, 62.0 ), 6.2, 0.001 );
}

/*  To test the integral function, we need a function
    to integrate */
class SinFunction : public RealFunction {
    double evaluate( double x );
};

double SinFunction::evaluate( double x ) {
    return sin(x);
}

static void testIntegral() {
    SinFunction integrand;
    double actual = integral( integrand, 1, 3, 1000 );
    double expected = -cos(3.0)+cos(1.0);
    ASSERT_APPROX_EQUAL( actual, expected, 0.000001);        
}

/**
 *  When you create a small class like this, using
 *  nested classes is easier.
 */
static void testIntegralVersion2() {

    class Sin : public RealFunction {
        double evaluate( double x ) {
            return sin(x);
        }
    };

    Sin integrand;
    double actual = integral( integrand, 1, 3, 1000 );
    double expected = -cos(3.0)+cos(1.0);
    ASSERT_APPROX_EQUAL( actual, expected, 0.000001);        
}


void testMatlib() {
    
    TEST( testLinspace );
    TEST( testMean );
    TEST( testStandardDeviation );
    TEST( testMin );
    TEST( testMax );
    TEST( testRanduniform );
    TEST( testRandn );
    TEST( testNormInv );
    TEST( testNormCdf );
    TEST( testPrctile );
    TEST( testIntegral );
    TEST( testIntegralVersion2 );
     
    TEST( testParallelSum );
}
