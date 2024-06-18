#include "Poisson.h"
#include "matlib.h"

using namespace std;

std::vector<double> generatePoisson_path(double start, double end, int n, double param){
    vector<double> x(n,0);
    
    double dt = (end-start)/n;
    vector<double> Z = randp(n,dt*param);
    for(int i = 0; i<n; i++){
        x[i+1] = x[i] + Z[i];
    }

    return x;
}

std::vector<double> generatePoisson_Martingale(double start, double end, int n, double param){
    vector<double> x(n,0);
    
    double dt = (end-start)/n;
    vector<double> Z = randp(n,dt*param);
    //print_window(Z, 0 , 10);
    for(int i = 0; i<n; i++){
        x[i+1] = x[i] - param*dt+ Z[i];
    }

    return x;
}

std::vector<double> generate_Merton_paths(double start, double end, int n, double s0, double param, double mu_j, double sigma_j, double r, double sigma, vector<double>& s){
    vector<double> x(n,0);
    //vector<double> s(n,0);
    
    double dt = (end-start)/n;
    x[0] = log(s0);
    s[0] = s0;
    
    double Eej = exp(mu_j + 0.5 * sigma_j * sigma_j);
    vector<double> ZPois = randp(n,param*dt);
    vector<double> Z = randn(n);
    vector<double> J = randn(n);
    for(int i =0; i<n;i++){
        J[i] = mu_j + sqrt(sigma_j)*J[i];
    }
    
    for(int i=0;i<n;i++){
        x[i+1] = x[i] + (r - param*(Eej-1) - 0.5*sigma*sigma)*dt + sigma*sqrt(dt)*Z[i] + J[i] * ZPois[i];
        s[i+1] = exp(x[i]);
    }

    return x;
}

static void testPoisson(){
    double param = 5;
    int n = 10000;
    vector<double> v = randp(n,param);
    
    std::cout << "testing mean: " << endl;
    ASSERT_APPROX_EQUAL(mean(v), param, 0.1);
    std::cout << "testing variance: " << endl;
    ASSERT_APPROX_EQUAL(standardDeviation(v)*standardDeviation(v), param, 0.1);
    
    hist("Poisson.html", v);
    //open_hist("Poisson.html");
    
    int n1 = 5;
    double start = 0;
    double end = 1;
    
    vector<double> v1 = randp(n1,param);
    double dt = (end - start)/n1;
    vector<double> times = linspace(dt, end, n1 );
    plot("Poisson_dist.html",times, v1 );
    //open_plot("Poisson_dist.html");
}

static void testProcess(){
    double start = 0;
    double end = 1;
    int n = 200;
    double param = 30;
    vector<double> path(n,0);
    vector<double> martingale_path(n,0);
    path = generatePoisson_path(start, end, n, param);
    martingale_path = generatePoisson_Martingale(start, end, n, param);
    
    double dt = (end - start)/n;
    vector<double> times = linspace(dt, end, n );
    plot("Poisson_Path.html",times, path);
    //open_plot("Poisson_Path.html");
    
    plot("Poisson_Martingale.html",times, martingale_path);
    //open_plot("Poisson_Martingale.html");
}

static void testMerton(){
    int n = 500;
    double start = 0;
    double end = 5;
    double param = 1;
    double mu_j = 0;
    double sigma_j = 0.5;
    double sigma = 0.2;
    double s0 = 100;
    double r = 0.05;
    
    vector<double> merton;
    vector<double> s(n,0.0);
    merton = generate_Merton_paths(start, end, n, s0, param, mu_j, sigma_j, r, sigma, s);
    
    double dt = (end - start)/n;
    vector<double> times = linspace(dt, end, n );
    
    plot("Merton_Process.html",times, merton);
    open_plot("Merton_Process.html");
    
    plot("expMerton_Process.html",times, s);
    open_plot("expMerton_Process.html");
}

void testPoissonProcess(){
    TEST(testPoisson);
    TEST(testProcess);
    TEST(testMerton);
}
