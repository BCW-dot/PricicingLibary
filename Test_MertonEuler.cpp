/*
#include "matlib.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "Test_MertonEuler.h"

//#include <math.h>
using namespace std;

double variance_jump_process(double lam, double sigJ, double muJ){
    return sqrt(lam * sigJ * sigJ + lam * muJ * muJ);
}

std::vector<double> init_condition(std::vector<double> x, double strike){
    vector<double> init_cond(x.size());
    for (int i = 0; i < int(x.size()); ++i) {
        init_cond[i] = payoff(exp(x[i]), strike); //European call option
    }
    return init_cond;
}
 

// Left boundary-condition (vanilla call option)
double boundary_left(double t, double x){
  return 0.0;  // Specifically for a CALL option
}

// Right boundary-condition (vanilla call option)
double boundary_right(double t, double maturity, double r, double strike, double x){
    return (exp(x) - strike * exp(-r*(maturity - t)));
}

std::vector<double> calc_boundary(double t, double maturity, double r, double strike, int extraP, std::vector<double> x){
    std::vector<double> boundary_vector(x.size());
    for (int i = x.size() - extraP - 1; i < x.size(); ++i) {
        boundary_vector[i] = boundary_right(t, maturity, r, strike, x[i]);// exp(x[i]) - K * exp(-r * j * dt);
    }
    for (int i = 0; i <= extraP; ++i) {
        boundary_vector[i] = 0;
    }
    return boundary_vector;
}

double calculate_integral(const std::vector<double>& V, const std::vector<double>& nu, int extraP, int index, double dx) {
    // Ensure index is within bounds
    int V_size = V.size();
    if (index < extraP || index >= V_size - extraP) {
        throw std::out_of_range("Index out of bounds for integral calculation");
    }

    // Calculate the integral as a sum
    double integral = 0.0;
    for (int k = -extraP; k <= extraP; ++k) {
        integral += nu[k + extraP] * V[index + k];
    }

    return integral;
}

std::vector<double> calculate_jump_integral(const std::vector<double>& V, const std::vector<double>& nu, int extraP, double dt) {
    int N = V.size();
    std::vector<double> V_jump(N, 0.0);

    // Calculate the truncated jump integral
    for (int i = extraP + 1; i < N - extraP - 1; ++i) {
        double sum = 0.0;
        for (int j = -extraP; j <= extraP; ++j) {
            int idx = i + j;
            if (idx >= 0 && idx < N) {
                sum += nu[extraP + j] * V[idx];
            }
        }
        V_jump[i] = dt * sum;
    }

    return V_jump;
}


void model_param(){
    // Market and option parameters
    double r = 0.1; //risk free rate
    r = 0.1;
    double sig = 0.2; //volatility
    sig = 0.2;
    double S0 = 100; //start stock price
    S0 = 100;
    
    double K = 100; //strike
    double Texpir = 1; //maturity

    // Jump parameters
    double lam = 0.8; //Poisson parameter
    double muJ = 0; //mean jump size
    double sigJ = 0.5; //variance jump size

    // Discretization parameters
    int Nspace = 5; //"space Dimension"
    int Ntime = 3; //Time dimension
    double S_max = 3 * K;
    double S_min = K / 3;
    double x_max = log(S_max);
    double x_min = log(S_min);

    double dev_X = variance_jump_process(lam, sigJ, muJ); //variance jump process var(sum_{0}^{N_t}Y_{i}

    double dx = (x_max - x_min) / (Nspace - 1);
    //std::cout << "dx has value: " << dx << endl;
    int extraP = int(floor(3 * dev_X / dx));
    vector<double> x = linspace2(x_min - extraP * dx, x_max + extraP * dx, Nspace + 2 * extraP);
    //print(x);
    double dt = Texpir / (Ntime - 1);
    //std::cout << "dt has value: " << dt << endl;
    
    //cout << "dx has the value: " << dx << endl;
    //cout << "Under this discretization there are " << extraP << " extra points" << endl;

    vector<double> Payoff = init_condition(x, K);
    vector<double> V_old = Payoff;
    print(V_old);
 
    
    //make sure to call this after init_cond 
    V_old = calc_boundary((Ntime-3)*dt, Texpir, r, K, extraP, x);
    //print(V_old);
 
    
    vector<double> V_new(x.size());
    //for (int j = 0; j < Ntime; ++j) {
        for (int i = x.size() - extraP - 1; i < x.size(); ++i) {
            V_new[i] = boundary_right((Ntime-3)*dt, Texpir, r, K, x[i]);// exp(x[i]) - K * exp(-r * j * dt);
        }
        for (int i = 0; i <= extraP; ++i) {
            V_new[i] = 0;
        }
    //}
    //print(V_new);
    

    vector<double> cdf=linspace(-(extraP + 1 + 0.5) * dx, (extraP + 1 + 0.5) * dx, 2 * (extraP + 2));
    for (int i = 0; i < cdf.size(); ++i) {
        cdf[i] = normal_cdf(-(extraP + 1 + 0.5) * dx + i * dx, muJ, sigJ);
    }
    vector<double> nu;//(extraP * 2 + 1);
    for (int i = 0; i < cdf.size()-1; ++i) {
        nu.push_back(lam * (cdf[i + 1] - cdf[i]));
    }
    //cout<<"nu: ";
    //print(nu);
    

    double lam_appr;
    lam_appr = accumulate(nu.begin(), nu.end(), 0.0);
    //cout << "Truncated jump activity: " << lam_appr << endl;

    //double m = lam * (exp(muJ + 0.5 * sigJ * sigJ) - 1);
    //cout << "True value: " << m << endl;

    double m_int = 0;
    for (double z = -(extraP + 1.5) * dx; z <= (extraP + 1.5) * dx; z += dx) {
        m_int += lam * (exp(z) - 1) * exp(-0.5 * pow((z - muJ) / sigJ, 2)) / (sigJ * sqrt(2 * M_PI)) * dx;
    }
    //cout << "Truncated value, using quad: " << m_int << endl;

    double m_appr = 0;
    for (int i = 0; i < nu.size(); ++i) {
        m_appr += (exp(-(extraP + 1) * dx + i * dx) - 1) * nu[i];
    }
    //cout << "Approximation value: " << m_appr << endl;
    //if (index < extraP || index >= V_size - extraP) {
    for(int j = extraP + 1; j < V_old.size() - extraP ; j++){
        double jump_contribution = calculate_integral(V_old, nu, extraP, j, dx);
        cout<< jump_contribution <<", ";
    }
    cout << endl;
    
    std::vector<double> jumps = calculate_jump_integral(V_old, nu, extraP, (Ntime-1)*dt);
    print(jumps);
    
    std::vector<double> jumps2 = calculate_jump_integral(V_old, nu, extraP, dt);
    print(jumps2);

}


*/

#include "Test_MertonEuler.h"
#include "matlib.h"
using namespace std;

/*
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>

class MertonPIDE {
private:
    double r, sig, S0, K, Texpir, lam, muJ, sigJ;
    int Nspace, Ntime, extraP;
    double S_max, S_min, x_max, x_min, dx, dt;
    std::vector<double> x, T;
    
    std::vector<double> init_condition(std::vector<double> x){
        vector<double> init_cond(x.size());
        for (int i = 0; i < int(x.size()); ++i) {
            init_cond[i] = payoff(exp(x[i]), K); //European call option
        }
        return init_cond;
    }
    
    // Left boundary-condition (vanilla call option)
    double boundary_left(double t, double x){
      return 0.0;  // Specifically for a CALL option
    }

    // Right boundary-condition (vanilla call option)
    double boundary_right(double t, double maturity, double r, double strike, double x){
        return (exp(x) - strike * exp(-r*(maturity - t)));
    }

    std::vector<double> calc_boundary(double t, double maturity, double r, double strike, int extraP, std::vector<double> x){
        std::vector<double> boundary_vector(x.size());
        for (int i = x.size() - extraP - 1; i < x.size(); ++i) {
            boundary_vector[i] = boundary_right(t, maturity, r, strike, x[i]);// exp(x[i]) - K * exp(-r * j * dt);
        }
        for (int i = 0; i <= extraP; ++i) {
            boundary_vector[i] = 0;
        }
        return boundary_vector;
    }

public:
    MertonPIDE(double r, double sig, double S0, double K, double Texpir,
               double lam, double muJ, double sigJ, int Nspace, int Ntime)
        : r(r), sig(sig), S0(S0), K(K), Texpir(Texpir), lam(lam), muJ(muJ), sigJ(sigJ),
          Nspace(Nspace), Ntime(Ntime) {
        
        //double X0 = std::log(S0);
        S_max = 3 * K;
        S_min = K / 3;
        x_max = std::log(S_max);
        x_min = std::log(S_min);

        double dev_X = std::sqrt(lam * sigJ * sigJ + lam * muJ * muJ);
        dx = (x_max - x_min) / (Nspace - 1);
        extraP = static_cast<int>(std::floor(3 * dev_X / dx));

        x.resize(Nspace + 2 * extraP);
        for (int i = 0; i < Nspace + 2 * extraP; ++i) {
            x[i] = x_min - extraP * dx + i * dx;
        }

        dt = Texpir / (Ntime - 1);
        T.resize(Ntime);
        for (int i = 0; i < Ntime; ++i) {
            T[i] = i * dt;
        }
    }

    double solve() {
        //cout <<"extraP: " << extraP <<endl;
        //cout<< "dx, dt: " << dx << ", " << dt << endl;
        //cout << "x: ";
        //print(x);
        //cout << "dt, T: " << dt << endl;
        //print(T);
        std::vector<double> V(Nspace + 2 * extraP);
        std::vector<double> V_new(Nspace + 2 * extraP);

        // Initialize payoff
        for (int i = 0; i < Nspace + 2 * extraP; ++i) {
            V[i] = std::max(std::exp(x[i]) - K, 0.0);
        }
        V = init_condition(x);
        //print(V);

        // Calculate nu
        std::vector<double> cdf=linspace(-(extraP + 1 + 0.5) * dx, (extraP + 1 + 0.5) * dx, 2 * (extraP + 2));
        for (int i = 0; i < cdf.size(); ++i) {
            cdf[i] = normal_cdf(-(extraP + 1 + 0.5) * dx + i * dx, muJ, sigJ);
        }
        std::vector<double> nu;//(extraP * 2 + 1);
        for (int i = 0; i < cdf.size()-1; ++i) {
            nu.push_back(lam * (cdf[i + 1] - cdf[i]));
        }
        //print(nu);
        double lam_appr = std::accumulate(nu.begin(), nu.end(), 0.0);
        cout << "lam_approx: " << lam_appr << endl;
        // Calculate m_appr
        double m_appr = 0;
        for (int i = 0; i < nu.size(); ++i) {
            m_appr += (exp(-(extraP + 1) * dx + i * dx) - 1) * nu[i];
        }
        //cout << "m_approx: " << m_appr << endl;
        
        // Construct tri-diagonal coefficients
        double sig2 = sig * sig;
        double dxx = dx * dx;
        double a = (dt / 2) * ((r - m_appr - 0.5 * sig2) / dx - sig2 / dxx);
        double b = 1 + dt * (sig2 / dxx + r + lam_appr);
        double c = -(dt / 2) * ((r - m_appr - 0.5 * sig2) / dx + sig2 / dxx);
         
        //cout << "sig2: " << sig2 << ", dxx: " << dxx << ", a: " << a << ", b: " << b << ", c: " << c << endl;
        //std::cout << "a: " << a << ", b: " << b << ", c: " << c << std::endl;

        // Backward iteration
        //cout << "V: " << endl;
        print(V);
        print(nu);
        cout << "running total number of time loops: " << Ntime-2 + 1<< endl;
        for (int n = Ntime-2; n >= 0; --n) {
            V_new = calc_boundary(n*dt, Texpir, r, K, extraP, x);
            
            
            // Fill J with nu values, shifted for tiem step n, this is wrong i think
            std::vector<double> J(Nspace + 2 * extraP,0.0);
            for (int j = Ntime-2-n; j < nu.size(); j++) {
                J[j] = nu[j-(Ntime-2-n)];
            }
                
            cout << "J_" << n << endl;
            print(J);
            
            // Interior points
            for (int i = extraP + 2; i < Nspace + extraP - 2; ++i) {
                double jump_integral = 0;
                for (int k = 0; k<J.size(); k++){
                    jump_integral += J[k]*V[k];
                }
                V_new[i] = (a * V[i-1] + b * V[i] + c * V[i+1] + dt * jump_integral);
            }
             
            
            for (int i = extraP + 2; i < Nspace + extraP - 2; ++i) {
                double jump_integral = 0.0;
                for (int k = 0; k < nu.size(); ++k) {
                    int j = i - extraP - 1 + k;
                    if (j >= 0 && j < V.size()) {
                        jump_integral += nu[k] * V[j];
                    }
                }
                V_new[i] = a * V[i-1] + b * V[i] + c * V[i+1] + dt * jump_integral;
            }
            V = V_new;
        }
        

        // Interpolate to find option price at S0
        double X0 = std::log(S0);
        int index = std::lower_bound(x.begin(), x.end(), X0) - x.begin();
        double alpha = (X0 - x[index-1]) / (x[index] - x[index-1]);
        print(V_new);
        return V[index-1] * (1 - alpha) + V[index] * alpha;
    }
};
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>

class MertonPIDE {
private:
    double r, sig, S0, K, Texpir, lam, muJ, sigJ;
    int Nspace, Ntime, extraP;
    double S_max, S_min, x_max, x_min, dx, dt;
    std::vector<double> x, T, nu;

    double normalCDF(double x) {
        return 0.5 * (1 + std::erf(x / std::sqrt(2)));
    }

    double normalPDF(double x) {
        return std::exp(-0.5 * x * x) / std::sqrt(2 * M_PI);
    }
    
    std::vector<double> init_condition(const std::vector<double>& x) {
        std::vector<double> init_cond(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            init_cond[i] = std::max(std::exp(x[i]) - K, 0.0); // European call option
        }
        return init_cond;
    }
    
    // Left boundary-condition (vanilla call option)
    double boundary_left(double t, double x) {
        return 0.0;  // Specifically for a CALL option
    }

    // Right boundary-condition (vanilla call option)
    double boundary_right(double t, double maturity, double r, double strike, double x) {
        return (std::exp(x) - strike * std::exp(-r * (maturity - t)));
    }

    std::vector<double> calc_boundary(double t, double maturity, double r, double strike, int extraP, const std::vector<double>& x) {
        std::vector<double> boundary_vector(x.size());
        for (size_t i = x.size() - extraP - 1; i < x.size(); ++i) {
            boundary_vector[i] = boundary_right(t, maturity, r, strike, x[i]);
        }
        for (int i = 0; i <= extraP; ++i) {
            boundary_vector[i] = boundary_left(t, x[i]);
        }
        return boundary_vector;
    }

public:
    MertonPIDE(double r, double sig, double S0, double K, double Texpir,
               double lam, double muJ, double sigJ, int Nspace, int Ntime)
        : r(r), sig(sig), S0(S0), K(K), Texpir(Texpir), lam(lam), muJ(muJ), sigJ(sigJ),
          Nspace(Nspace), Ntime(Ntime) {
        
        S_max = 3 * K;
        S_min = K / 3;
        x_max = std::log(S_max);
        x_min = std::log(S_min);

        double dev_X = std::sqrt(lam * sigJ * sigJ + lam * muJ * muJ);
        dx = (x_max - x_min) / (Nspace - 1);
        extraP = static_cast<int>(std::floor(3 * dev_X / dx));

        x.resize(Nspace + 2 * extraP);
        for (int i = 0; i < Nspace + 2 * extraP; ++i) {
            x[i] = x_min - extraP * dx + i * dx;
        }

        dt = Texpir / (Ntime - 1);
        T.resize(Ntime);
        for (int i = 0; i < Ntime; ++i) {
            T[i] = i * dt;
        }

        // Calculate nu
        std::vector<double> y(2 * extraP + 1);
        for (int i = 0; i < static_cast<int>(y.size()); ++i) {
            y[i] = -extraP * dx + i * dx;
        }

        nu.resize(y.size());
        for (size_t i = 0; i < nu.size(); ++i) {
            nu[i] = lam * dt * normalPDF((y[i] - muJ) / sigJ) / sigJ;
        }

        // Normalize nu to ensure the total probability is lam * dt
        double sum_nu = std::accumulate(nu.begin(), nu.end(), 0.0);
        for (size_t i = 0; i < nu.size(); ++i) {
            nu[i] *= lam * dt / sum_nu;
        }
    }

    double solve() {
        std::vector<double> V(Nspace + 2 * extraP);
        std::vector<double> V_new(Nspace + 2 * extraP);

        // Initialize payoff
        V = init_condition(x);

        // Construct tri-diagonal coefficients
        double sig2 = sig * sig;
        double dxx = dx * dx;
        double a = (dt / 2) * ((r - 0.5 * sig2) / dx - sig2 / dxx);
        double b = 1 + dt * (sig2 / dxx + r + lam);
        double c = -(dt / 2) * ((r - 0.5 * sig2) / dx + sig2 / dxx);

        // Backward iteration
        for (int n = Ntime - 2; n >= 0; --n) {
            V_new = calc_boundary(n * dt, Texpir, r, K, extraP, x);

            // Interior points
            for (int i = extraP + 2; i < Nspace + extraP - 2; ++i) {
                double jump_integral = 0.0;
                for (size_t k = 0; k < nu.size(); ++k) {
                    int j = i - extraP - 1 + k;
                    if (j >= 0 && j < static_cast<int>(V.size())) {
                        jump_integral += nu[k] * V[j];
                    }
                }
                V_new[i] = a * V[i-1] + b * V[i] + c * V[i+1] + dt * jump_integral;
            }
            std::swap(V, V_new);
        }

        // Interpolate to find option price at S0
        double X0 = std::log(S0);
        auto it = std::lower_bound(x.begin(), x.end(), X0);
        int index = std::distance(x.begin(), it);
        double alpha = (X0 - x[index-1]) / (x[index] - x[index-1]);
        return V[index-1] * (1 - alpha) + V[index] * alpha;
    }
};

static void test_Merton_PIDE() {
    double r = 0.1;
    double sig = 0.2;
    double S0 = 100;
    double K = 100;
    double Texpir = 1;
    double lam = 0.8;
    double muJ = 0;
    double sigJ = 0.5;
    int Nspace = 5;
    int Ntime = 3;

    MertonPIDE merton(r, sig, S0, K, Texpir, lam, muJ, sigJ, Nspace, Ntime);
    //merton.printGrid();
    double option_price = merton.solve();

    std::cout << "Option price: " << option_price << std::endl;
}

void testMertonEuler(){
    //TEST(model_param);
    TEST(test_Merton_PIDE) ;
}

