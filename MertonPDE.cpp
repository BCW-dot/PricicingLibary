#include "MertonPDE.h"
#include "matlib.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

//#include <math.h>
using namespace std;


MertonPDE::MertonPDE(std::shared_ptr<PathIndependentOption> option, std::shared_ptr<BlackScholesModel> model,
                     double poisson_parameter, double mean_jump_size, double variance_jump_size, double m_hat, double lamda_hat, double nu) :
    option(option),
    model(model),
    poisson_parameter(poisson_parameter),
    mean_jump_size(mean_jump_size),
    variance_jump_size(variance_jump_size),
    m_hat(m_hat),
    lamda_hat(lamda_hat),
    nu(nu)
{}

// Diffusion coefficient
double MertonPDE::diff_coeff(double t, double x) const {
  double vol = model->getVolatility();
  return 0.5*vol*vol;  // \frac{1}{2} \sigma^2
}

// Convection coefficient
double MertonPDE::conv_coeff(double t, double x) const {
  return model->getRiskFreeRate() - 0.5*model->getVolatility()*model->getVolatility() - m_hat;  // r-0.5sig^2 - m_hat
}

// Zero-term coefficient
double MertonPDE::zero_coeff(double t, double x) const {
  return -(model->getRiskFreeRate() + lamda_hat);  // -r
}

// Source coefficient
double MertonPDE::source_coeff(double t, double x) const {
  return 0.0;
}

/*
double MertonPDE::integral_item(double t, double x) const{
    double sum =0;
    for(int i =0; i<n; i++){
        sum += nu * 
    }
}
*/
// Left boundary-condition (vanilla call option)
double MertonPDE::boundary_left(double t, double x) const {
  return 0.0;  // Specifically for a CALL option
}

//This needs to be modified for MertonJumps
// Right boundary-condition (vanilla call option)
double MertonPDE::boundary_right(double t, double x) const {
  // This is via Put-Call Parity and works for a call option
  return exp(x) - option->getStrike() * exp(-model->getRiskFreeRate()*(option->getMaturity() - t));
}

// Initial condition (vanilla call option)
double MertonPDE::init_cond(double x) const {
  return option->payoff(x);
}
/*
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

double normal_cdf(double x, double mu = 0, double sigma = 1) {
    return 0.5 * erfc(-(x - mu) / (sigma * sqrt(2)));
}

double payoff(double S, double K) {
    return max(S - K, 0.0);
}
 */

static void testPDIE() {
    // Market and option parameters
    double r = 0.1;
    double sig = 0.2;
    double S0 = 100;
    
    double K = 100;
    double Texpir = 1;

    // Jump parameters
    double lam = 0.8;
    double muJ = 0;
    double sigJ = 0.5;

    // Discretization parameters
    int Nspace = 5;
    int Ntime = 3;
    double S_max = 3 * K;
    double S_min = K / 3;
    double x_max = log(S_max);
    double x_min = log(S_min);

    double dev_X = sqrt(lam * sigJ * sigJ + lam * muJ * muJ);

    double dx = (x_max - x_min) / (Nspace - 1);
    //std::cout << "dx has value: " << dx << endl;
    int extraP = int(floor(3 * dev_X / dx));
    vector<double> x = linspace2(x_min - extraP * dx, x_max + extraP * dx, Nspace + 2 * extraP);
    //print(x);
    double dt = Texpir / (Ntime - 1);
    //std::cout << "dt has value: " << dt << endl;
    
    //cout << "dx has the value: " << dx << endl;
    //cout << "Under this discretization there are " << extraP << " extra points" << endl;

    vector<double> Payoff(x.size());
    for (int i = 0; i < x.size(); ++i) {
        Payoff[i] = max(exp(x[i]) - K, 0.0);
    }

    vector<double> V = Payoff;
    vector<double> V_new(x.size());

    for (int j = 0; j < Ntime; ++j) {
        for (int i = x.size() - extraP - 1; i < x.size(); ++i) {
            V[i] = exp(x[i]) - K * exp(-r * j * dt);
        }
        for (int i = 0; i <= extraP; ++i) {
            V[i] = 0;
        }
    }
    

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
    

    double lam_appr = accumulate(nu.begin(), nu.end(), 0.0);
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

    double sig2 = sig * sig;
    double dxx = dx * dx;
    double a = (dt / 2) * ((r - m_appr - 0.5 * sig2) / dx - sig2 / dxx);
    double b = 1 + dt * (sig2 / dxx + r + lam_appr);
    double c = -(dt / 2) * ((r - m_appr - 0.5 * sig2) / dx + sig2 / dxx);
    
    //cout<< "sig2: " << sig2 << ", dxx: " << dxx << ", a: " << a << ", b: " << b << ", c: ";

    for (int i = Ntime - 2; i >= 0; --i) {
            for (int j = extraP + 1; j < x.size() - extraP - 1; ++j) {
                double V_jump = V[j];
                for (int k = 0; k < nu.size(); ++k) {
                    V_jump += dt * nu[k] * V[j + k - extraP - 1];
                }
                V_new[j] = a * V[j - 1] + b * V[j] + c * V[j + 1] - V_jump;
            }

            for (int j = extraP + 1; j < x.size() - extraP - 1; ++j) {
                V[j] = V_new[j];
            }
        }
    double option_price = 0;
    for (int i = 0; i < x.size(); ++i) {
        if (x[i] >= log(S0) && x[i - 1] <= log(S0)) {
            option_price = V[i - 1] + (log(S0) - x[i - 1]) * (V[i] - V[i - 1]) / (x[i] - x[i - 1]);
            break;
        }
    }

    cout << "Option price: " << option_price << endl;
}

void testMertonPDE(){
    TEST( testPDIE );
}

