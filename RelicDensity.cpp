#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_odeiv2.h>

const double c = 3.0e5; // Speed of light in km/s
const double Mp = 1.22e19; // Planck mass in GeV
const double mNu = 0.0; // ν mass
const double me = 0.0005; // electron mass
const double mMu = 0.105; // μ mass
const double mTau = 1.800; // τ mass
const double DensityDM = 0.1129;
const double GeVtocm2 = pow(1.0 / 5.06e13, 2); // cm^-2
const double GeVtog = 1.0 / 1.78e-24; // g
const double cf = GeVtocm2 * GeVtog;
const double DensityFactor = 1.5e8; // Product of s0 h^2/Subscript[ρ, c]



// Cross section function
double sigmaVtoZZApprox(double g, double m) {
    const double Pi = 3.14159265358979323846;
    return pow(g, 4) / (16.0 * Pi * pow(m, 2));
}

// Degree of freedom parameter
const double a = 10.2;
const double b = 2.349;
const double c_value = 0.252;

double gsToHalf(double T) {
    return a / (1.0 + exp(-b * (T - c_value)));
}

// Functions
double Yeq(double x) {
    return 0.145 * pow(x, 1.5) * exp(-x);
}


struct ODEParameters {
    double mass;
    double gprime;
};

int BoltzmannP(double x, const double y[], double f[], void* params) {
    ODEParameters* parameters = static_cast<ODEParameters*>(params);

    double mass = parameters -> mass;
    double gprime = parameters-> gprime;

    double Yeq_val = Yeq(x);
    double sigmaV = Mp * sqrt(M_PI / 45.0) * gsToHalf(mass / x) * sigmaVtoZZApprox(gprime, mass) * mass;

    f[0] = -sigmaV * (y[0] * y[0] - Yeq_val * Yeq_val) / (x * x);

    return GSL_SUCCESS;
}

int main() {
    const double x_start = 1.0;
    const double x_end = 100.0;

    std::vector<double> DMmassArray;

/*
    for (double i = 0.1; i <= 100.0; i += 0.01) {
        DMmassArray.push_back(i);
    }

    std::vector<double> gArray;
    for (double g = 0.001; g <= 0.3; g += 0.001) {
        gArray.push_back(g);
    }

    std::vector<std::vector<double>> RelicArray;
    std::vector<std::vector<double>> GetParams;

    gsl_odeiv2_system system = {BoltzmannP, nullptr, 1, nullptr};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    for (double gprime : gArray) {
        std::cout << "g_chi = " << gprime << std::endl;

        for (double mass : DMmassArray) {
            // Set initial conditions
            double x = x_start;
            double y[1] = {Yeq(x_start)};

            // Set parameters for the ODE
            ODEParameters parameters = {mass, gprime};
            system.params = &parameters;

            // Integrate the ODE
            while (x < x_end) {
                double x_next = x + 1.0; // Adjust as needed
                int status = gsl_odeiv2_driver_apply(driver, &x, x_next, y);

                if (status != GSL_SUCCESS) {
                    std::cerr << "Integration failed" << std::endl;
                    break;
                }

                // Save or process the results as needed
            }

            // Save Arrays
            // Implement saving of parameters and relic density arrays here
        }
    }

    gsl_odeiv2_driver_free(driver);

*/
    return 0;
}
