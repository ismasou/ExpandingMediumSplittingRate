#include "ChebyshevTheta.h"
#include "Params.h"
#include <cmath>
#include <cstdlib>
#include <fmt/core.h>

// Numerical integrals using Chebyshev quadrature
double numericalInt1(double p, double q) {
    auto fct = [&p, &q](double theta){
        double k2 = p * p - 2.0 * p * q * std::cos(theta) + q * q;
        return WaveFct::CRate(k2);
    };
    double sum = thetaIntegrate(fct);
    return sum / (2.0 * M_PI);
}

double numericalInt2(double p, double q) {
    auto fct = [&p, &q](double theta){
        double cosT = std::cos(theta);
        double k2 = p * p - 2.0 * p * q * cosT + q * q;
        return WaveFct::CRate(k2) * cosT;
    };
    double sum = thetaIntegrate(fct);
    return sum / (2.0 * M_PI);
}

int main() {
    srand(42);
    const int Ntests = 10;
    const double tol = 1e-8;
    int passed = 0;

    fmt::print("Testing theta integrals: numerical vs analytical\n\n");

    for (int i = 0; i < Ntests; i++) {
        double p = 0.1 + 5.0 * rand() / RAND_MAX;
        double q = 0.1 + 5.0 * rand() / RAND_MAX;

        double num1 = numericalInt1(p, q);
        double ana1 = WaveFct::CRateAng(p * p, q * q, 1.0);
        double err1 = std::abs(num1 - ana1) / (std::abs(ana1) +
                                               1e-15); // relative error

        double num2 = numericalInt2(p, q);
        double ana2 = WaveFct::CRateAng2(p * p, q * q, 1.0);
        double err2 = std::abs(num2 - ana2) / (std::abs(ana2) +
                                               1e-15); // relative error

        bool ok1 = err1 < tol;
        bool ok2 = err2 < tol;

        std::string color1 = ok1 ? "\033[1;32m" : "\033[1;31m"; // Green if both pass, red otherwise
        std::string color2 = ok2 ? "\033[1;32m" : "\033[1;31m";

        fmt::print("Test {:2d}: p={:.3f} q={:.3f} \n", i + 1, p, q);
        fmt::print("{}  Int1: num={:.8f} ana={:.8f} err={:.2e} {}\n", color1, num1, ana1, err1,
                   ok1 ? "OK" : "FAIL");
        fmt::print("{}  Int2: num={:.8f} ana={:.8f} err={:.2e} {}\n\033[0m", color2, num2, ana2,
                   err2, ok2 ? "OK" : "FAIL");


        if (ok1 && ok2) { passed++; }
    }

    fmt::print("\nPassed: {}/{}\n", passed, Ntests);
    return (passed == Ntests) ? 0 : 1;
}
