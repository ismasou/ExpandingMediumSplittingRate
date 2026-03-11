#include "Polynomial.h"
#pragma once
#include <iostream>
#include <cmath>
#include "Params.h"


#define CHEBYSHEV
// #define TRAPEZOID
// #define GAUSSKONROD

#ifdef CHEBYSHEV
namespace Polynomial {
using namespace WaveFct;

LEGENDRE GaussPoints, Gausst, pPoints, IntegrationWeights, pWeights;
NT time;
const double tWeight = M_PI / NpM1;

inline double Cardinal(int i, double x) {

    return std::cos(NpM1 *std::acos(x)) * std::sin(Gausst[i]) /
           (NpM1 *std::sin(NpM1 *Gausst[i]) * (x - GaussPoints[i]));
}

inline double pTox(double p) {
    return 0.25 / M_PI *std::acos(p) - 1.0;
}

double pMin = 1e-2;
double pMax = 1e2;
double yMin = std::atan(pMin);
double yMax = std::atan(pMax);
inline double xTop(double x) {
    return std::tan(0.25 * M_PI * (x + 1.0));
}

inline double Jacobian(double p) {
    return 0.25 * (1.0 + p *p) * M_PI;
}

inline double timeFct(int i) {
    double tMin = 4.0 * dt;

    return tMin + (tmax - tMin) * std::pow(double(i) / (Nt - 1), 1.5);


}



void Setup() {

    for (int i = 0; i < Np; i++) {
        Gausst[i] = 0.5 * M_PI * (2 * i + 1) / NpM1;
        GaussPoints[i] = std::cos(Gausst[i]);
        IntegrationWeights[i] = M_PI / NpM1 *std::sin(Gausst[i]);
        pPoints[i] = xTop(GaussPoints[i]);
        pWeights[i] = Jacobian(pPoints[i]) * IntegrationWeights[i];
    }

    for (int i = 0; i < Nt; i++) {
        time[i] = timeFct(i);
    }
}

void TestIntegral() {
    // Let's integrate cos(x)^2 from 0 to pi
    auto f0 = [](double x) { return std::cos(x) * std::cos(x); };

    double res0 = 0.0;

    for (int i = 0; i < Np; i++) {
        res0 += tWeight *f0(Gausst[i]);
    }

    std::cout << "Integral of cos(x)^2 from 0 to pi: " << res0 << std::endl;
    std::cout << "True Value: " << M_PI / 2.0 << std::endl;
    std::cout << "Error: " << std::abs(res0 - M_PI / 2.0) << std::endl;

    // Let's integrate simple x^2/sqrt(1-x^2) from -1 to 1
    auto f4 = [](double x) { return x *x / std::sqrt(1.0 - x *x); };

    double res4 = 0.0;

    for (int i = 0; i < Np; i++) {
        res4 += IntegrationWeights[i] * f4(GaussPoints[i]);
    }

    std::cout << "Integral of x^2 from -1 to 1: " << res4 << std::endl;
    std::cout << "True Value: " << M_PI / 2.0 << std::endl;
    std::cout << "Error: " << std::abs(res4 - M_PI / 2.0) << std::endl;

    // Let's integrate simple x^2 from -1 to 1
    auto f1 = [](double x) { return x *x; };

    double res1 = 0.0;

    for (int i = 0; i < Np; i++) {
        res1 += IntegrationWeights[i] * f1(GaussPoints[i]);
    }

    std::cout << "Integral of x^2 from -1 to 1: " << res1 << std::endl;
    std::cout << "True Value: " << 2.0 / 3.0 << std::endl;
    std::cout << "Error: " << std::abs(res1 - 2.0 / 3.0) << std::endl;

    // Let's integrate x^2/(x^4+1) from 0 to infinity
    auto f2 = [](double x) { return x *x / (x *x *x *x + 1.0); };

    double TrueValue = 0.5 * M_PI / std::sqrt(2.0);

    // auto f2 = [](double x) { return 1.0 / (x * x + 1.0); };
    // double TrueValue = M_PI / 2.0;

    double res = 0.0;

    for (int i = 0; i < Np; i++) {
        res += pWeights[i] * f2(pPoints[i]);
    }

    std::cout << "Integral x^2/(x^4+1) from 0 to infinity: " << res << std::endl;
    std::cout << "True Value: " << TrueValue << std::endl;
    std::cout << "Error: " << std::abs(res - TrueValue) << std::endl;

    auto f3 = [](double x) { return 1.0 / (x *x + 1.0); };

    TrueValue = M_PI / 2.0;

    res = 0.0;

    for (int i = 0; i < Np; i++) {
        res += pWeights[i] * f3(pPoints[i]);
    }

    std::cout << "Integral 1/(x^2+1) from 0 to infinity: " << res << std::endl;
    std::cout << "True Value: " << TrueValue << std::endl;
    std::cout << "Error: " << std::abs(res - TrueValue) << std::endl;
}

}

#endif

#ifdef TRAPEZOID

#include <boost/math/quadrature/gauss_kronrod.hpp>
namespace Polynomial {

LEGENDRE GaussPoints, Gausst, pPoints, IntegrationWeights, pWeights;
const double tWeight = M_PI / NpM1;

inline double pTox(double p) {
    return 0.25 / M_PI *std::acos(p) - 1.0;
}

double pMin = 0.0;
double pMax = 2000.0;
inline double xTop(double x) {
    return std::tan(0.25 * M_PI * (x + 1.0));
    // return 0.5 * (x + 1.0) * (pMax - pMin) + pMin;
    // return (1.0 + x) / (1.0 - x);
}

inline double Jacobian(double p) {
    return 0.25 * (1.0 + p *p) * M_PI;
    // return 0.5 * (pMax - pMin);
    // return 0.5 * (1.0 + p) * (1.0 + p);
}



void Setup() {

    double xMin = 1e-3 - 1.0;
    double xMax = 1.0 - 1e-3;

    for (int i = 0; i < Np; i++) {
        double dx = (xMax - xMin) / NpM1;
        GaussPoints[i] = xMin + i *dx;
        IntegrationWeights[i] = ((i == 0 || i == NpM1) ? 0.5 : 1.0) * dx;
        pPoints[i] = xTop(GaussPoints[i]);
        pWeights[i] = Jacobian(pPoints[i]) * IntegrationWeights[i];
    }

}

void TestIntegral() {
    // Let's integrate cos(x)^2 from 0 to pi
    auto f0 = [](double x) { return std::cos(x) * std::cos(x); };

    double res0 = 0.0;

    for (int i = 0; i < Np; i++) {
        res0 += tWeight *f0(Gausst[i]);
    }

    std::cout << "Integral of cos(x)^2 from 0 to pi: " << res0 << std::endl;
    std::cout << "True Value: " << M_PI / 2.0 << std::endl;
    std::cout << "Error: " << std::abs(res0 - M_PI / 2.0) << std::endl;

    // Let's integrate simple x^2/sqrt(1-x^2) from -1 to 1
    auto f4 = [](double x) { return x *x / std::sqrt(1.0 - x *x); };

    double res4 = 0.0;

    for (int i = 0; i < Np; i++) {
        res4 += IntegrationWeights[i] * f4(GaussPoints[i]);
    }

    std::cout << "Integral of x^2 from -1 to 1: " << res4 << std::endl;
    std::cout << "True Value: " << M_PI / 2.0 << std::endl;
    std::cout << "Error: " << std::abs(res4 - M_PI / 2.0) << std::endl;

    // Let's integrate simple x^2 from -1 to 1
    auto f1 = [](double x) { return x *x; };

    double res1 = 0.0;

    for (int i = 0; i < Np; i++) {
        res1 += IntegrationWeights[i] * f1(GaussPoints[i]);
    }

    std::cout << "Integral of x^2 from -1 to 1: " << res1 << std::endl;
    std::cout << "True Value: " << 2.0 / 3.0 << std::endl;
    std::cout << "Error: " << std::abs(res1 - 2.0 / 3.0) << std::endl;

    // Let's integrate x^2/(x^4+1) from 0 to infinity
    auto f2 = [](double x) { return x *x / (x *x *x *x + 1.0); };

    double TrueValue = 0.5 * M_PI / std::sqrt(2.0);

    // auto f2 = [](double x) { return 1.0 / (x * x + 1.0); };
    // double TrueValue = M_PI / 2.0;

    double res = 0.0;

    for (int i = 0; i < Np; i++) {
        res += pWeights[i] * f2(pPoints[i]);
    }

    std::cout << "Integral x^2/(x^4+1) from 0 to infinity: " << res << std::endl;
    std::cout << "True Value: " << TrueValue << std::endl;
    std::cout << "Error: " << std::abs(res - TrueValue) << std::endl;

    auto f3 = [](double x) { return 1.0 / (x *x + 1.0); };

    TrueValue = M_PI / 2.0;

    res = 0.0;

    for (int i = 0; i < Np; i++) {
        res += pWeights[i] * f3(pPoints[i]);
    }

    std::cout << "Integral 1/(x^2+1) from 0 to infinity: " << res << std::endl;
    std::cout << "True Value: " << TrueValue << std::endl;
    std::cout << "Error: " << std::abs(res - TrueValue) << std::endl;
}

}

#endif

#ifdef GAUSSKONROD

#include <boost/math/quadrature/gauss_kronrod.hpp>
namespace Polynomial {

LEGENDRE GaussPoints, Gausst, pPoints, IntegrationWeights, pWeights;
const double tWeight = M_PI / NpM1;

inline double pTox(double p) {
    return 0.25 / M_PI *std::acos(p) - 1.0;
}

double pMin = 0.0;
double pMax = 2000.0;
inline double xTop(double x) {
    return std::tan(0.25 * M_PI * (x + 1.0));
    // return 0.5 * (x + 1.0) * (pMax - pMin) + pMin;
    // return (1.0 + x) / (1.0 - x);
}

inline double Jacobian(double p) {
    return 0.25 * (1.0 + p *p) * M_PI;
    // return 0.5 * (pMax - pMin);
    // return 0.5 * (1.0 + p) * (1.0 + p);
}



void Setup() {
    boost::math::quadrature::gauss_kronrod<double, Np> integrator;
    auto const &helper1 = integrator.abscissa();
    auto const &helper1W = integrator.weights();


    // Populate the positive values of the Stieltjes-quadrature
    for (int i = Np / 2; i < Np; i++)
    {
        GaussPoints[i]  = helper1[i - Np / 2];
        IntegrationWeights[i] = helper1W[i - Np / 2];
    }

    // Populate the negative values of the Stieltjes-quadrature
    int iMax = (Np % 2 == 0) ? Np / 2 - 1 : Np / 2;

    for (int i = 0; i < iMax; i++)
    {
        GaussPoints [Np / 2 - 1 - i]  = -GaussPoints[Np / 2 +
                                        1 + i];
        IntegrationWeights[Np / 2 - 1 - i]  =  IntegrationWeights[Np / 2 +
                                               1 + i];
    }

    for (int i = 0; i < Np; i++) {
        pPoints[i] = xTop(GaussPoints[i]);
        pWeights[i] = Jacobian(pPoints[i]) * IntegrationWeights[i];
    }

}

void TestIntegral() {
    // Let's integrate cos(x)^2 from 0 to pi
    auto f0 = [](double x) { return std::cos(x) * std::cos(x); };

    double res0 = 0.0;

    for (int i = 0; i < Np; i++) {
        res0 += tWeight *f0(Gausst[i]);
    }

    std::cout << "Integral of cos(x)^2 from 0 to pi: " << res0 << std::endl;
    std::cout << "True Value: " << M_PI / 2.0 << std::endl;
    std::cout << "Error: " << std::abs(res0 - M_PI / 2.0) << std::endl;

    // Let's integrate simple x^2/sqrt(1-x^2) from -1 to 1
    auto f4 = [](double x) { return x *x / std::sqrt(1.0 - x *x); };

    double res4 = 0.0;

    for (int i = 0; i < Np; i++) {
        res4 += IntegrationWeights[i] * f4(GaussPoints[i]);
    }

    std::cout << "Integral of x^2 from -1 to 1: " << res4 << std::endl;
    std::cout << "True Value: " << M_PI / 2.0 << std::endl;
    std::cout << "Error: " << std::abs(res4 - M_PI / 2.0) << std::endl;

    // Let's integrate simple x^2 from -1 to 1
    auto f1 = [](double x) { return x *x; };

    double res1 = 0.0;

    for (int i = 0; i < Np; i++) {
        res1 += IntegrationWeights[i] * f1(GaussPoints[i]);
    }

    std::cout << "Integral of x^2 from -1 to 1: " << res1 << std::endl;
    std::cout << "True Value: " << 2.0 / 3.0 << std::endl;
    std::cout << "Error: " << std::abs(res1 - 2.0 / 3.0) << std::endl;

    // Let's integrate x^2/(x^4+1) from 0 to infinity
    auto f2 = [](double x) { return x *x / (x *x *x *x + 1.0); };

    double TrueValue = 0.5 * M_PI / std::sqrt(2.0);

    // auto f2 = [](double x) { return 1.0 / (x * x + 1.0); };
    // double TrueValue = M_PI / 2.0;

    double res = 0.0;

    for (int i = 0; i < Np; i++) {
        res += pWeights[i] * f2(pPoints[i]);
    }

    std::cout << "Integral x^2/(x^4+1) from 0 to infinity: " << res << std::endl;
    std::cout << "True Value: " << TrueValue << std::endl;
    std::cout << "Error: " << std::abs(res - TrueValue) << std::endl;

    auto f3 = [](double x) { return 1.0 / (x *x + 1.0); };

    TrueValue = M_PI / 2.0;

    res = 0.0;

    for (int i = 0; i < Np; i++) {
        res += pWeights[i] * f3(pPoints[i]);
    }

    std::cout << "Integral 1/(x^2+1) from 0 to infinity: " << res << std::endl;
    std::cout << "True Value: " << TrueValue << std::endl;
    std::cout << "Error: " << std::abs(res - TrueValue) << std::endl;
}

}

#endif

// int main() {
//     Polynomial::Setup();
//     Polynomial::TestIntegral();
//     return 0;
// }
