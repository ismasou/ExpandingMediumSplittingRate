#ifndef THETA_QUADRATURE_H
#define THETA_QUADRATURE_H

#include <cmath>
#include <array>

namespace ThetaQuad {

// Quadrature for theta integral from 0 to 2*pi
// Exploits symmetry: integrand depends on cos(theta) which is even
// So we integrate from 0 to pi and multiply weights by 2
//
// Switch between methods by defining THETA_CHEBYSHEV before including
// Default: Trapezoidal (optimal for periodic functions)

constexpr int Ntheta = 32;

#ifdef THETA_CHEBYSHEV
// Chebyshev-Gauss nodes on [0, pi]: theta_k = pi*(2k+1)/(2N)
inline std::array<double, Ntheta> makeThetaPoints() {
    std::array<double, Ntheta> pts{};

    for (int k = 0; k < Ntheta; k++) {
        pts[k] = M_PI * (2.0 * k + 1.0) / (2.0 * Ntheta);
    }

    return pts;
}

#else
// Trapezoidal rule on [0, pi]: theta_k = pi*k/(N-1)
inline std::array<double, Ntheta> makeThetaPoints() {
    std::array<double, Ntheta> pts{};
    for (int k = 0; k < Ntheta; k++) pts[k] = M_PI * k / (Ntheta - 1);
    return pts;
}
#endif

// Weight: pi/(N-1) * 2 = 2*pi/(N-1) (factor of 2 from cos(theta) symmetry)
inline const std::array<double, Ntheta> makeThetaWeights() {
    std::array<double, Ntheta> w{};
    double h = 2.0 * M_PI / (Ntheta - 1);

    for (int k = 0; k < Ntheta; k++) {
        w[k] = h;
    }

    w[0] *= 0.5;
    w[Ntheta - 1] *= 0.5;
    return w;
}

inline const std::array<double, Ntheta> thetaPoints = makeThetaPoints();
// Precompute cos(theta) to avoid calling std::cos inside hot loops
inline const std::array<double, Ntheta> cosThetaPoints = []{
    std::array<double, Ntheta> a{};
    for (int k = 0; k < Ntheta; ++k) a[k] = std::cos(thetaPoints[k]);
    return a;
}();
inline const std::array<double, Ntheta> sinThetaPoints = []{
    std::array<double, Ntheta> a{};
    for (int k = 0; k < Ntheta; ++k) a[k] = std::sin(thetaPoints[k]);
    return a;
}();
// Weight: (pi/(N-1)) * 2 = 2*pi/(N-1) (factor of 2 from cos(theta) symmetry)
inline const double thetaWeight = 2.0 * M_PI / (Ntheta - 1);

template <typename Func>
inline double Integrate(Func f) {
    double sum = 0.0;

    for (int k = 1; k < Ntheta - 1; k++) {
        sum += thetaWeight * f(thetaPoints[k]);
    }

    sum += 0.5 * thetaWeight * (f(thetaPoints[0]) + f(thetaPoints[Ntheta - 1]));
    return sum;
}

template <size_t n, typename Func>
std::array<double, n> Integrate(Func f) {
    std::array<double, n> sum = {};

    for (int k = 1; k < Ntheta - 1; k++) {
        auto res = f(thetaPoints[k], cosThetaPoints[k], sinThetaPoints[k]);

        for (size_t i = 0; i < n; i++) {
            sum[i] += thetaWeight * res[i];
        }
    }

    auto res1 = f(thetaPoints[0], cosThetaPoints[0], sinThetaPoints[0]);
    auto res2 = f(thetaPoints[Ntheta - 1], cosThetaPoints[Ntheta - 1], sinThetaPoints[Ntheta - 1]);

    for (size_t i = 0; i < n; i++) {
        sum[i] += 0.5 * thetaWeight * (res1[i] + res2[i]);
    }

    return sum;
}
}

#endif
