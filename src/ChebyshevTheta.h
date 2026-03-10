#ifndef THETA_QUADRATURE_H
#define THETA_QUADRATURE_H

#include <cmath>
#include <array>

// Quadrature for theta integral from 0 to 2*pi
// Exploits symmetry: integrand depends on cos(theta) which is even
// So we integrate from 0 to pi and multiply weights by 2
//
// Switch between methods by defining THETA_CHEBYSHEV before including
// Default: Trapezoidal (optimal for periodic functions)

constexpr int Ntheta = 128;

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
    for (int k = 0; k < Ntheta; k++) {
        pts[k] = M_PI * k / (Ntheta - 1);
    }
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
// Weight: (pi/(N-1)) * 2 = 2*pi/(N-1) (factor of 2 from cos(theta) symmetry)
inline const double thetaWeight = 2.0 * M_PI / (Ntheta - 1);

template <typename Func>
inline double thetaIntegrate(Func f) {
    double sum = 0.0;
    for (int k = 1; k < Ntheta-1; k++) {
        sum += thetaWeight * f(thetaPoints[k]);
    }
    sum += 0.5 * thetaWeight * (f(thetaPoints[0]) + f(thetaPoints[Ntheta-1]));
    return sum;
}
#endif
