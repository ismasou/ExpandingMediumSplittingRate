#ifndef CHEBYSHEV_THETA_H
#define CHEBYSHEV_THETA_H

#include <cmath>
#include <array>

// Chebyshev quadrature for theta integral from 0 to 2*pi
// Uses Chebyshev-Gauss nodes: theta_k = pi*(2k+1)/N, weight = 2*pi/N
constexpr int Ntheta = 32;

inline std::array<double, Ntheta> makeThetaPoints() {
    std::array<double, Ntheta> pts{};
    for (int k = 0; k < Ntheta; k++) {
        pts[k] = M_PI * (2.0 * k + 1.0) / Ntheta;
    }
    return pts;
}

inline std::array<double, Ntheta> makeThetaWeights() {
    std::array<double, Ntheta> wts{};
    double w = 2.0 * M_PI / Ntheta;
    for (int k = 0; k < Ntheta; k++) {
        wts[k] = w;
    }
    return wts;
}

inline const std::array<double, Ntheta> thetaPoints = makeThetaPoints();
inline const std::array<double, Ntheta> thetaWeights = makeThetaWeights();

#endif
