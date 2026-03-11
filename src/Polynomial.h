#pragma once

#include <array>
#include <cmath>
#include <complex>

#define Nt 32
// #define Np 512
const int Np = 512;
#define NpM1 (Np - 1)
typedef std::array<double, Np> LEGENDRE;
typedef std::array<std::complex<double>, Np> CLEGENDRE;

typedef std::array<std::complex<double>, Nt> CNT;
typedef std::array<double, Nt> NT;

typedef std::array<NT, Np> LEGENDRE2;
typedef std::array<CNT, Np> CLEGENDRE2;

namespace Polynomial {

extern LEGENDRE GaussPoints, Gausst, pPoints, IntegrationWeights, pWeights;
extern const double tWeight;
extern NT time;


inline double Cardinal(int i, double x);

inline double pTox(double p);

inline double xTop(double x);

inline double Jacobian(double p);
// int timeIndex(double t);

void Setup() ;

}
