#include "Polynomial.h"
#include <memory>
#define NQUADRATURE 701
#include "DoubleIntegral.h"
#include <omp.h>
#include "Params.h"
#include <fmt/os.h>
double C_1, C_z, C_zB, m1, mz, mzB;

namespace WaveFct {
double T0, P, z, zB, g, mDSqr, lambda, t0, t1, dt, tmax, alpha;
std::string Header, shortName, longName;
std::unique_ptr<fmt::ostream> RateFile;
LEGENDRE PsiVals, CintVals;
double Rate;
std::unique_ptr<DoubleGauusKronrod> doubleIntegrator;

double Integrand(const double q, const double p,
                 const int ip, const int iq) {
    const double p2 = p * p;
    const double q2 = q * q;

    const double z2 = z * z;
    const double zB2 = zB * zB;

    const double Prefactor1 = C_1 * CRateAng(p2, q2, 1.0)
                              + C_z * CRateAng(p2, q2, z2)
                              + C_zB * CRateAng(p2, q2, zB2);
    const double Prefactor2 = C_1 * CRateAng2(p2, q2, 1.0)
                              + C_z * CRateAng2(p2, q2, z2)
                              + C_zB * CRateAng2(p2, q2, zB2);
    auto pPsi = Prefactor1 * PsiVals[ip]
                - Prefactor2 * (p / q) * PsiVals[iq];

    if (std::abs(pPsi) < 1e-15) {
        return 0.0;
    }

    auto res = q * pPsi;
    return res;
}

void Initialize() {
    for (int i = 0; i < Np; i++) {
        double p = Polynomial::pPoints[i];
        double p2 = p * p;
        PsiVals[i] = p2 / DeltaE(p2);
    }
}

void Output(int i, double t) {
    std::fstream file;
    file.open("./Output/Evolution/output-" + std::to_string(i) + ".dat",
              std::ios::out);
    file << "p Re(Psi) Im(Psi) Re(Cint) Im(Cint) RateIntegrand" << std::endl;

    for (int ip = 0; ip < Np; ip++) {
        double p = Polynomial::pPoints[ip];
        double p2 = p * p;
        double logExp = -DeltaE(p2) * t;
        auto RateIntegrand = p * real(PsiVals[ip] *
                                      std::complex<double>(std::cos(logExp), std::sin(logExp)));
        file << Polynomial::pPoints[ip] << " "
             << PsiVals[ip] << " "
             << CintVals[ip] << " "
             << RateIntegrand << std::endl;
    }

    file.close();
}

void ComputeRate(double dt, double t) {
    double ra = 0.0;

    int step = 1;

    #pragma omp parallel for reduction(+:ra)

    for (int ip = step; ip < Np; ip += step) {
        double p = Polynomial::pPoints[ip];
        double p2 = p * p;
        double logExp = DeltaE(p2) * t;
        double Integrandq = 0.0;

        #pragma omp parallel for reduction(+:Integrandq)

        for (int iq = 0; iq < Np; iq++) {
            double q = Polynomial::pPoints[iq];

            Integrandq += Integrand(q, p, ip, iq) * Polynomial::pWeights[iq];
        }

        ra += p * Polynomial::pWeights[ip] * Integrandq * (1.0 - std::cos(
                  logExp)) / DeltaE(p2);
    }

    Rate = ra / (4.0 * M_PI * M_PI);

}

void Evolve() {

    std::string fname = "Output/Opacity/opacity" + shortName;
    RateFile = std::make_unique<fmt::ostream>(fmt::output_file(fname + ".dat",
               fmt::file::WRONLY | fmt::file::TRUNC));
    RateFile->print("{}", Header);
    RateFile->print("t Rate \n");
    Rate = 0.0;
    double t = 0.0;
    double dt = 1e-1;
    int i = 0;


    RateFile->print("{} 0.0\n", t);

    while (t <= tmax) {
        t += dt;
        i++;
        ComputeRate(dt, t);
        RateFile->print("{} {}\n", t, Rate);
    }

    RateFile->close();
}


void Setup() {
    doubleIntegrator = std::make_unique<DoubleGauusKronrod>();
    Initialize();
    Evolve();

}



void Cleanup() {
}

}



int main(int argc, char *argv[]) {

    WaveFct::Setup(argc, argv);
    Polynomial::Setup();
    WaveFct::Setup();
    return 0;
}
