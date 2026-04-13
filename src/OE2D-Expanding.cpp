#include "Polynomial.h"
#include <iostream>
#include <omp.h>
#include <fstream>
#include <fmt/core.h>
#include <fmt/os.h>
#define BJORKEN
#include "Params.h"
#define EULER
double C_1, C_z, C_zB, m1, mz, mzB;

namespace WaveFct {
std::fstream RateFile;
CLEGENDRE2 PsiVals, CintVals;
NT Rate = {}, RateL = {}, RateR = {};

double T0, P, z, zB, g, mDSqr, lambda, t0, t1, dt, tmax, alpha;
std::string Header, shortName, longName;

void Integrand(const double q, const double p,
               const int ip, const int iq, int itMin, double Deltat, CNT &Cint) {
    const double Jac = 1.0 / (2.0 * M_PI);
    const double p2 = p *p;
    const double q2 = q *q;

    const double z2 = z *z;
    const double zB2 = zB *zB;

    for (int it = itMin; it < Nt; it++) {
        double t = Polynomial::time[it];
        double t1 = t - Deltat;

        const double logExp = (p2 - q2) * Deltat;
        auto dExp = std::complex<double>(std::cos(logExp), std::sin(logExp));
        const double Prefactor1 = C_1 *CRateAng(p2, q2, 1.0, t1)
                                  + C_z *CRateAng(p2, q2, z2, t1)
                                  + C_zB *CRateAng(p2, q2, zB2, t1);
        const double Prefactor2 = C_1 *CRateAng2(p2, q2, 1.0, t1)
                                  + C_z *CRateAng2(p2, q2, z2, t1)
                                  + C_zB *CRateAng2(p2, q2, zB2, t1);
        auto pPsi = q *Prefactor1 *PsiVals[ip][it]
                    - Prefactor2 *p *dExp *PsiVals[iq][it];

        if (std::abs(pPsi) < 1e-15) {
            continue;
        }

        Cint[it] += pPsi *Jac *Polynomial::pWeights[iq];
        #ifdef DEBUG

        #pragma omp critical

        if (!std::isfinite(Cint[it].real()) || !std::isfinite(Cint[it].imag())) {
            std::cout << "Error: " << Cint[it] << std::endl;
            throw std::runtime_error("Error: Integrand is not finite");
        }

        #endif
    }
}

void Integrate(double Deltat, CLEGENDRE2 &CintVals, int itMin) {

    #pragma omp parallel for

    for (int ip = 0; ip < Np; ip++) {
        double p = Polynomial::pPoints[ip];
        CintVals[ip] = {};

        for (int iq = 0; iq < Np; iq++) {
            double q = Polynomial::pPoints[iq];

            Integrand(q, p, ip, iq, itMin, Deltat, CintVals[ip]);
        }
    }
}

void Initialize() {

    const int Niter = 1e4;
    const double dx = 1e-3;

    for (int ip = 0; ip < Np; ip++) {
        double p = Polynomial::pPoints[ip];
        double p2 = p *p;

        for (int it = 0; it < Nt; it++) {
            double tm = Polynomial::time[it];
            std::complex<double> F = 0.0;
            double x = 0.0;

            for (int ix = 0; ix < Niter; ix++) {
                double t = 1.0 / (x + dx);

                if (t < tm) {
                    break;
                }

                double t2 = t *t;
                F = (F / dx + t2) /
                    (1.0 / dx + std::complex<double>(0.0, DeltaE(p2, t) * t2));



                x += dx;
            }

            PsiVals[ip][it] = -p2 *F;
            #ifdef DEBUG

            #pragma omp critical

            if (!std::isfinite(PsiVals[ip][it].real())
                    || !std::isfinite(PsiVals[ip][it].imag())) {
                std::cout << "Error: " << PsiVals[ip][it] << std::endl;
                throw std::runtime_error("Error: Initial Psi is not finite");
            }

            #endif
        }
    }

    Integrate(0.0, CintVals, 0);
    PsiVals = CintVals;
}


void ComputeRate(double dt, double Deltat, int itMin) {

    for (int it = itMin; it < Nt; it++) {
        double t = Polynomial::time[it];
        double ra = 0.0;
        double DeltatR = Deltat + dt;

        if (DeltatR >= t) {
            DeltatR = t;
        }

        #pragma omp parallel for reduction(+:ra)

        for (int ip = 0; ip < Np; ip++) {
            double p = Polynomial::pPoints[ip];
            double p2 = p *p;
            double logExp = -DeltaE(p2, t, DeltatR);
            ra += p *Polynomial::pWeights[ip] * real(PsiVals[ip][it] *
                  std::complex<double>(std::cos(logExp), std::sin(logExp)));

            #ifdef DEBUG

            #pragma omp critical

            if (!std::isfinite(ra)) {
                std::cerr << "t: " << t << " Deltat: " << Deltat << " dt: " << dt << " ra:"
                          << ra << std::endl;
                throw std::runtime_error("Error: Rate is not finite");
            }

            #endif
        }

        RateR[it] = ra / (2.0 * M_PI);

        Rate[it] += 0.5 * dt * (RateL[it] + RateR[it]);

        RateL[it] = RateR[it];
    }

}

void Evolve() {

    // Strip last part of string

    std::string fname = "Output/Opacity/bjorken" + longName + ".dat";
    fmt::ostream RateFile = fmt::output_file(fname);

    RateFile.print("{}", Header);
    RateFile.print("t Rate\n");
    RateFile.print("0.0 0.0\n");
    double Deltat = 0.0;
    int i = 0;


    int itDone = 0;

    while (Deltat <= tmax + dt) {
        if (i % 100 == 0) {
            std::cout << "Deltat: " << Deltat << "\n";
            //     Integrate(t, CintVals);
            //     Output(i / 100, t);
        }

        #ifdef DEBUG

        if (Deltat >= Polynomial::time[itDone]) {
            std::cerr << "itDone: " << itDone << " Deltat: " << Deltat << " time[itDone]: "
                      <<
                      Polynomial::time[itDone] <<
                      std::endl;
            throw std::runtime_error("Error: Time index is not correct");
        }

        #endif

        ComputeRate(dt, Deltat, itDone);
        Deltat += dt;
        i++;

        for (int it = itDone; it < Nt; it++) {
            if (Deltat >= Polynomial::time[it]) {
                RateFile.print("{} {}\n", Polynomial::time[it], Rate[it]);
                // std::cout << "t: " << Polynomial::time[it] << " Rate: " << Rate[it] <<
                //           std::endl;
                // std::cout << "Deltat: " << Deltat << std::endl;
                itDone = it + 1;
            }
        }

    }

    RateFile.close();
}


void Setup() {
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
