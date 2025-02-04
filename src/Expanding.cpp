#include "Polynomial.h"
#include <omp.h>
#include <fmt/os.h>
#define BJORKEN
#include "Params.h"
#define EULER
double C_1, C_z, C_zB, m1, mz, mzB;
std::string process = "G_GG";

namespace WaveFct {
CLEGENDRE2 PsiVals, CintVals;
NT Rate = {}, RateL = {}, RateR = {};

double T0, P, z, zB, g, mDSqr, lambda, t0, t1, dt, tmax, alpha;
std::string Header, shortName, longName;

/**
 * @brief This function calculates the integrand for a given set of parameters.
 *
 * @param q The q parameter for the integrand.
 * @param p The p parameter for the integrand.
 * @param ip The index of p in the PsiVals array.
 * @param iq The index of q in the PsiVals array.
 * @param itMin The minimum time index for the calculation.
 * @param Deltat The time step for the calculation.
 * @param Cint The array to store the calculated integrand values.
 *
 * The function calculates the integrand for each time step from itMin to Nt.
 * It uses the Polynomial::time array for the time values and the Polynomial::pWeights array for the weights.
 * The calculated integrand is added to the corresponding element of the Cint array.
 * If the absolute value of the calculated integrand is less than 1e-15, it is ignored.
 * In debug mode, the function checks if the calculated integrand is finite and throws a runtime error if it is not.
 */
void Integrand(const double q, const double p,
               const int ip, const int iq, int itMin, double Deltat, CNT &Cint) {
    const double Jac = 1.0 / (2.0 * M_PI);
    const double p2 = p * p;
    const double q2 = q * q;

    const double z2 = z * z;
    const double zB2 = zB * zB;

    for (int it = itMin; it < Nt; it++) {
        double t = Polynomial::time[it];
        double t1 = t - Deltat;

        const double logExp = (p2 - q2) * Deltat;
        auto dExp = std::complex<double>(std::cos(logExp), std::sin(logExp));
        const double Prefactor1 = C_1 * CRateAng(p2, q2, 1.0, t1)
                                  + C_z * CRateAng(p2, q2, z2, t1)
                                  + C_zB * CRateAng(p2, q2, zB2, t1);
        const double Prefactor2 = C_1 * CRateAng2(p2, q2, 1.0, t1)
                                  + C_z * CRateAng2(p2, q2, z2, t1)
                                  + C_zB * CRateAng2(p2, q2, zB2, t1);
        auto pPsi = q * Prefactor1 * PsiVals[ip][it]
                    - Prefactor2 * p * dExp * PsiVals[iq][it];

        if (std::abs(pPsi) < 1e-15) {
            continue;
        }

        Cint[it] += pPsi * Jac * Polynomial::pWeights[iq];
        #ifdef DEBUG

        #pragma omp critical

        if (!std::isfinite(Cint[it].real()) || !std::isfinite(Cint[it].imag())) {
            fmt::println(stderr, "Error: {}", Cint[it]);
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
        double p2 = p * p;

        for (int it = 0; it < Nt; it++) {
            double tm = Polynomial::time[it];
            std::complex<double> F = 0.0;
            double x = 0.0;

            for (int ix = 0; ix < Niter; ix++) {
                double t = 1.0 / (x + dx);

                if (t < tm) {
                    break;
                }

                double t2 = t * t;
                F = (F / dx + t2) /
                    (1.0 / dx + std::complex<double>(0.0, DeltaE(p2, t) * t2));


                x += dx;
            }

            PsiVals[ip][it] = -p2 * F;
            #ifdef DEBUG

            #pragma omp critical

            if (!std::isfinite(PsiVals[ip][it].real())
                    || !std::isfinite(PsiVals[ip][it].imag())) {
                fmt::println(stderr, "Error: {}", PsiVals[ip][it]);
                throw std::runtime_error("Error: Initial Psi is not finite");
            }

            #endif
        }
    }

    Integrate(0.0, CintVals, 0);
    PsiVals = CintVals;

}

void UpdatePsiVals(double dt, double Deltat, int itMin) {
    #ifdef EULER
    // Euler
    Integrate(Deltat, CintVals, itMin);

    for (int ip = 0; ip < Np; ip++) {
        for (int it = itMin; it < Nt; it++) {
            PsiVals[ip][it] -= dt * lambda * CintVals[ip][it];
        }
    }

    #endif

    // RK4
    #ifdef RK4
    CLEGENDRE k1, k2, k3, k4;
    double c2 = 0.5;
    double c3 = 0.5;
    double c4 = 1.0;
    double a21 = 0.5;
    double a32 = 0.5;
    double a43 = 1.0;
    double b1 = 1.0 / 6.0;
    double b2 = 1.0 / 3.0;
    double b3 = 1.0 / 3.0;
    double b4 = 1.0 / 6.0;

    CLEGENDRE PsiValsTemp = PsiVals;

    Integrate(t, k1);

    for (int i = 0; i < Np; i++) {
        PsiVals[i] = PsiValsTemp[i] - a21 * dt * lambda * k1[i];
    }

    Integrate(t + c2 * dt, k2);

    for (int i = 0; i < Np; i++) {
        PsiVals[i] = PsiValsTemp[i] - a32 * dt * lambda * k2[i];
    }

    Integrate(t + c3 * dt, k3);

    for (int i = 0; i < Np; i++) {
        PsiVals[i] = PsiValsTemp[i] - a43 * dt * lambda * k3[i];
    }

    Integrate(t + c4 * dt, k4);

    for (int i = 0; i < Np; i++) {
        PsiVals[i] = PsiValsTemp[i] - dt * lambda * (b1 * k1[i] + b2 * k2[i] + b3 *
                     k3[i] + b4 * k4[i]);
    }

    #endif


}

void Output(int i) {
    fmt::ostream file = fmt::output_file("Output/Expanding/bjorken-" +
                                         std::to_string(i) + ".dat");

    file.print("t p Re(Psi) Im(Psi) Re(Cint) Im(Cint) dE dEConst\n");

    for (int it = 0; it < Nt; it++) {
        double t = Polynomial::time[it];

        for (int ip = 0; ip < Np; ip++) {
            double p = Polynomial::pPoints[ip];
            double p2 = p * p;

            file.print("{} {} {} {} {} {} {} {} \n", Polynomial::time[it],
                       Polynomial::pPoints[ip], PsiVals[ip][it].real(), PsiVals[ip][it].imag(),
                       CintVals[ip][it].real(), CintVals[ip][it].imag(), DeltaE(p2, t), DeltaE(p2));
        }

        file.print("\n");
    }

    file.close();
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
            double p2 = p * p;
            double logExp = -DeltaE(p2, t, DeltatR);
            ra += p * Polynomial::pWeights[ip] * real(PsiVals[ip][it] *
                  std::complex<double>(std::cos(logExp), std::sin(logExp)));

            #ifdef DEBUG

            #pragma omp critical

            if (!std::isfinite(ra)) {
                fmt::println(stderr, "t: {} Deltat: {} dt: {} ra: {}", t, Deltat, dt, ra);
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
    std::string fname = "Output/Expanding/rate" + longName;
    fmt::ostream RateFile = fmt::output_file(fname + ".dat");
    // open(fname + ".dat", std::ios::out);
    RateFile.print("{}", Header);
    RateFile.print("t Rate \n0.0 0.0 \n");
    double Deltat = 0.0;
    int i = 0;


    int itDone = 0;

    while (Deltat <= tmax + dt) {

        #ifdef DEBUG

        if (Deltat >= Polynomial::time[itDone]) {
            fmt::println(stderr, "itDone: {} Deltat: {} time[itDone]: {}",
                         itDone, Deltat, Polynomial::time[itDone]);
            throw std::runtime_error("Error: Time index is not correct");
        }

        #endif

        UpdatePsiVals(dt, Deltat, itDone);
        ComputeRate(dt, Deltat, itDone);
        Deltat += dt;
        i++;

        for (int it = itDone; it < Nt; it++) {
            if (Deltat >= Polynomial::time[it]) {
                RateFile.print("{:.8} {:.8}\n", Polynomial::time[it], Rate[it]);
                RateFile.flush();
                itDone = it + 1;
                fmt::println(stderr, "itDone: {}", itDone);
            }
        }

    }

    RateFile.close();
}


void Setup() {
    Initialize();
    Output(0);
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
