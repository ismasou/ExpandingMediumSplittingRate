#include "Polynomial.h"
#include "ChebyshevTheta.h"
#include <memory>
#define NQUADRATURE 701
#include <omp.h>
#include "Params.h"
#include <fmt/os.h>
#include "CRateTable.h"
#define EULER
// #define RK4
#define DEBUG

double C_1, C_z, C_zB, m1, mz, mzB;
std::string process = "G_GG";
namespace WaveFct {
double T0, P, z, zB, g, mDSqr, lambda, t0, t1, dt, tmax, alpha;
std::string Header, longName, shortName;
std::unique_ptr<fmt::ostream> RateFile;
CLEGENDRE PsiVals, PsiValsOld, CintVals;
double Rate, RateAverage, RateOld, RateL = 0.0, RateR;
int NAverage = 0;
const int NAv = 20;
bool RateNotConverged = true;


double CRate(double q2, double theta){
    if (q2 == 0.0) {
	return 0.0;
    }
    // return 1.0 / (q2 * (q2 + 1.0));
    return CRateInterp(q2, theta);
}

inline std::complex<double> Integrand(const double q, const double p,
                                      const double Deltat, const int ip, const int iq) {


    const double Jac = 1.0 / (2.0 * M_PI);
    const double p2 = p * p;
    const double q2 = q * q;

    const double z2 = z * z;
    const double zB2 = zB * zB;

    const double logExp = (p2 - q2) * Deltat;
    auto dExp = std::complex<double>(std::cos(logExp), std::sin(logExp));


    auto integrand = [&](double theta, double cosTheta, double sinTheta) {

        double k2 = p2 - 2.0 * p *q *cosTheta + q2;

        double C = C_1  * CRate(k2, theta)
                 + C_z  * CRate(k2 / z2,  theta) / z2
                 + C_zB * CRate(k2 / zB2, theta) / zB2;

        return std::array<double, 2> {C,
                                      C * cosTheta
                                     };
    };

    auto [thetaInt1, thetaInt2] = ThetaQuad::Integrate<2>(integrand);
    thetaInt1 /= (2.0 * M_PI);
    thetaInt2 /= (2.0 * M_PI);

    auto pPsi = q * thetaInt1 * PsiVals[ip]
                - thetaInt2 * p * dExp * PsiVals[iq];

    if (std::abs(pPsi) < 1e-15) {
        return 0.0;
    }

    auto res = Jac * pPsi * Polynomial::pWeights[iq];

    if (!std::isfinite(res.real()) || !std::isfinite(res.imag())) {
        fmt::println(stderr, "p: {} q: {} ip: {} iq: {} Deltat: {} res: {} + i {}", p,
                     q, ip, iq, Deltat, res.real(), res.imag());
        throw std::runtime_error("Error: Integrand is not finite");
    }

    return res;
}

void Integrate(double Deltat, CLEGENDRE &CintVals) {

    #pragma omp parallel for

    for (int ip = 0; ip < Np; ip++) {
        double p = Polynomial::pPoints[ip];
        CintVals[ip] = 0.0;

        for (int iq = 0; iq < Np; iq++) {
            double q = Polynomial::pPoints[iq];

            CintVals[ip] += Integrand(q, p, Deltat, ip, iq);
        }
    }
}

void Initialize() {
    for (int i = 0; i < Np; i++) {
        double p = Polynomial::pPoints[i];
        double p2 = p * p;
        PsiVals[i] = {0.0, p2 / DeltaE(p2)};
    }

    Integrate(0.0, CintVals);

    for (int i = 0; i < Np; i++) {
        PsiVals[i] = CintVals[i];
    }
}

void UpdatePsiVals(double dt, double Deltat) {
    #ifdef EULER
    // Euler
    Integrate(Deltat, CintVals);

    std::copy(PsiVals.begin(), PsiVals.end(), PsiValsOld.begin());
    #pragma omp parallel for

    for (int i = 0; i < Np; i++) {
        PsiVals[i] -= dt * lambda * CintVals[i];
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

    Integrate(Deltat, k1);

    for (int i = 0; i < Np; i++) {
        PsiVals[i] = PsiValsTemp[i] - a21 * dt * lambda * k1[i];
    }

    Integrate(Deltat + c2 * dt, k2);

    for (int i = 0; i < Np; i++) {
        PsiVals[i] = PsiValsTemp[i] - a32 * dt * lambda * k2[i];
    }

    Integrate(Deltat + c3 * dt, k3);

    for (int i = 0; i < Np; i++) {
        PsiVals[i] = PsiValsTemp[i] - a43 * dt * lambda * k3[i];
    }

    Integrate(Deltat + c4 * dt, k4);

    for (int i = 0; i < Np; i++) {
        PsiVals[i] = PsiValsTemp[i] - dt * lambda * (b1 * k1[i] + b2 * k2[i] + b3 *
                     k3[i] + b4 * k4[i]);
    }

    #endif


}

void Output(int i, double t) {
    fmt::ostream file = fmt::output_file("./Output/Evolution/output-" +
                                         std::to_string(i) + ".dat",
                                         fmt::file::WRONLY | fmt::file::CREATE);
    file.print("p Re(Psi) Im(Psi) Re(Cint) Im(Cint) RateIntegrand \n");

    for (int ip = 0; ip < Np; ip++) {
        double p = Polynomial::pPoints[ip];
        double p2 = p * p;
        double logExp = -DeltaE(p2) * t;
        auto RateIntegrand = p * real(PsiVals[ip] *
                                      std::complex<double>(std::cos(logExp), std::sin(logExp)));
        file.print("{} {} {} {} {} {}\n",
                   Polynomial::pPoints[ip], PsiVals[ip].real(), PsiVals[ip].imag(),
                   CintVals[ip].real(), CintVals[ip].imag(), RateIntegrand);
    }

    file.close();
}

void ComputeRate(double dt, double Deltat) {
    double ra = 0.0;

    #pragma omp parallel for reduction(+:ra)

    for (int ip = 0; ip < Np; ip++) {
        double p = Polynomial::pPoints[ip];
        double p2 = p * p;
        double dE = DeltaE(p2);

        double logExp = -dE * Deltat;
        ra += p * Polynomial::pWeights[ip] * real(PsiVals[ip] *
              std::complex<double>(std::cos(logExp), std::sin(logExp)));
    }

    RateR = ra / (2.0 * M_PI);

    Rate += 0.5 * dt * (RateL + RateR);

    RateL = RateR;

}

void Evolve() {

    std::string fname = "Output/Static2D/rate" + shortName;

    RateFile = std::make_unique<fmt::ostream>(fmt::output_file(fname + ".dat",
               fmt::file::WRONLY | fmt::file::CREATE));
    RateFile->print("{}", Header);
    RateFile->print("t Rate \n");
    RateFile->flush();
    Rate = 0.0;
    double Deltat = 0.0;
    int i = 0;


    while (Deltat <= tmax) {
        #ifdef DEBUG

        if (i % 100 == 0) {
            Integrate(Deltat, CintVals);
            Output(i / 100, Deltat);
        }

        #endif

        if (RateNotConverged) {
            UpdatePsiVals(dt, Deltat);
            ComputeRate(dt, Deltat);
        }

        if (i % 10 == 0) {
            RateFile->print("{} {} {}\n", Deltat, Rate, RateOld - Rate);
            RateFile->flush();
            // Lets check how big is the wave function
            double norm = 0.0;

            for (int ip = 0; ip < Np; ip++) {
                norm += std::norm(PsiVals[ip]) * Polynomial::pWeights[ip];
            }

            if (norm < 1e-9) {
                fmt::println(stderr, "Wave function converged to zero");
                break;
            }
        }

        RateOld = Rate;

        Deltat += dt;
        i++;
    }

    RateFile->close();
}


void Setup() {
    InitCRateTable("/pc2/users/h/hion0024/QCDKinetic/Data/Cqperp/OUTPUT10/CqperpXgg_gg_0.txt", std::sqrt(mDSqr), 3, g, T0);
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
