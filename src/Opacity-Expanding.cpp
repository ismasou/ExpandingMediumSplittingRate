#include "Polynomial.h"
#include <memory>
#define NQUADRATURE 701
#include "DoubleIntegral.h"
#include <omp.h>
#define BJORKEN
#include "Params.h"
#define EULER
double C_1, C_z, C_zB, m1, mz, mzB;
std::string process = "G_GG";

namespace WaveFct {
std::fstream RateFile;
double T0, P, z, zB, g, mDSqr, lambda, t0, t1, dt, tmax, alpha;
std::string Header, shortName, longName;
CLEGENDRE PsiVals;
double Rate = 0.0;
std::unique_ptr<DoubleGauusKronrod> doubleIntegrator;

std::complex<double> Integrand(const double q, const double p,
                               const int ip, const int iq, double t) {
    const double p2 = p * p;
    const double q2 = q * q;

    const double z2 = z * z;
    const double zB2 = zB * zB;

    const double Prefactor1 = C_1 * CRateAng(p2, q2, 1.0, t)
                              + C_z * CRateAng(p2, q2, z2, t)
                              + C_zB * CRateAng(p2, q2, zB2, t);
    const double Prefactor2 = C_1 * CRateAng2(p2, q2, 1.0, t)
                              + C_z * CRateAng2(p2, q2, z2, t)
                              + C_zB * CRateAng2(p2, q2, zB2, t);
    auto pPsi = q * Prefactor1 * PsiVals[ip]
                - Prefactor2 * p * PsiVals[iq];

    if (std::abs(pPsi) < 1e-15) {
        return 0.0;
    }

    auto res = pPsi;
    return res;
}

void Initialize() {

    const int Niter = 1e4;
    const double dx = 1e-3;

    for (int i = 0; i < Np; i++) {
        double p = Polynomial::pPoints[i];
        double p2 = p * p;
        std::complex<double> F = 0.0;
        std::complex<double> FOld = 0.0;
        double x = 0.0;
        double tol = 1e-15;

        for (int ix = 0; ix < Niter; ix++) {
            double t = 1.0 / (x + dx);
            double t2 = t * t;
            F = (F / dx + t2) /
                (1.0 / dx + std::complex<double>(0.0, DeltaE(p2, t) * t2));

            if (ix > 10 && std::abs((F - FOld) / F) < tol) {
                break;
            }

            FOld = F;
        }

        PsiVals[i] = F * p2;
    }
}

void Output(int i, double t) {
    std::fstream file;
    file.open("./Output/Evolution/output-" + std::to_string(i) + ".dat",
              std::ios::out);
    file << "p Re(Psi) Im(Psi) RateIntegrand dE dEConst" <<
         std::endl;

    for (int ip = 0; ip < Np; ip++) {
        double p = Polynomial::pPoints[ip];
        double p2 = p * p;
        double logExp = -DeltaE(p2) * t;
        auto RateIntegrand = p * real(PsiVals[ip] *
                                      std::complex<double>(std::cos(logExp), std::sin(logExp)));
        file << Polynomial::pPoints[ip] << " "
             << PsiVals[ip].real() << " " << PsiVals[ip].imag() << " "
             << RateIntegrand << " "
             << DeltaE(p2, t) << " "
             << DeltaE(p2)
             << std::endl;
    }

    file.close();
}

std::complex<double> IntegrateConverge(double p, double t, double dt) {

    double ta = t - dt;
    double p2 = p * p;
    double dE = DeltaE(p2);
    auto Fct = [&t, &dt, &p2](double x) {
        double t1 = t - dt + x * dt;
        double logExp = DeltaE(p2, t, t1);
        // double logExp = 0.0;
        // double sMin = t - t1;
        // double ds = t1;
        // for (int is = 0; is < GAUSSKONROD::N; is++) {
        //     double s = sMin + 0.5 * ds * (1.0 + GAUSSKONROD::x[is]);
        //     logExp += 0.5 * ds * GAUSSKONROD::w[is] * DeltaE(p2, s);
        // }

        return std::complex<double> {-dt *std::cos(-logExp), -dt *std::sin(-logExp)};
    };

    auto Shanks = [](double A0, double A1, double A2) {
        return A2 - (A2 - A1) * (A2 - A1) / (A2 - 2.0 * A1 + A0);
    };


    auto IntegrandOld = std::complex<double> {0.0, 0.0};
    auto TrueIntegrale = std::complex<double> {0.0, -(std::cos(dE * t) - std::cos(dE * ta)) / dE};
    // double TrueIntegrale = (1.0 - std::cos(dE * t)) / dE;
    // return TrueIntegrale;

    // Adaptive Trapzoidal Rule

    // We first calculate the integral with 1 point dx = 1/2
    // We do a change of variables from t1= [t-dt, t] =?
    // Int = 0.5 (F(a) + 2F(a+dx) + F(b)) * dx
    int N = 2;
    double dx = 1.0 / double(N);
    auto Integrale = (0.5 * Fct(0.0) + Fct(dx) + 0.5 * Fct(1.0)) * dx;

    double tol = 1e-7 * dE;

    static const int SnN = 3;
    std::array<std::complex<double>, SnN> An = {0.0};
    std::complex<double> Sn = 0.0;


    while (std::abs(Sn / IntegrandOld - 1.0) > tol || N < 1e2) {
        IntegrandOld = Sn;
        N *= 2;
        dx = 1.0 / double(N);
        Integrale *= 0.5;

        // USE TRAPEZOIDAL RULE with increasing number of points
        for (int ix = 1; ix < N; ix += 2) {

            double x = ix * dx;

            auto EXP = Fct(x);
            Integrale += EXP * dx;
        }


        Sn = Integrale;

        // An[0] = An[1];
        // An[1] = An[2];
        // An[2] = Integrale;
        // Sn = An[2] + (An[2] - An[1]) * (An[2] - An[1]) / (An[2] - 2.0 * An[1] + An[0]);
        // // if (!std::isfinite(Sn)) {
        // //     throw std::runtime_error("Error: Integral is not finite");
        // // }

        // double Sn0 = An[2] + (An[2] - An[1]) * (An[2] - An[1]) / (An[2] - 2.0 * An[1] + An[0]);
        // double Sn1 = An[3] + (An[3] - An[2]) * (An[3] - An[2]) / (An[3] - 2.0 * An[2] + An[1]);
        // double Sn2 = An[4] + (An[4] - An[3]) * (An[4] - An[3]) / (An[4] - 2.0 * An[3] + An[2]);
        // Sn = Sn2 + (Sn2 - Sn1) * (Sn2 - Sn1) / (Sn2 - 2.0 * Sn1 + Sn0);
        #pragma omp critical

        if (N > 1e8 || !std::isfinite(Sn.real()) || !std::isfinite(Sn.imag())){
            std::cout << "N: " << N << std::endl;
            std::cout << "p: " << p << " t: " << t << std::endl;
            std::cout << "dE: " << dE << std::endl;
            std::cout << "Integrand: " << Integrale << std::endl;
            std::cout << "Sn: " << Sn << std::endl;
            std::cout << "True Integrand: " << TrueIntegrale << std::endl;
            std::cout << "Tolerance: " << std::abs(Sn / IntegrandOld - 1.0) << std::endl;
            std::cout << "Error: " << std::abs(Sn / TrueIntegrale - 1.0) << std::endl;
            throw std::runtime_error("Error: Integral did not converge");
        }
    }

    Integrale = Sn;

    // std::cout << "N: " << N << std::endl;
    // std::cout << "Integrand: " << Integrand << std::endl;
    // std::cout << "Tolerance: " << std::abs(Integrand/IntegrandOld - 1.0) << std::endl;
    // double TrueIntegrand = -(1.0 - std::cos(dE * t)) / dE;
    // std::cout << "Error: " << std::abs(Integrand - TrueIntegrand) << std::endl;

    // N = 10;
    // dx = t / double(N - 1);
    // Integrand = (0.5 * 0.0 + 0.5 * t) * dx;
    //
    // for (int it1 = 1; it1 < N - 1; it1++) {
    //
    //     double t1 = it1 * dx;
    //
    //
    //     Integrand += t1 * dx;
    // }


    // #pragma omp critical
    // if (std::abs(Integrand / TrueIntegrale - 1.0) > 1e-3 * dE) {
    //     std::cout << "p: " << p << " t: " << t << std::endl;
    //     std::cout << "tIntegrand: " << Integrand << " TrueIntegrale: " << TrueIntegrale <<
    //               std::endl;
    //     std::cout << "Tol: " << std::abs(Integrand / IntegrandOld - 1.0) << std::endl;
    //     std::cout << "N: " << N << std::endl;
    //     throw std::runtime_error("Error: Integrand is not correct");
    // }

    return Integrale;
}

void ComputeRate(double t, double dt) {
    double ra = 0.0;


    // double OldRate = 0.0;
    int step = 1;

    #pragma omp parallel for reduction(+:ra)

    for (int ip = step; ip < Np; ip += step) {
        double p = Polynomial::pPoints[ip];
        auto tIntegrand = IntegrateConverge(p, t, dt);

        std::complex<double> qIntegrand = 0.0;
        // double p2 = p * p;
        // double logExp = DeltaE(p2) * t;
        // double tIntegrand = -(1.0 - std::cos(logExp)) / DeltaE(p2);

        for (int iq = 0; iq < Np; iq++) {
            double q = Polynomial::pPoints[iq];

            qIntegrand += Polynomial::pWeights[iq] * Integrand(q, p, ip, iq, t);

        }

        ra += p * Polynomial::pWeights[ip] * std::real(qIntegrand * tIntegrand);
        // double dp = p - Polynomial::pPoints[ip - step];
        // double NewRate = -p * qIntegrand;
        // ra += (OldRate + NewRate) / 2.0 * dp;
        // OldRate = NewRate;
    }

    Rate += ra / (4.0 * M_PI * M_PI);

}
void Evolve() {

    // Strip last part of string
    std::string fname = "Output/Opacity/opacity-bjorken" + longName;
    RateFile.open(fname + ".dat", std::ios::out);
    RateFile << "t Rate" << std::endl;
    Rate = 0.0;
    double t = 0.0;
    double dt = 1e-1;
    int i = 0;

    Rate = 0.0;
    RateFile << Header;
    RateFile << t << " " << Rate << std::endl;

    while (t <= tmax) {
        t += dt;
        i++;
        ComputeRate(t, dt);
        RateFile << t << " " << Rate << std::endl;
        // std::cout << "t: " << t << " Rate: " << Rate << std::endl;
        // dt *= 2.0;
    }

    RateFile.close();
}


void Setup() {
    doubleIntegrator = std::make_unique<DoubleGauusKronrod>();
    Initialize();
    Output(0, 0.0);
    Evolve();

    // std::cout << "Testing Integral " << IntegrateConverge(1.0, 1.0) << std::endl;
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
