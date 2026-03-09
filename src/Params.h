#include <cmath>
#include <cstring>
#include <fmt/base.h>
#include <iostream>
#include <omp.h>
#include <ostream>
const double CA = 3.0;
const double CF = 4.0 / 3.0;
#define THERMAL_MASS
#ifdef THERMAL_MASS
    const double mq = 1.0 / CA / 1.5;
    const double mg = 0.5;
#else
    const double mq = 0.0;
    const double mg = 0.0;
#endif
const double Nc = 3.0;
const double Nf = 3.0;

#define GREEN "[\033[1;32m"
#define RESET "\033[0m]"
#define BJORKEN
// #define DEBUG
extern double C_1, C_z, C_zB, m1, mz, mzB;
extern std::string process;

namespace WaveFct {
// Parameters //
extern double T0, P, z, zB, g, mDSqr, lambda, t0, t1, dt, tmax, alpha;
extern std::string Header, longName, shortName;

inline void Setup(int argc, char *argv[]) {
    std::string PStripedOfZeros, T0StripedOfZeros, zStripedOfZeros,
        t0StripedOfZeros, t1StripedOfZeros, alphaStripedOfZeros;
    // Parameters //
    T0 = 0.01;
    P = 4.0;
    z = 8.0 / 16.0;
    alpha = 1e-2;
    tmax = 20.0;
    double tmaxfm = 0.0;
    double InvGevTofm = 0.197326;
    double t0fm = 0.3;
    double t1fm = 0.0;
    std::string customName = "";
    std::string process = "G_GG";

    if (argc > 1) {
        for (int ia = 1; ia < argc; ia += 2) {
            if (std::strcmp(argv[ia], "-P") == 0) {
                P = std::stod(argv[ia + 1]);
            }

            else if (std::strcmp(argv[ia], "-z") == 0) {
                z = std::stod(argv[ia + 1]);
            }

            else if (std::strcmp(argv[ia], "-T") == 0) {
                T0 = std::stod(argv[ia + 1]);
            }

            else if (std::strcmp(argv[ia], "-a") == 0) {
                alpha = std::stod(argv[ia + 1]);
            }

            else if (std::strcmp(argv[ia], "-tmax") == 0) {
                tmax = std::stod(argv[ia + 1]);
            }

            else if (std::strcmp(argv[ia], "-tmaxfm") == 0) {
                tmaxfm = std::stod(argv[ia + 1]);
            }

            else if (std::strcmp(argv[ia], "-t0") == 0) {
                t0fm = std::stod(argv[ia + 1]);
            }

            else if (std::strcmp(argv[ia], "-t1") == 0) {
                t1fm = std::stod(argv[ia + 1]);
            }


            else if (std::strcmp(argv[ia], "-CN") == 0) {
                customName = std::string(argv[ia + 1]);
            }

            else if (std::strcmp(argv[ia], "-Process") == 0) {
                process = std::string(argv[ia + 1]);
            }

            else {
                fmt::println(stderr, "{} is not a valid option", argv[ia]);
                fmt::println(stderr, "Usage \n"
                                     "-Process: Process (G_GG, Q_GQ, G_QQ) \n"
                                     "-T: Temperature \n"
                                     "-P: Parent momentum \n"
                                     "-z: Momentum fraction \n"
                                     "-tmax: Maximum time \n"
                                     "-tmaxfm: Maximum time in fm \n"
                                     "-t0: Starting time of Bjorken expansion \n"
                                     "-a: Exponent \n"
                                     "-t1: The time the hard parton enters the medium \n");
                exit(1);
            }
        }

    }

    if (process == "G_GG") {
        C_1 = CA / 2.0;
        C_z = CA / 2.0;
        C_zB = CA / 2.0;
        m1 = mg;
        mz = mg;
        mzB = mg;
    }

    else if (process == "Q_GQ") {
        C_1 = CA / 2.0;
        C_z = (2.0 * CF - CA) / 2.0;
        C_zB = CA / 2.0;
        m1 = mq;
        mz = mg;
        mzB = mq;
    }

    else if (process == "G_QQ") {
        C_1 = (2.0 * CF - CA) / 2.0;
        C_z = CA / 2.0;
        C_zB = CA / 2.0;
        m1 = mg;
        mz = mq;
        mzB = mq;
    }

    else {
        fmt::print(stderr, "Process {} is not valid", process);
        exit(1);
    }

    zB = 1.0 - z;
    g = std::sqrt(0.3 * 4.0 * M_PI);
    mDSqr = g * g * T0 * T0 * (Nc / 3.0 + Nf / 6.0);
    lambda = (g * g * T0) * (2.0 * P * z * (1.0 - z)) / mDSqr;
    t0 = (t0fm / InvGevTofm) * mDSqr / (2.0 * P * z * (1.0 - z));
    t1 = (t1fm / InvGevTofm) * mDSqr / (2.0 * P * z * (1.0 - z));
    dt = 1e-3;

    if (tmaxfm != 0.0) {
        double tmaxtry = (tmaxfm / InvGevTofm) * mDSqr / (2.0 * P * z * (1.0 - z));
        tmax = tmaxtry;
        dt = tmax * 1e-4;

    }

    Header = "# T0 " + std::to_string(T0) +
             " P " + std::to_string(P) +
             " z " + std::to_string(z) +
             " g " + std::to_string(g) +
             " alpha " + std::to_string(alpha) +
             " t0 " + std::to_string(t0fm) +
             " t1 " + std::to_string(t1fm) +
             " Process " + process + "\n";

    fmt::println(stderr, "# T0 = {}{}{} P = {}{}{} "
                 "z = {}{}{} g = {}{}{} alpha = {}{}{} "
                 " t0 = {}{}{} t1 = {}{}{} Process = {}{}{} tmax = {}{}{} dt = {}{}{}",
                 GREEN,  T0,  RESET, GREEN,  P, RESET, GREEN,  z, RESET, GREEN,  g, RESET, GREEN,
                 alpha,
                 RESET, GREEN, t0fm, RESET, GREEN, t1fm, RESET, GREEN, process,
                 RESET, GREEN, tmax, RESET, GREEN, dt, RESET);

    fmt::println(stderr, "Number Of OMP Threads = {}{}{}", GREEN,
                 omp_get_max_threads(), RESET);

    PStripedOfZeros = std::to_string(P).substr(0,
                      std::to_string(P).find_last_not_of('0') + 1);
    T0StripedOfZeros = std::to_string(T0).substr(0,
                       std::to_string(T0).find_last_not_of('0') + 1);
    zStripedOfZeros = std::to_string(z).substr(0,
                      std::to_string(z).find_last_not_of('0') + 1);

    t0StripedOfZeros = std::to_string(t0fm).substr(0,
                       std::to_string(t0fm).find_last_not_of('0') + 1);

    t1StripedOfZeros = std::to_string(t1fm).substr(0,
                       std::to_string(t1fm).find_last_not_of('0') + 1);

    alphaStripedOfZeros = std::to_string(alpha).substr(0,
                          std::to_string(alpha).find_last_not_of('0') + 1);

    if (PStripedOfZeros[PStripedOfZeros.size() - 1] == '.') {
        PStripedOfZeros.pop_back();
    }

    if (T0StripedOfZeros[T0StripedOfZeros.size() - 1] == '.') {
        T0StripedOfZeros.pop_back();
    }

    if (zStripedOfZeros[zStripedOfZeros.size() - 1] == '.') {
        zStripedOfZeros.pop_back();
    }

    if (alphaStripedOfZeros[alphaStripedOfZeros.size() - 1] == '.') {
        alphaStripedOfZeros.pop_back();
    }

    if (t0StripedOfZeros[t0StripedOfZeros.size() - 1] == '.') {
        t0StripedOfZeros.pop_back();
    }

    if (t1StripedOfZeros[t1StripedOfZeros.size() - 1] == '.') {
        t1StripedOfZeros.pop_back();
    }


    if (customName != "") {
        shortName = customName;
        longName = customName;

    } else {
        shortName = "_P" + PStripedOfZeros + "_z" + zStripedOfZeros + "_T" +
                    T0StripedOfZeros;
        longName = "_P" + PStripedOfZeros
                   + "_z" + zStripedOfZeros + "_T" + T0StripedOfZeros
                   + "_alpha" + alphaStripedOfZeros
                   + "_t0" + t0StripedOfZeros
                   + "_t1" + t1StripedOfZeros;
    }


}

inline double Temp(double t) {

    #ifdef BJORKEN
    // BJORKEN EXPANSION
    return std::pow(t0 / (t + t1 + t0), alpha / 3.0);
    #else
    return 1.0;
    #endif
}

inline double CRate(double q2) {
    return 1.0 / (q2 * (q2 + 1));
}

inline double CRateAng(double p2, double q2, double z2) {
    if (p2 == q2) {
        return 0.0;
    }

    return
        1.0 / std::abs(p2 - q2)
        - 1.0 / std::sqrt((z2 + p2 + q2) * (z2 + p2 + q2) - 4.0 * p2 * q2);
}

inline double CRateAng2(double p2, double q2, double z2) {
    if (p2 == q2 || p2 == 0.0 || q2 == 0.0) {
        return 0.0;
    }

    return
        (0.5 / std::sqrt(p2 * q2)) *
        ((p2 + q2) / std::abs(p2 - q2) -
         (z2 + p2 + q2) / std::sqrt((z2 + p2 + q2) * (z2 + p2 + q2) - 4.0 * p2 * q2));
}

inline double CRate(double p2, double q2, double z2, double zB2) {
    return C_1 * CRateAng(p2, q2, 1.0)
           + C_z * CRateAng(p2, q2, z2)
           + C_zB * CRateAng(p2, q2, zB2);
}

inline double CRate2(double p2, double q2, double z2, double zB2) {
    return C_1 * CRateAng2(p2, q2, 1.0)
           + C_z * CRateAng2(p2, q2, z2)
           + C_zB * CRateAng2(p2, q2, zB2);
}

inline double CRateAng(double p2, double q2, double z2, double t) {
    if (p2 == q2) {
        return 0.0;
    }

    double TT = Temp(t);
    double T2 = TT * TT;
    double T3 = 1.0 / TT;
    p2 /= T2;
    q2 /= T2;

    return T3 *
           (1.0 / std::abs(p2 - q2)
            - 1.0 / std::sqrt((z2 + p2 + q2) * (z2 + p2 + q2) - 4.0 * p2 * q2));
}

inline double CRateAng2(double p2, double q2, double z2, double t) {
    if (p2 == q2 || p2 == 0.0 || q2 == 0.0) {
        return 0.0;
    }

    double TT = Temp(t);
    double T2 = TT * TT;
    double T3 = 1.0 / TT;
    p2 /= T2;
    q2 /= T2;
    return
        (T3 * 0.5 / std::sqrt(p2 * q2)) *
        ((p2 + q2) / std::abs(p2 - q2) -
         (z2 + p2 + q2) / std::sqrt((z2 + p2 + q2) * (z2 + p2 + q2) - 4.0 * p2 * q2));
}

// Collision kernel C(k) = 1/(k^2 * (k^2 + 1)) where k^2 = |p - q|^2
inline double CKernel(double k2, double mSqr, double t) {
    double TT = Temp(t);
    double T2 = TT * TT;
    double T3 = 1.0 / TT;
    k2 /= T2;
    return 1.0 / (k2 * (k2 + 1.0));
}

inline double Meff() {
    return zB * mz + z * mzB - z * zB * m1;
}

inline double Meff(double z, double zB) {
    return zB * mz + z * mzB - z * zB * m1;
}

inline double DeltaE(double p2) {
    return p2 + Meff();
}

inline double DeltaE(double p2, double t) {
    double Tt = Temp(t);
    return p2 + Meff() * Tt * Tt;
}

inline double DeltaE(double p2, double t, double Deltat) {

    #ifdef BJORKEN
    // BJORKEN EXPANSION
    // (a*t0**u*((t + t0)**u*(t + t0 - t1) - (t + t0)*(t + t0 - t1)**u)) / (((t + t0)*(t + t0 - t1))**u*(-1 + u))
    static const double u = alpha * 2.0 / 3.0;
    return p2 * Deltat
           + Meff() * std::pow(t0, u) *
           (std::pow(t + t0, u) * (t + t0 - Deltat)
            - (t + t0) * std::pow(t + t0 - Deltat, u)) /
           (std::pow((t + t0) * (t + t0 - Deltat), u) * (-1 + u));

    #else
    return (p2 + Meff()) * Deltat;
    #endif
}

}
