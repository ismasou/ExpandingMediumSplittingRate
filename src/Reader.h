#pragma once

#include <vector>
#include <gsl/gsl_interp2d.h>
#include <omp.h>
class Reader {
private:
    // Scaling functions
    std::vector<double> qVals, thetaVals, omegaVals, muVals,
        TStartVals, EpsVals;
    std::vector<std::vector<double>> CVals;
    size_t Nq, Ntheta, Nomega;
    std::vector<gsl_interp_accel *> qacc, thetacc, omegaacc;
    std::vector<gsl_interp2d *> Cinterp;

public:

    Reader() {
        //TODO: Setup code to read files and initialize the grids
        // For now we use dummy values
        Nq = 16;
        Ntheta = 16;
        Nomega = 16;
        qVals.resize(Nq);
        thetaVals.resize(Ntheta);
        omegaVals.resize(Nomega);
        CVals.resize(Nomega);
        for (size_t i = 0; i < Nomega; ++i) {
            CVals[i].resize(Nq * Ntheta);
        }
        muVals.resize(Nomega);
        TStartVals.resize(Nomega);
        EpsVals.resize(Nomega);
        size_t numThreads = omp_get_max_threads();

        for (size_t i = 0; i < Nq; ++i) {
            qVals[i] = 0.1 * (i + 1);
        }

        for (size_t i = 0; i < Ntheta; ++i) {
            thetaVals[i] = 0.1 * (i + 1);
        }

        for (size_t i = 0; i < Nomega; ++i) {
            omegaVals[i] = 0.1 * (i + 1);
            muVals[i] = 0.01 * (i + 1);
            TStartVals[i] = 0.1 * (i + 1);
            EpsVals[i] = 0.001 * (i + 1);
        }

        Cinterp.resize(numThreads);
        qacc.resize(numThreads);
        thetacc.resize(numThreads);
        omegaacc.resize(numThreads);

        for (size_t tID = 0; tID < numThreads; ++tID) {
            qacc[tID] = gsl_interp_accel_alloc();
            thetacc[tID] = gsl_interp_accel_alloc();
            omegaacc[tID] = gsl_interp_accel_alloc();
            Cinterp[tID] = gsl_interp2d_alloc(gsl_interp2d_bilinear, Nq, Ntheta);
        }
    }

    ~Reader() {
        for (size_t tID = 0; tID < omegaacc.size(); ++tID) {
            gsl_interp_accel_free(qacc[tID]);
            gsl_interp_accel_free(thetacc[tID]);
            gsl_interp_accel_free(omegaacc[tID]);
            gsl_interp2d_free(Cinterp[tID]);
        }
    }

    inline double Cqperp(double qperp, double theta, double omega) {
        size_t tID = omp_get_thread_num();

        if (omega <= omegaVals.front()) {
            return gsl_interp2d_eval(Cinterp[tID], qVals.data(), thetaVals.data(),
                                     CVals[0].data(), qperp, theta, qacc[tID], thetacc[tID]);
        }
        if (omega >= omegaVals.back()) {
            // Before the initial step
            return 0.0;
        }

        size_t i0 = gsl_interp_accel_find(omegaacc[tID], omegaVals.data(), Nomega, omega);
        size_t i1 = i0 + 1;

        double C0 = gsl_interp2d_eval(Cinterp[tID], qVals.data(), thetaVals.data(),
                                      CVals[i0].data(), qperp, theta, qacc[tID], thetacc[tID]);
        double C1 = gsl_interp2d_eval(Cinterp[tID], qVals.data(), thetaVals.data(),
                                      CVals[i1].data(), qperp, theta, qacc[tID], thetacc[tID]);

        double w = (omega - omegaVals[i0]) / (omegaVals[i1] - omegaVals[i0]);
        return (1.0 - w) * C0 + w * C1;
    }

};
