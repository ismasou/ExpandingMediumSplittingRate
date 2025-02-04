#ifndef DOUBLEINTEGRAL_H
#define DOUBLEINTEGRAL_H

#include <functional>
#include <vector>
#include <fstream>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#ifndef NQUADRATURE
    #define NQUADRATURE 51
#endif

class DoubleGauusKronrod {
public:

    static const int Nquadrature = NQUADRATURE; // !!! Only takes odd numbers !!! //
    static const int NLegendre = (Nquadrature - 1) / 2;

    DoubleGauusKronrod() {
        auto helper  = boost::math::quadrature::gauss<double, NLegendre>::abscissa();
        auto helperW = boost::math::quadrature::gauss<double, NLegendre>::weights();
        auto helper1  =
            boost::math::quadrature::gauss_kronrod<double, Nquadrature>::abscissa();
        auto helper1W =
            boost::math::quadrature::gauss_kronrod<double, Nquadrature>::weights();

        // Populate the positive values of the Gauss-quadrature
        for (int i = NLegendre / 2; i < NLegendre; i++)
        {
            GaussLegendrePoints[i] = helper[i - NLegendre / 2];
            GaussLegendreWeights[i] = helperW[i - NLegendre / 2];
        }

        // Populate the negative values of the Gauss-quadrature
        for (int i = 0; i < NLegendre / 2; i++)
        {
            GaussLegendrePoints [NLegendre / 2 - 1 - i] = -GaussLegendrePoints[NLegendre / 2
                + 1 + i];
            GaussLegendreWeights[NLegendre / 2 - 1 - i] = GaussLegendreWeights[NLegendre / 2
                + 1 + i];
        }

        // Populate the positive values of the Stieltjes-quadrature
        for (int i = Nquadrature / 2; i < Nquadrature; i++)
        {
            StieltjesPoints[i]  = helper1[i - Nquadrature / 2];
            StieltjesWeights[i] = helper1W[i - Nquadrature / 2];
        }

        // Populate the negative values of the Stieltjes-quadrature
        for (int i = 0; i < Nquadrature / 2; i++)
        {
            StieltjesPoints [Nquadrature / 2 - 1 - i]  = -StieltjesPoints [Nquadrature / 2 +
                1 + i];
            StieltjesWeights[Nquadrature / 2 - 1 - i]  =  StieltjesWeights[Nquadrature / 2 +
                1 + i];
        }

        // Populate the Double Weights for faster cubature
        for (size_t i = 0; i < Nquadrature; i++)
        {
            for (size_t j = 0; j < Nquadrature; j++)
            {
                StieltjesDoubleWeights[j + i * Nquadrature] = StieltjesWeights[i] *
                    StieltjesWeights[j];
            }
        }

        for (size_t i = 0; i < NLegendre; i++)
        {
            for (size_t j = 0; j < NLegendre; j++)
            {
                GaussLegendreDoubleWeights[j + i * NLegendre] = GaussLegendreWeights[i] *
                    GaussLegendreWeights[j];
            }
        }

        return;
    }

    double IntegralNdim4(double (*Integrand)(double *, size_t, void *),
                         std::array<double, 4>  const &a, std::array<double, 4>  const &b,
                         void *params) {

        for (size_t i = 0; i < a.size(); i++) {
            if (a[i] > b[i]) {
                std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
                throw std::runtime_error("Lower bound is larger than upper bound");
            }
        }

        double Sum = 0.0;
        size_t N = a.size();

        double x[4];

        for (size_t i0 = 0; i0 < Nquadrature; i0++) {
            x[0] = 0.5 * (b[0] + a[0]) + 0.5 * (b[0] - a[0]) * StieltjesPoints[i0];
            double w0 = StieltjesWeights[i0] * (b[0] - a[0]) * 0.5;

            for (size_t i1 = 0; i1 < Nquadrature; i1++) {
                x[1] = 0.5 * (b[1] + a[1]) + 0.5 * (b[1] - a[1]) * StieltjesPoints[i1];
                double w1 = StieltjesWeights[i1] * (b[1] - a[1]) * 0.5;

                for (size_t i2 = 0; i2 < Nquadrature; i2++) {
                    x[2] = 0.5 * (b[2] + a[2]) + 0.5 * (b[2] - a[2]) * StieltjesPoints[i2];
                    double w2 = StieltjesWeights[i2] * (b[2] - a[2]) * 0.5;

                    for (size_t i3 = 0; i3 < Nquadrature; i3++) {
                        x[3] = 0.5 * (b[3] + a[3]) + 0.5 * (b[3] - a[3]) * StieltjesPoints[i3];
                        double w3 = StieltjesWeights[i3] * (b[3] - a[3]) * 0.5;

                        Sum += w0 * w1 * w2 * w3 * Integrand(x, N, params);
                    }
                }
            }
        }

        return Sum;
    }

    double IntegralNdim(double (*Integrand)(double *, size_t, void *),
                        std::array<double, 3> const & a, std::array<double, 3> const & b,
                        void *params) {

        for (size_t i = 0; i < a.size(); i++) {
            if (a[i] > b[i]) {
                std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
                throw std::runtime_error("Lower bound is larger than upper bound");
            }
        }

        double Sum = 0.0;
        size_t N = a.size();

        double x[3];

        for (size_t i0 = 0; i0 < Nquadrature; i0++) {
            x[0] = 0.5 * (b[0] + a[0]) + 0.5 * (b[0] - a[0]) * StieltjesPoints[i0];
            double w0 = StieltjesWeights[i0] * (b[0] - a[0]) * 0.5;

            for (size_t i1 = 0; i1 < Nquadrature; i1++) {
                x[1] = 0.5 * (b[1] + a[1]) + 0.5 * (b[1] - a[1]) * StieltjesPoints[i1];
                double w1 = StieltjesWeights[i1] * (b[1] - a[1]) * 0.5;

                for (size_t i3 = 0; i3 < Nquadrature; i3++) {
                    x[2] = 0.5 * (b[2] + a[2]) + 0.5 * (b[2] - a[2]) * StieltjesPoints[i3];
                    double w2 = StieltjesWeights[i3] * (b[2] - a[2]) * 0.5;

                    Sum += w0 * w1 * w2 * Integrand(x, N, params);
                }
            }

        }

        return Sum;
    }

    double DoubleIntegral(std::function<double(double, double)> const &Integrand,
                          double a, double b, double a1, double b1, double &Error, double epsabs)
    {
        // Integrate a 2d std::function \int_{a}^{b} dx \int_{a1}^{b1}dy f(x,y)  //
        // The cubature is done using NLegendre points for Gauss-Legendre and (2*NLegendre + 1) points for Gauss-Stieltjes //
        // and the error is computed from the difference of the two //

        double kronrod_result = 0;
        double gauss_result = 0;
        double Jacobian = 0.25 * (b - a) * (b1 - a1);

        for (size_t i = 1; i < Nquadrature; i += 2)
        {
            double x = 0.5 * (b + a) + 0.5 * (b - a) * StieltjesPoints[i];

            for (size_t j = 1; j < Nquadrature; j += 2)
            {
                double y = 0.5 * (b1 + a1) + 0.5 * (b1 - a1) * StieltjesPoints[j];
                double w = StieltjesDoubleWeights[j + i * Nquadrature];
                double fct = Integrand(x, y);
                kronrod_result += w * fct;
                gauss_result += GaussLegendreDoubleWeights[(j / 2) + NLegendre * (i / 2)] * fct;
            }

            for (size_t j = 0; j < Nquadrature; j += 2)
            {
                double y = 0.5 * (b1 + a1) + 0.5 * (b1 - a1) * StieltjesPoints[j];
                double w = StieltjesDoubleWeights[j + i * Nquadrature];
                double fct = Integrand(x, y);
                kronrod_result += w * fct;
            }
        }

        for (size_t i = 0; i < Nquadrature; i += 2)
        {
            double x = 0.5 * (b + a) + 0.5 * (b - a) * StieltjesPoints[i];

            for (size_t j = 0; j < Nquadrature; j++)
            {
                double y = 0.5 * (b1 + a1) + 0.5 * (b1 - a1) * StieltjesPoints[j];
                double w = StieltjesDoubleWeights[j + i * Nquadrature];
                double fct = Integrand(x, y);
                kronrod_result += w * fct;

            }
        }

        kronrod_result *= Jacobian;
        gauss_result *= Jacobian;

        Error = std::abs(kronrod_result - gauss_result);

        if (std::abs(kronrod_result) > 1e-15) { Error /= std::abs(kronrod_result); }

        if (Error > epsabs) {
            std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
            std::cerr << "Error: " << Error << " epsabs: " << epsabs << "\n";
            throw std::runtime_error("Estimate of the cubature error " + std::to_string(
                                         Error) + " is larger than the target " + std::to_string(
                                         epsabs) + "\n Try larger target or using more Gauss points");
        }

        return kronrod_result;
    }

    double SingleIntegral(std::function<double(double, double)> const &Integrand,
                          double t, double a, double b, double &Error, double epsabs)
    {
        // Integrate a 1d with parameter t std::function \int_{a}^{b} dx f(x,y)  //
        // The cubature is done using NLegendre points for Gauss-Legendre and (2*NLegendre + 1) points for Gauss-Stieltjes //
        // and the error is computed from the difference of the two //

        double kronrod_result = 0;
        double gauss_result = 0;
        double Jacobian = 0.5 * (b - a);

        for (size_t i = 1; i < Nquadrature; i += 2)
        {
            double x = 0.5 * (b + a) + 0.5 * (b - a) * StieltjesPoints[i];
            double w = StieltjesWeights[i];
            double fct = Integrand(x, t);
            kronrod_result += w * fct;
            gauss_result += GaussLegendreWeights[(i / 2)] * fct;
        }

        for (size_t i = 0; i < Nquadrature; i += 2)
        {
            double x = 0.5 * (b + a) + 0.5 * (b - a) * StieltjesPoints[i];
            double w = StieltjesWeights[i];
            double fct = Integrand(x, t);
            kronrod_result += w * fct;
        }

        kronrod_result *= Jacobian;
        gauss_result *= Jacobian;

        Error = std::abs(kronrod_result - gauss_result);

        if (std::abs(kronrod_result) > 1e-15) { Error /= std::abs(kronrod_result); }

        if (Error > epsabs) {
            std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
            throw std::runtime_error("Estimate of the quadrature error " + std::to_string(
                                         Error) + " is larger than the target " + std::to_string(
                                         epsabs) + "\n Try larger target or using more Gauss points");
        }

        return kronrod_result;
    }

    ~DoubleGauusKronrod() {

    }

    std::array<double, NLegendre> GaussLegendrePoints;
    std::array<double, NLegendre> GaussLegendreWeights;
    std::array<double, NLegendre *NLegendre> GaussLegendreDoubleWeights;
    std::array<double, Nquadrature> StieltjesPoints;
    std::array<double, Nquadrature> StieltjesWeights;
    std::array<double, Nquadrature *Nquadrature> StieltjesDoubleWeights;
private:
};

#endif // !DEBUG
