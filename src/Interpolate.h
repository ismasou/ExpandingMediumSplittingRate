#pragma once

#include <gsl/gsl_spline2d.h>
#include <omp.h>
#include <vector>

class Interpolate{
public:

    Interpolate(size_t nx, size_t ny){
        xacc.resize(omp_get_max_threads());
        yacc.resize(omp_get_max_threads());
        for(size_t i = 0; i < xacc.size(); i++){
            xacc[i] = gsl_interp_accel_alloc();
            yacc[i] = gsl_interp_accel_alloc();
        }
        x.resize(nx);
        y.resize(ny);
        z.resize(nx*ny);
        spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
    }
    ~Interpolate(){
        for(size_t i = 0; i < xacc.size(); i++){
            gsl_interp_accel_free(xacc[i]);
        }
        gsl_spline2d_free(spline);
    }

    void setX(size_t i, double xi){
        x[i] = xi;
    }
    void setY(size_t j, double yj){
        y[j] = yj;
    }
    void setValues(size_t i, size_t j, double zij){
        gsl_spline2d_set(spline, z.data(), i, j, zij);
    }

    void init(){
        gsl_spline2d_init(spline, x.data(), y.data(), z.data(), x.size(), y.size());
    }

    double operator()(double x, double y){
        size_t thread_id = omp_get_thread_num();
        return gsl_spline2d_eval(spline, x, y, xacc[thread_id], yacc[thread_id]);
    }


private:
    std::vector<gsl_interp_accel *> xacc, yacc;
    std::vector<double> x, y, z;
    gsl_spline2d *spline;
};
