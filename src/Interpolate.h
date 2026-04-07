#pragma once

#include <gsl/gsl_interp2d.h>
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
        xvec.resize(nx);
        yvec.resize(ny);
        zvec.resize(nx*ny);
        spline = gsl_interp2d_alloc(gsl_interp2d_bilinear, nx, ny);
    }
    ~Interpolate(){
        for(size_t i = 0; i < xacc.size(); i++){
            gsl_interp_accel_free(xacc[i]);
            gsl_interp_accel_free(yacc[i]);
        }
        gsl_interp2d_free(spline);
    }

    void setX(size_t i, double xi){
        xvec[i] = xi;
    }
    void setY(size_t j, double yj){
        yvec[j] = yj;
    }
    void setValues(size_t i, size_t j, double zij){
        gsl_interp2d_set(spline, zvec.data(), i, j, zij);
    }

    void init(){
        gsl_interp2d_init(spline, xvec.data(), yvec.data(), zvec.data(), xvec.size(), yvec.size());
    }

    double operator()(double x, double y){
        size_t thread_id = omp_get_thread_num();
        return gsl_interp2d_eval(spline, xvec.data(), yvec.data(), zvec.data(), x, y, xacc[thread_id], yacc[thread_id]);
    }


private:
    std::vector<gsl_interp_accel *> xacc, yacc;
    std::vector<double> xvec, yvec, zvec;
    gsl_interp2d *spline;
};
