#pragma once
#include <cmath>


/**
 * @brief H2O in wt%. See eq. (1) from Yan Lavallée, Thermal vesiculation during volcanic eruptions, 10.1038/nature16153
 * @param p pressure in Pa
 * @param T temperature in C
 * @param scale scale coefficient
 * @return 
 */
inline double h2o_wt_lavallee2015(double p, double T, double scale = 1.0){
    double pmpa = p * 1e-6;
    double psqrt = std::sqrt(pmpa);
    double ppsqrt = pmpa * psqrt;
    double h2o = (354.94*psqrt + 9.623*pmpa - 1.5233*ppsqrt) / (T + 273.15) + 0.0012439*ppsqrt;
    return scale * h2o;
}