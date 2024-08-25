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
    return scale * h2o * 0.01;
}


/**
 * @brief Return ideal gas density
 * @param p pressure in Pa
 * @param T temperature in C
 * @param R specific gas constant (for water vapor R = 461.52)
 * @return 
 */
inline double h2o_vapor_density(double p, double T, double R = 461.52){
    return p / (R * (T + 273.15));
}


/**
 * @brief return water saturated melt density
 * @param rhom0 melt phase density
 * @param rhow0 dissolved water density
 * @param gamma mass concentration of water into melt
 * @return 
 */
inline double h2o_sat_melt_density(double rhom0, double rhow0, double gamma){
    return rhom0 * rhow0 / (gamma * rhom0 + (1.0 - gamma) * rhow0);
}