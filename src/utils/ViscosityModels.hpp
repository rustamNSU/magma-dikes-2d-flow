#pragma once
#include <string>
#include <cmath>


/**
 * @brief Vogel–Fulcher–Tammann equation for temperature dependent viscosity with constant parameters
 * 
 * @param T [K]
 * @param A 
 * @param B 
 * @param C [K]
 * @return double 
 */
inline double VFT_constant_viscosity(double T, double A, double B, double C){
    constexpr double ln10 = 2.302585092994045684;
    return std::exp(ln10 * (A + B / (T - C)));
}