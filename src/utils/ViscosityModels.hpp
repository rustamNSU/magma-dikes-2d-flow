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
    return std::exp(A + B / (T - C));
}