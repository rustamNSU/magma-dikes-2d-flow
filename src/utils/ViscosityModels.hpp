#pragma once
#include <string>
#include <cmath>


/* Vogel–Fulcher–Tammann equation for temperature dependent viscosity with constant parameters */
inline double VFT_constant_viscosity(double T, double A, double B, double C){
    return std::exp(A + B / (T - C));
}