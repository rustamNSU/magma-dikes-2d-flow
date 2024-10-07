#pragma once


#include <cmath>
#include <array>


inline double parab_derivative(double x, const std::array<double, 3>& X, const std::array<double, 3>& Y){
    return Y[0] * (2*x - X[1] - X[2]) / ((X[0]-X[1]) * (X[0]-X[2])) +
           Y[1] * (2*x - X[2] - X[0]) / ((X[1]-X[2]) * (X[1]-X[0])) +
           Y[2] * (2*x - X[0] - X[1]) / ((X[2]-X[0]) * (X[2]-X[1]));
}