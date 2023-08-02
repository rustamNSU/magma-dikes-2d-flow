#pragma once
#include <vector>

namespace Utils{
    /**
     * @brief Return mean from a to b of piecewise constant function.
     * If b > X(end), return Y(end)
     * 
     * @param X
     * @param Y Y[i] from X[i] to X[i+1]
     * @param a 
     * @param b
     */
    double mean_piecewise(
        const std::vector<double> &X, 
        const std::vector<double> &Y,
        double a, double b
    );
}// Utils