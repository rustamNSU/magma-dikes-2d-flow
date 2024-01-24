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


    /**
     * @brief Tridiagonal matrix algorithm for NxN system. See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
     * 
     * @param a N-vector, a[0] = 0
     * @param b N-vector
     * @param c N-vector, c[N-1] = 0
     * @param rhs N-vector
     * @return std::vector<double> 
     */
    std::vector<double> tridiagonal_solver(
        const std::vector<double> a,
        const std::vector<double> b,
        const std::vector<double> c,
        const std::vector<double> rhs
    );
}// Utils