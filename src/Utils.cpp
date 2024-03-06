#include "Utils.hpp"
#include <cmath>
#include <stdexcept>


/**
 * @brief Return mean from a to b of piecewise constant function.
 * If b > X(end), return Y(end)
 * 
 * @param X
 * @param Y Y[i] from X[i] to X[i+1]
 * @param a 
 * @param b
 */
double Utils::mean_piecewise(
    const std::vector<double> &X, 
    const std::vector<double> &Y,
    double a, double b
){
    if (a > X.back()){
        return Y.back();
    }
    if (a < X.front()){
        throw std::invalid_argument("received wrong arguments, a < X.front()");
    }
    std::vector<int> inner_ind;
    std::vector<double> inner_x;

    /* Check if X[i] in [a, b) */
    for (int i = 0; i < X.size(); ++i){
        if (X[i] >= a && X[i] < b){
            inner_ind.push_back(i);
            inner_x.push_back(X[i]);
        }
    }

    /* If inner_ind is empty than X[j] < a < b < X[j+1] */
    if (inner_ind.size() == 0){
        int j = -1;
        for (int i = 0; i < X.size(); ++i){
            if (a > X[i]){
                ++j;
            } else{
                break;
            }
        }
        if (j == -1){
            return Y.front();
        }
        return Y[j];
    }

    double s = Y[inner_ind.back()] * (b - inner_x.back());
    double xLeft = a;
    for (int i = 0; i < inner_ind.size(); ++i){
        double xRight = inner_x[i];
        int element = inner_ind[i] - 1;
        s += (xRight - xLeft) * Y[element];
        xLeft = xRight;
    }
    double result = s / (b - a);
    return result;
}


/**
 * @brief Tridiagonal matrix algorithm for NxN system. See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 * 
 * @param a N-vector, a[0] = 0
 * @param b N-vector
 * @param c N-vector, c[N-1] = 0
 * @param rhs N-vector
 * @return std::vector<double> 
 */
std::vector<double> Utils::tridiagonal_solver(
    const std::vector<double> a,
    const std::vector<double> b,
    const std::vector<double> c,
    const std::vector<double> rhs
){
    int n = a.size();
    if (n != b.size() ||
        n != c.size() ||
        n != rhs.size() ||
        n < 2
    ){
        throw std::invalid_argument("received wrong arguments\n");
    }
    auto cp = a;
    auto dp = a;
    auto x = a;
    cp[0] = c[0] / b[0];
    for (int i = 1; i < n - 1; i++){
        cp[i] = c[i] / (b[i] - a[i] * cp[i-1]);
    }
    dp[0] = rhs[0] / b[0];
    for (int i = 1; i < n; i++){
        dp[i] = (rhs[i] - a[i] * dp[i-1]) / (b[i] - a[i] * cp[i-1]);
    }
    x[n-1] = dp[n-1];
    for (int i = n-2; i >= 0; i--){
        x[i] = dp[i] - cp[i] * x[i+1];
    }
    return x;
}