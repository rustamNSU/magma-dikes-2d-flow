#include <iostream>
#include <cmath>
#include "Interp.hpp"


int main(int argc, char ** argv){
    std::array<double, 3> X = {-1, 0, 1};
    std::array<double, 3> Y = {0.5, 0, 0.5};
    std::cout << "dfdx(-2) = " << parab_derivative(-2, X, Y) << std::endl;
    std::cout << "dfdx(-1) = " << parab_derivative(-1, X, Y) << std::endl;
    std::cout << "dfdx(0) = " << parab_derivative(0, X, Y) << std::endl;
    std::cout << "dfdx(1) = " << parab_derivative(1, X, Y) << std::endl;
    std::cout << "dfdx(2) = " << parab_derivative(2, X, Y) << std::endl;
    return 0;
}