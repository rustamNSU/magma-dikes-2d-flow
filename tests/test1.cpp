#include <iostream>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>
#include <filesystem>
#include "ViscosityModels.hpp"
#include <cmath>


int main(int argc, char ** argv){
    std::array<double, 11> composition = {
        62.4, 0.55, 20.01, 0.03, 0.02, 3.22, 9.08, 3.52, 0.93, 0.12, 0.0
    };
    double wh2o = 0.02;
    double T = 1000.0;
    GiordanoViscosity model;
    model.setComposition(composition);
    double mu = model.calculateViscosity(wh2o, T);
    std::cout << "log(mu) = " << std::log10(mu) << std::endl;
    std::cout << "log(mu) = " << std::log10(model.calculateViscosity(0.0, T)) << std::endl;
    std::cout << "log(mu) = " << std::log10(model.calculateViscosity(0.0696, 900)) << std::endl;
    return 0;
}