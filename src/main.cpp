#include <iostream>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>
#include <filesystem>
#ifdef USE_OMP
#include <Eigen/Core>
#endif

#include "InputData.hpp"
#include "Mesh.hpp"
#include "ReservoirData.hpp"
#include "Elasticity.hpp"
#include "DikeData.hpp"
#include "DikeModel2d.hpp"

using json = nlohmann::json;
using Eigen::VectorXd;
using Eigen::MatrixXd;
namespace fs = std::filesystem;



int main(int argc, char ** argv){
    #ifdef USE_OMP
    Eigen::initParallel();
    #endif
    const std::filesystem::path workdir = std::filesystem::current_path();
    std::string runner_path = argv[1];
	std::ifstream f(runner_path);
    json runner = json::parse(f);
	f.close();
    auto input_paths = runner["input_paths"].get<std::vector<std::string>>();
    for (const auto& input_path : input_paths){
        DikeModel2d model(input_path);
        model.run();
    }
    return 0;
}