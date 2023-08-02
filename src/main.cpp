#include "Mesh.hpp"
#include "Elasticity.hpp"
#include "DikeData.hpp"
#include "InputData.hpp"
#include "MassBalance.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using Eigen::VectorXd;
using Eigen::MatrixXd;


int main(int argc, char ** argv){
    std::string input_path = argv[1];
	std::ifstream f(input_path);
    json input_json = json::parse(f);
	f.close();

    Mesh mesh(10, -1.0, 1.0);
    Elasticity elasticity(1.0, 0.2, &mesh);
    DikeData dike(&mesh);
    VectorXd vec = VectorXd::Ones(mesh.size());
    dike.setDensity(vec);
    dike.setWidth(vec);
    dike.setViscosity(vec);
    DikeData old_dike = dike;
    dike.setTime(10.0);
    InputData input(input_json);
    MassBalance mass_balance(
        &input,
        &elasticity,
        &mesh,
        nullptr
    );
    mass_balance.setNewTimestepData(&dike, &old_dike);
    mass_balance.solve();

    std::cout << mesh.getx() << std::endl;
    std::cout << mesh.getxl() << std::endl;
    std::cout << mesh.getxr() << std::endl;
    std::cout << elasticity.getMatrix() << std::endl;
    return 0;
}