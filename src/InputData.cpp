#include "InputData.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using json = nlohmann::json;


InputData::InputData(const json& input, Mesh* mesh) : mesh(mesh){
    int n = mesh->size();
    E = input["reservoirProperties"]["E"];
    nu = input["reservoirProperties"]["nu"];
    g = input["reservoirProperties"]["g"];
    KIc = input["reservoirProperties"]["KIc"];
    density_model = input["reservoirProperties"]["densityModel"];
    if (density_model == "constant_density"){
        auto density_properties = input["reservoirProperties"]["constantDensity"];
        rhoR = VectorXd::Constant(n, density_properties["rho"]);
        calculateLithostaticPressure();
    }
}


const VectorXd& InputData::getPlith() const{
    return plith;
}


double InputData::getg() const{
    return g;
}


void InputData::calculateLithostaticPressure(){
    int n = mesh->size();
    plith = VectorXd::Zero(n);
    auto x = mesh->getx();
    auto xl = mesh->getxl();
    auto xr = mesh->getxr();
    double upper_weight = 0.0;
    for (int i = n-1; i >= 0; --i){
        double dx = xr[i] - xl[i];
        plith[i] = rhoR[i] * g * dx / 2.0 + upper_weight;
        upper_weight += rhoR[i] * g * dx;
    }
    return;
}