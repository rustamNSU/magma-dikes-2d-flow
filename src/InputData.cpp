#include "InputData.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using json = nlohmann::json;


InputData::InputData(const json& input, Mesh* mesh) : mesh(mesh){
    this->input = input;
    int n = mesh->size();
    E = input["reservoirProperties"]["E"];
    nu = input["reservoirProperties"]["nu"];
    g = input["reservoirProperties"]["g"];
    KIc = input["reservoirProperties"]["KIc"];
    density_model = input["reservoirProperties"]["densityModel"];
    reservoir_temperature_model = input["reservoirProperties"]["temperatureModel"];
    if (density_model == "constant_density"){
        auto density_properties = input["reservoirProperties"]["constantDensity"];
        reservoir_density = VectorXd::Constant(n, density_properties["rho"]);
    }
    calculateLithostaticPressure();
    calculateReservoirTemperature();
}


const VectorXd& InputData::getLithostaticPressure() const{
    return lithostatic_pressure;
}


const VectorXd& InputData::getReservoirTemperature() const{
    return reservoir_temperature;
}


double InputData::getGravityAcceleration() const{
    return g;
}


/* return <E, nu, KIc> */
std::tuple<double, double, double> InputData::getElasticityParameters() const{
    return std::make_tuple(E, nu, KIc);
}


void InputData::calculateLithostaticPressure(){
    int n = mesh->size();
    lithostatic_pressure = VectorXd::Zero(n);
    auto x = mesh->getx();
    auto xl = mesh->getxl();
    auto xr = mesh->getxr();
    double upper_weight = 0.0;
    for (int i = n-1; i >= 0; --i){
        double dx = xr[i] - xl[i];
        lithostatic_pressure[i] = reservoir_density[i] * g * dx / 2.0 + upper_weight;
        upper_weight += reservoir_density[i] * g * dx;
    }
    return;
}


void InputData::calculateReservoirTemperature(){
    if (reservoir_temperature_model == "constant_gradient"){
        auto temperature_properties = input["reservoirProperties"]["constantTemperatureGradient"];
        double dT = temperature_properties["dT"];
        auto x = mesh->getx();
        reservoir_temperature = -dT * x;
    }
}