#include "ReservoirData.hpp"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using nlohmann::json;


ReservoirData::ReservoirData(Mesh* mesh, json properties) :
    mesh(mesh),
    properties(properties)
{
    int nx = mesh->size();
    ny = properties["numberOfLayers"];
    L = properties["reservoirWidth"];
    C = properties["specificHeatCapacity"];
    k = properties["thermalConductivity"];
    E = properties["E"];
    nu = properties["nu"];
    g = properties["g"];
    KIc = properties["KIc"];
    if (properties["meshRefinementAlgorithm"] == "cosine"){
        generateCosineRefinementMesh();
    }

    density_model = properties["densityModel"];
    temperature_model = properties["temperatureModel"];
    calculateDensity();
    calculateReservoirTemperature();
    temperature = MatrixXd::Zero(nx, ny);
    for (int iy = 0; iy < ny; iy++){
        temperature.col(iy) = initial_temperature;
    }
}


const VectorXd& ReservoirData::getInitialTemperature() const{
    return initial_temperature;
}


const MatrixXd& ReservoirData::getTemperature() const{
    return temperature;
}


void ReservoirData::setTemperature(MatrixXd&& T){
    temperature = std::move(T);
    return;
}


void ReservoirData::generateCosineRefinementMesh(){
    constexpr double pi = 3.14159265358979323846;
    VectorXd Y = VectorXd::LinSpaced(ny+1, 0.0, 1.0);
    yb = Y.unaryExpr([&](double x){
        return L * (1.0 - std::cos(0.5 * pi * x));
    });
    dy = yb(Eigen::seq(1, ny)) - yb(Eigen::seq(0, ny-1));
    yc = 0.5 * (yb(Eigen::seq(1, ny)) + yb(Eigen::seq(0, ny-1)));
}



/* return <E, nu, KIc> */
std::tuple<double, double, double> ReservoirData::getElasticityParameters() const{
    return std::make_tuple(E, nu, KIc);
}


void ReservoirData::calculateDensity(){
    int nx = mesh->size();
    if (density_model == "constant_density"){
        auto density_properties = properties["constantDensity"];
        density = VectorXd::Constant(nx, density_properties["rho"]);
    }
    calculateLithostaticPressure();
}


void ReservoirData::calculateLithostaticPressure(){
    int nx = mesh->size();
    lithostatic_pressure = VectorXd::Zero(nx);
    auto x = mesh->getx();
    auto xl = mesh->getxl();
    auto xr = mesh->getxr();
    double upper_weight = -xr(Eigen::last) * g * density(Eigen::last);
    for (int i = nx-1; i >= 0; --i){
        double dx = xr[i] - xl[i];
        lithostatic_pressure[i] = density[i] * g * dx / 2.0 + upper_weight;
        upper_weight += density[i] * g * dx;
    }
    return;
}


void ReservoirData::calculateReservoirTemperature(){
    if (temperature_model == "constant_gradient"){
        auto temperature_properties = properties["constantTemperatureGradient"];
        double dT = temperature_properties["dT"];
        double Tmax = temperature_properties["maximum_temperature"];
        double Tmin = temperature_properties["minimum_temperature"];

        auto x = mesh->getx();
        Eigen::VectorXd tmpT = -dT * x;
        initial_temperature = tmpT.unaryExpr([&](double x){
            if (x > Tmax){
                return Tmax;
            }
            else if (x < Tmin){
                return Tmin;
            }
            else{
                return x;
            }
        });
    }
}


double ReservoirData::getGravityAcceleration() const{
    return g;
}


const Eigen::VectorXd& ReservoirData::getLithostaticPressure() const{
    return lithostatic_pressure;
}