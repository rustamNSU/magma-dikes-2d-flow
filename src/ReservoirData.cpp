#include "ReservoirData.hpp"
#include <cmath>
#include <highfive/H5Easy.hpp>

using Eigen::ArrayXXd;
using Eigen::ArrayXd;
using Eigen::ArrayXi;
using nlohmann::json;
using H5Easy::File;
using H5Easy::dump;


ReservoirData::ReservoirData(Mesh* mesh, json properties) :
    mesh(mesh),
    properties(properties)
{
    int nx = mesh->size();
    ny = properties["number_of_layers"];
    L = properties["reservoir_width"];
    C = properties["specific_heat_capacity"];
    k = properties["thermal_conductivity"];
    E = properties["youngs_modulus"];
    nu = properties["poisson_ratio"];
    g = properties["gravitational_acceleration"];
    KIc = properties["KIc"];
    if (properties["mesh_refinement_algorithm"] == "cosine") {
        generateCosineRefinementMesh();
    }

    density_model = properties["density_model"];
    temperature_model = properties["temperature_model"];
    calculateDensity();
    calculateReservoirTemperature();
    temperature = ArrayXXd::Zero(nx, ny);
    for (int iy = 0; iy < ny; iy++) {
        temperature.col(iy) = initial_temperature;
    }
}


const ArrayXd& ReservoirData::getInitialTemperature() const{
    return initial_temperature;
}


const ArrayXXd& ReservoirData::getTemperature() const{
    return temperature;
}


void ReservoirData::setTemperature(ArrayXXd&& T){
    temperature = std::move(T);
}


void ReservoirData::generateCosineRefinementMesh(){
    constexpr double pi = 3.14159265358979323846;
    ArrayXd Y = ArrayXd::LinSpaced(ny+1, 0.0, 1.0);
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
    if (density_model["name"] == "constant_density"){
        density = ArrayXd::Constant(nx, density_model["rho"]);
    }
    calculateLithostaticPressure();
}


void ReservoirData::calculateLithostaticPressure(){
    int nx = mesh->size();
    lithostatic_pressure = ArrayXd::Zero(nx);
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
    if (temperature_model["name"] == "constant_temperature_gradient"){
        double dT = temperature_model["temperature_gradient"];
        double Tmax = temperature_model["maximum_temperature"];
        double Tmin = temperature_model["minimum_temperature"];

        auto x = mesh->getx();
        Eigen::ArrayXd tmpT = -dT * x;
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


const Eigen::ArrayXd& ReservoirData::getLithostaticPressure() const{
    return lithostatic_pressure;
}


void ReservoirData::saveData(const std::string& savepath) const{
    File file(savepath, File::Overwrite);
    dump(file, "ny", ny);
    dump(file, "L", L);
    dump(file, "yc", yc);
    dump(file, "yb", yb);
    dump(file, "dy", dy);
    dump(file, "initial_temperature", initial_temperature);
    dump(file, "temperature", temperature);
    dump(file, "time", time);
}