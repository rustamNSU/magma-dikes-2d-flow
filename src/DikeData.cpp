#include "DikeData.hpp"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;


DikeData::DikeData(Mesh* mesh) : mesh(mesh){
    int n = mesh->size();
    width = VectorXd::Zero(n);
    density = VectorXd::Zero(n);
    pressure = VectorXd::Zero(n);
    overpressure = VectorXd::Zero(n);
    viscosity = VectorXd::Zero(n);
    mobility = VectorXd::Zero(n+1);
    time = 0.0;
}


double DikeData::getTime() const{
    return time;
}


const VectorXd& DikeData::getMobility() const{
    return mobility;
}


const VectorXd& DikeData::getDensity() const{
    return density;
}


const VectorXd& DikeData::getWidth() const{
    return width;
}


const VectorXd& DikeData::getPressure() const{
    return pressure;
}


void DikeData::setTime(double time){
    this->time = time;
    return;
}


void DikeData::setDensity(const VectorXd& vec){
    density = vec;
    return;
}


void DikeData::setWidth(const VectorXd& vec){
    width = vec;
    return;
}


void DikeData::setViscosity(const VectorXd& vec){
    viscosity = vec;
    return;
}


void DikeData::setPressure(const VectorXd& vec){
    pressure = vec;
    return;
}


void DikeData::calculateMobility(){
    int n = mesh->size();
    mobility.fill(0.0);
    for (int i = 1; i < n; ++i){
        double ml = width[i-1] <= MIN_MOBILITY_WIDTH ? 0.0 : std::pow(width[i-1], 3) / 12.0 / viscosity[i-1];
        double mr = width[i] <= MIN_MOBILITY_WIDTH ? 0.0 : std::pow(width[i], 3) / 12.0 / viscosity[i];
        mobility[i] = (ml + mr) / 2.0;
    }
    return;
}