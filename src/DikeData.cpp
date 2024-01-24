#include "DikeData.hpp"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;


DikeData::DikeData(Mesh* mesh, int ny) : mesh(mesh), ny(ny){
    int n = mesh->size();
    width = VectorXd::Zero(n);
    density = VectorXd::Zero(n);
    pressure = VectorXd::Zero(n);
    overpressure = VectorXd::Zero(n);

    yc = MatrixXd::Zero(n, ny);
    yb = MatrixXd::Zero(n, ny+1);
    temperature = MatrixXd::Zero(n, ny);
    viscosity = MatrixXd::Zero(n, ny);
    mobility = MatrixXd::Zero(n, ny);
    total_mobility = VectorXd::Zero(n);
    total_face_mobility = VectorXd::Zero(n+1);
    face_flux = MatrixXd::Zero(n+1, ny);
    total_face_flux = VectorXd::Zero(n+1);
    time = 0.0;
}


double DikeData::getTime() const{
    return time;
}


int DikeData::getLayersNumber() const{
    return ny;
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


const MatrixXd& DikeData::getYc() const{
    return yc;
}


const MatrixXd& DikeData::getYb() const{
    return yb;
}


const MatrixXd& DikeData::getTemperature() const{
    return temperature;
}


const MatrixXd& DikeData::getMobility() const{
    return mobility;
}


const VectorXd& DikeData::getTotalMobility() const{
    return total_mobility;
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
    for (int ix = 0; ix < mesh->size(); ix++){
        double hw = width[ix] / 2.0;
        yb.row(ix) = VectorXd::LinSpaced(ny+1, 0, hw);
        yc.row(ix) = (yb(ix, Eigen::seq(1, ny)) + yb(ix, Eigen::seq(0, ny-1))) / 2.0;
    }
    return;
}


void DikeData::setPressure(const VectorXd& vec){
    pressure = vec;
    return;
}


void DikeData::setTemperature(const MatrixXd& mat){
    temperature = mat;
    return;
}


void DikeData::setViscosity(const MatrixXd& mat){
    viscosity = mat;
    return;
}