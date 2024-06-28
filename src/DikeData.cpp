#include "DikeData.hpp"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using nlohmann::json;


DikeData::DikeData(Mesh* mesh, const json& alg_properties) :
    meshX(mesh),
    algorithm_properties(alg_properties)
{
    ny = algorithm_properties["numberOfLayers"];
    int nx = meshX->size();
    width = VectorXd::Zero(nx);
    hw = VectorXd::Zero(nx);
    pressure = VectorXd::Zero(nx);
    overpressure = VectorXd::Zero(nx);

    yb = VectorXd::LinSpaced(ny+1, 0.0, 1.0);;
    yc = 0.5 * (
        yb(Eigen::seq(0, Eigen::last - 1)) + 
        yb(Eigen::seq(1, Eigen::last))
    );
    
    density = MatrixXd::Zero(nx, ny);
    temperature = MatrixXd::Zero(nx, ny);
    viscosity = MatrixXd::Zero(nx, ny);
    qx = MatrixXd::Zero(nx+1, ny);
    A = MatrixXd::Zero(nx+1, ny);
    C = MatrixXd::Zero(nx+1, ny);
    qy = MatrixXd::Zero(nx, ny+1);
    mobility = MatrixXd::Zero(nx, ny);
    Qx = VectorXd::Zero(nx+1);
}