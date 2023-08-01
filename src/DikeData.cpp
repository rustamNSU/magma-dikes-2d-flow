#include "DikeData.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;


DikeData::DikeData(Mesh* mesh) : mesh(mesh){
    int n = mesh->size();
    width = VectorXd::Zero(n);
    rho = VectorXd::Zero(n);
    pressure = VectorXd::Zero(n);
}