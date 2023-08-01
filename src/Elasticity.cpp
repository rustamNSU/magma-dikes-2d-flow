#include "Elasticity.hpp"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

Elasticity::Elasticity(
    double E,
    double nu,
    Mesh* mesh
) : E(E), nu(nu), mesh(mesh)
{
    Ep = E / (1.0 - nu*nu);
    generateMatrix();
}


void Elasticity::generateMatrix(){
    int n = mesh->size();
    auto x = mesh->getx();
    auto xl = mesh->getxl();
    auto xr = mesh->getxr();
    matrix.resize(n, n);
    for (int r = 0; r < n; ++r){
        for (int c = 0; c < n; ++c){
            matrix(r,c) = -Ep / (4.0 * M_PI) *
                (1.0 / (x(r) - xr(c)) -
                 1.0 / (x(r) - xl(c)));
        }
    }
    return;
}


const MatrixXd& Elasticity::get_matrix() const{
    return matrix;
}