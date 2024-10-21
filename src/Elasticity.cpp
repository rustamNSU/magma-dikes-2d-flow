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
    matrix = MatrixXd::Zero(n, n);
    A = MatrixXd::Zero(n, n);
    B = MatrixXd::Zero(n, n);
    for (int r = 0; r < n; ++r){
        for (int c = 0; c < n; ++c){
            matrix(r,c) = -Ep / (4.0 * M_PI) *
                    (1.0 / (x(r) - xr(c)) -
                    1.0 / (x(r) - xl(c)));
            if (std::abs(r-c) <= step){
                A(r,c) = matrix(r,c);
            }
            else{
                B(r,c) = matrix(r,c);
            }
        }
    }
    return;
}


const MatrixXd& Elasticity::getMatrix() const{
    return matrix;
}


void Elasticity::setCohesiveParameters(double dc, double dm, double Gc){
    this->dc = dc;
    this->dm = dm;
    this->Gc = Gc;
    sigmac = 2 * Gc / dc;
}