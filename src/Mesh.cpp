#include "Mesh.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;


Mesh::Mesh(const nlohmann::json& properties){
    n = properties["n"];
    xmin = properties["xmin"];
    xmax = properties["xmax"];
    dx = (xmax - xmin) / n;
    xl = VectorXd::LinSpaced(n, xmin, xmax-dx);
    xr = VectorXd::LinSpaced(n, xmin+dx, xmax);
    x = (xl + xr) / 2.0;
}


const VectorXd& Mesh::getx() const{
    return x;
}

const VectorXd& Mesh::getxl() const{
    return xl;
}

const VectorXd& Mesh::getxr() const{
    return xr;
}


int Mesh::size() const{
    return n;
}


double Mesh::getdx() const{
    return dx;
}