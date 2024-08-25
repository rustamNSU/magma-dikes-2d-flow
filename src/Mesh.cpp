#include "Mesh.hpp"
using Eigen::ArrayXd;


Mesh::Mesh(const nlohmann::json& properties){
    n = properties["n"];
    xmin = properties["xmin"];
    xmax = properties["xmax"];
    dx = (xmax - xmin) / n;
    xl = ArrayXd::LinSpaced(n, xmin, xmax-dx);
    xr = ArrayXd::LinSpaced(n, xmin+dx, xmax);
    x = (xl + xr) / 2.0;
}


const ArrayXd& Mesh::getx() const{
    return x;
}

const ArrayXd& Mesh::getxl() const{
    return xl;
}

const ArrayXd& Mesh::getxr() const{
    return xr;
}


int Mesh::size() const{
    return n;
}


double Mesh::getdx() const{
    return dx;
}