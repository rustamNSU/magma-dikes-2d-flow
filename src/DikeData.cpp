#include "DikeData.hpp"
#include <cmath>
#include <highfive/H5Easy.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using nlohmann::json;
using H5Easy::File;
using H5Easy::dump;


DikeData::DikeData(Mesh* mesh, const json& alg_properties) :
    meshX(mesh),
    algorithm_properties(alg_properties)
{
    ny = algorithm_properties["numberOfLayers"];
    int nx = meshX->size();
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


void DikeData::save(const std::string path) const{
    File file(path, File::Overwrite);
    dump(file, "mesh/xc", meshX->getx());
    dump(file, "mesh/xl", meshX->getxl());
    dump(file, "mesh/xr", meshX->getxr());

    dump(file, "ny", ny);
    dump(file, "yc", yc);
    dump(file, "yb", yb);
    dump(file, "halfwidth", hw);
    dump(file, "width", getWidth());
    dump(file, "pressure", pressure);
    dump(file, "overpressure", overpressure);
    dump(file, "density", density);
    dump(file, "viscosity", viscosity);
    dump(file, "temperature", temperature);
    dump(file, "Twall", Twall);
    dump(file, "qx", qx);
    dump(file, "qy", qy);
    dump(file, "time", time);
}