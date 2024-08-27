#include "DikeData.hpp"
#include <cmath>
#include <highfive/H5Easy.hpp>

using Eigen::ArrayXXd;
using Eigen::ArrayXd;
using nlohmann::json;
using H5Easy::File;
using H5Easy::dump;


DikeData::DikeData(Mesh* mesh, const json& alg_properties) :
    meshX(mesh),
    algorithm_properties(alg_properties)
{
    ny = algorithm_properties["numberOfLayers"];
    int nx = meshX->size();
    hw = ArrayXd::Zero(nx);
    pressure = ArrayXd::Zero(nx);
    overpressure = ArrayXd::Zero(nx);

    yb = ArrayXd::LinSpaced(ny+1, 0.0, 1.0);;
    yc = 0.5 * (
        yb(Eigen::seq(0, Eigen::last - 1)) + 
        yb(Eigen::seq(1, Eigen::last))
    );
    
    density = ArrayXXd::Zero(nx, ny);
    rhom = ArrayXXd::Zero(nx, ny);
    rhoc = ArrayXXd::Zero(nx, ny);
    rhog = ArrayXXd::Zero(nx, ny);
    rhom_liquid = ArrayXXd::Zero(nx, ny);
    temperature = ArrayXXd::Zero(nx, ny);
    viscosity = ArrayXXd::Zero(nx, ny);
    Tliquidus = ArrayXXd::Zero(nx, ny);
    Tsolidus = ArrayXXd::Zero(nx, ny);
    betaeq = ArrayXXd::Zero(nx, ny);
    beta = ArrayXXd::Zero(nx, ny);
    alpha = ArrayXXd::Zero(nx, ny);
    gamma = ArrayXXd::Zero(nx, ny);
    qx = ArrayXXd::Zero(nx+1, ny);
    qy = ArrayXXd::Zero(nx, ny+1);
    A = ArrayXXd::Zero(nx+1, ny);
    C = ArrayXXd::Zero(nx+1, ny);
    shear_heat = ArrayXXd::Zero(nx, ny);
    mobility = ArrayXXd::Zero(nx, ny);
    Qx = ArrayXd::Zero(nx+1);
    Mx = ArrayXd::Zero(nx+1);
    Twall = ArrayXd::Zero(nx);
    G = ArrayXd::Zero(nx+1);
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
    dump(file, "temperature", temperature);
    dump(file, "viscosity", viscosity);
    dump(file, "betaeq", betaeq);
    dump(file, "beta", beta);
    dump(file, "alpha", alpha);
    dump(file, "gamma", gamma);
    dump(file, "Twall", Twall);
    dump(file, "dpdz", G);
    dump(file, "qx", qx);
    dump(file, "qy", qy);
    dump(file, "time", time);
}