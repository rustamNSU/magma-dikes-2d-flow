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
    ny = algorithm_properties["number_of_layers"].get<int>();
    MIN_WIDTH = algorithm_properties["min_width"].get<double>();
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
    tau = ArrayXXd::Zero(nx, ny);
    alpha = ArrayXXd::Zero(nx, ny);
    gamma = ArrayXXd::Zero(nx, ny);
    wth2o = ArrayXXd::Zero(nx, ny);
    wtco2 = ArrayXXd::Zero(nx, ny);
    xh2od = ArrayXXd::Zero(nx, ny);
    xh2og = ArrayXXd::Zero(nx, ny);
    qx = ArrayXXd::Zero(nx+1, ny);
    qy = ArrayXXd::Zero(nx, ny+1);
    mx = ArrayXXd::Zero(nx+1, ny);
    my = ArrayXXd::Zero(nx, ny+1);
    A = ArrayXXd::Zero(nx+1, ny);
    C = ArrayXXd::Zero(nx+1, ny);
    shear_heat = ArrayXXd::Zero(nx, ny);
    mobility = ArrayXXd::Zero(nx, ny);
    Qx = ArrayXd::Zero(nx+1);
    Mx = ArrayXd::Zero(nx+1);
    Twall = ArrayXd::Zero(nx);
    G = ArrayXd::Zero(nx+1);
    magma_to_rock_heat_flux = ArrayXd::Zero(nx);
    open_elements = std::vector<bool>(nx, false);
    tip_element = 0;
    front = meshX->getxr()[tip_element];
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
    dump(file, "mx", mx);
    dump(file, "my", my);
    dump(file, "time", time);
}


void DikeData::updateOpenElements(){
    tip_element = 0;
    for (int i = 0; i < meshX->size(); i++){
        if (hw[i] > MIN_WIDTH){
            open_elements[i] = true;
            tip_element = i;
        }
        else{
            open_elements[i] = false;
        }
    }
    front = meshX->getxr()(tip_element);
}


void DikeData::setMagmaStateAfterTip(){
    setMagmaStateAfterTip(tip_element, std::min(tip_element + 1, meshX->size()-1));
    setMagmaStateAfterTip(tip_element, std::min(tip_element + 2, meshX->size()-1));
    setMagmaStateAfterTip(tip_element, std::min(tip_element + 3, meshX->size()-1));
}


void DikeData::setMagmaStateAfterTip(int ntip, int nout){
    if (ntip == nout) return;
    density.row(nout) = density.row(ntip);
    viscosity.row(nout) = viscosity.row(ntip);
    beta.row(nout) = beta.row(ntip);
    betaeq.row(nout) = betaeq.row(ntip);
    tau.row(nout) = tau.row(ntip);
    alpha.row(nout) = alpha.row(ntip);
    gamma.row(nout) = gamma.row(ntip);
    wth2o.row(nout) = wth2o.row(ntip);
    wtco2.row(nout) = wtco2.row(ntip);
    xh2od.row(nout) = xh2od.row(ntip);
    xh2og.row(nout) = xh2og.row(ntip);
    temperature.row(nout) = temperature.row(ntip);
    return;
}


double DikeData::getTotalMass() const{
    return 2.0 * getElementsVolume().cwiseProduct(density.rowwise().mean()).sum();
}


ArrayXd DikeData::getElementsHalfMass() const{
    return getElementsVolume().cwiseProduct(density.rowwise().mean());
}


ArrayXd DikeData::getElementsHalfEnergy(double Cm) const{
    ArrayXXd rhoT = density * temperature;
    ArrayXd energy = Cm * getElementsVolume().cwiseProduct(rhoT.rowwise().mean());
    return energy;
}


double DikeData::calculateMassBalanceError(const Eigen::ArrayXd &Mold, double dt) const{
    double error = 0.0;
    ArrayXd Mnew = getElementsHalfMass();
    for (int i = 0; i < meshX->size(); i++){
        double res = std::abs(Mnew[i] - Mold[i] - dt * Mx[i] + dt * Mx[i+1]);
        double err = res / std::max(MIN_WIDTH * 2000 * meshX->getdx(), Mnew[i]);
        error = error > err ? error : err;
    }
    return error;
}