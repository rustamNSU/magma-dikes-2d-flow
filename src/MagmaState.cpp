#include "MagmaState.hpp"
#ifdef USE_OMP
#include <omp.h>
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using json = nlohmann::json;


MagmaState::MagmaState(Mesh* mesh, json&& properties) : 
    mesh(mesh),
    properties(properties)
{
    density_model = this->properties["densityModel"].get<std::string>();
    viscosity_model = this->properties["viscosityModel"].get<std::string>();
    thermal_conductivity = this->properties["thermalConductivity"];
    specific_heat = this->properties["specificHeatCapacity"];
    latent_heat = this->properties["latentHeat"];

    if (density_model == "constantDensity"){
        rho = this->properties["constantDensity"]["rho"];
    }
    if (viscosity_model == "constantViscosity"){
        viscosity_properties = this->properties["constantViscosity"];
    }
    else if (viscosity_model == "linearViscosity"){
        viscosity_properties = this->properties["linearViscosity"];
    }
    else if (viscosity_model == "vftConstantViscosity" || "averageVftConstantViscosity" || "vftConstantViscosityCrystallization"){
        viscosity_properties = this->properties["vftConstantViscosity"];
    }
}


void MagmaState::updateDensity(DikeData* dike) const{
    if (density_model == "constantDensity"){
        dike->density.fill(rho);
    }
}


void MagmaState::updateViscosity(DikeData* dike) const{
    if (viscosity_model == "vftConstantViscosityCrystallization"){
        int nx = mesh->size();
        int ny = dike->getLayersNumber();
        const double A = viscosity_properties["A"].get<double>();
        const double B = viscosity_properties["B"].get<double>();
        const double C = viscosity_properties["C"].get<double>();
        double mu_max = viscosity_properties["muMaxLimit"].get<double>();
        double K = 273.15;
        const auto& T = dike->temperature;
        const auto& beta = dike->beta;
        auto& viscosity = dike->viscosity;
        // #pragma omp parallel for
        const auto& hw = dike->hw;
        for (int ix = 0; ix < nx; ix++){
            int il = std::max(0, ix-1);
            int ir = std::min(nx-1, ix+1);
            if (std::max({hw[il], hw[ix], hw[ir]}) < 1e-13 && ix > std::min(10, nx)){
                viscosity.row(ix).fill(mu_max);
            }
            else{
                for (int iy = 0; iy < ny; iy++){
                    double theta = theta_coef(beta(ix, iy));
                    // double mu_melt = T(ix, iy) + K > C ? VFT_constant_viscosity(T(ix, iy) + K, A, B, C) : mu_max;
                    double mu_melt = VFT_constant_viscosity(T(ix, iy) + K, A, B, C);
                    viscosity(ix, iy) = std::min(theta*mu_melt, mu_max);
                }
            }
        }
        // viscosity = T.binaryExpr()
    }
    else if (viscosity_model == "constantViscosity"){
        dike->viscosity.fill(viscosity_properties["mu"]);
    }
    else if (viscosity_model == "linearViscosity"){
        double mu1 = viscosity_properties["much"];
        double mu2 = viscosity_properties["musurf"];
        auto mu_col = VectorXd::LinSpaced(mesh->size(), mu1, mu2);
        for (int icol = 0; icol < dike->getLayersNumber(); icol++){
            dike->viscosity.col(icol) = mu_col;
        }
    }
    else if (viscosity_model == "vftConstantViscosity"){
        int nx = dike->meshX->size();
        int ny = dike->ny;
        double A = viscosity_properties["A"];
        double B = viscosity_properties["B"];
        double C = viscosity_properties["C"];
        double mu_max = viscosity_properties["muMaxLimit"];
        double K = 273.15;
        auto func = [=](double T){
            double mu = T + K > C ? VFT_constant_viscosity(T + K, A, B, C) : mu_max; 
            return std::min(mu, mu_max);
        };
        dike->viscosity = dike->temperature.unaryExpr(func);
    }
    else if (viscosity_model == "averageVftConstantViscosity"){
        int nx = mesh->size();
        int ny = dike->getLayersNumber();
        double A = viscosity_properties["A"];
        double B = viscosity_properties["B"];
        double C = viscosity_properties["C"];
        double mu_max = viscosity_properties["muMaxLimit"];
        double K = 273.15;
        auto func = [=](double T){
            double mu = T + K > C ? VFT_constant_viscosity(T + K, A, B, C) : mu_max; 
            return std::min(mu, mu_max);
        };
        MatrixXd Tavg = dike->temperature;
        for (int ix = 0; ix < nx; ix++){
            Tavg.row(ix).fill(dike->temperature.row(ix).mean());
        }
        dike->viscosity = Tavg.unaryExpr(func);
    }
}


void MagmaState::updateEquilibriumCrystallization(DikeData* dike) const{
    int nx = dike->meshX->size();
    int ny = dike->ny;
    // #pragma omp parallel for
    const auto& hw = dike->hw;
    for (int ix = 0; ix < nx; ix++){
        int il = std::max(0, ix-1);
        int ir = std::min(nx-1, ix+1);
        if (std::max({hw[il], hw[ix], hw[ir]}) < 1e-13 && ix > std::min(10, nx)) continue;
        for (int iy = 0; iy < ny; iy++){
            dike->betaeq(ix, iy) = beta_equilibrium(dike->pressure(ix), dike->temperature(ix, iy));
            dike->Tliquidus(ix, iy) = liquidus_temperature(dike->pressure(ix));
            dike->Tsolidus(ix, iy) = solidus_temperature(dike->pressure(ix));
        }
    }
    return;
}


double MagmaState::getThermalConductivity() const{
    return thermal_conductivity;
}


double MagmaState::getSpecificHeat() const{
    return specific_heat;
}