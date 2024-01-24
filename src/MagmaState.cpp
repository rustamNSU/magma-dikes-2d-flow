#include "MagmaState.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using json = nlohmann::json;


MagmaState::MagmaState(const json& magma_properties, Mesh* mesh) : 
    magma_properties(magma_properties),
    mesh(mesh)
{
    density_model = std::string(magma_properties["densityModel"]);
    viscosity_model = std::string(magma_properties["viscosityModel"]);
    thermal_conductivity = magma_properties["thermalConductivity"];
    specific_heat = magma_properties["specificHeat"];

    if (density_model == "constant_density"){
        rho = magma_properties["constantDensity"]["rho"];
    }
    if (viscosity_model == "constant_viscosity"){
        viscosity_properties = magma_properties["constantViscosity"];
    }
    else if (viscosity_model == "linear_viscosity"){
        viscosity_properties = magma_properties["linearViscosity"];
    }
}


void MagmaState::updateDensity(DikeData* dike) const{
    if (density_model == "constant_density"){
        dike->density.fill(rho);
    }
}


void MagmaState::updateViscosity(DikeData* dike) const{
    if (viscosity_model == "constant_viscosity"){
        dike->viscosity.fill(viscosity_properties["mu"]);
    }
    else if (viscosity_model == "linear_viscosity"){
        double mu1 = viscosity_properties["much"];
        double mu2 = viscosity_properties["musurf"];
        auto mu_col = VectorXd::LinSpaced(mesh->size(), mu1, mu2);
        for (int icol = 0; icol < dike->getLayersNumber(); icol++){
            dike->viscosity.col(icol) = mu_col;
        }
    }
    else if (viscosity_model == "VFT_constant_viscosity"){
        int nx = mesh->size();
        int ny = dike->getLayersNumber();
        double A = viscosity_properties["A"];
        double B = viscosity_properties["B"];
        double C = viscosity_properties["C"];
        double mu_max = viscosity_properties["max_viscosity"];
        auto func = [=](double T){
            double mu = T > C ? VFT_constant_viscosity(T, A, B, C) : mu_max; 
            return std::min(mu, mu_max);
        };
        dike->viscosity = dike->temperature.unaryExpr(func);
    }
}


double MagmaState::getThermalConductivity() const{
    return thermal_conductivity;
}


double MagmaState::getSpecificHeat() const{
    return specific_heat;
}