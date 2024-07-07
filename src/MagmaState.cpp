#include "MagmaState.hpp"

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
    else if (viscosity_model == "vftConstantViscosity" || "averageVftConstantViscosity"){
        viscosity_properties = this->properties["vftConstantViscosity"];
    }
}


void MagmaState::updateDensity(DikeData* dike) const{
    if (density_model == "constantDensity"){
        dike->density.fill(rho);
    }
}


void MagmaState::updateViscosity(DikeData* dike) const{
    if (viscosity_model == "constantViscosity"){
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
    for (int ix = 0; ix < nx; ix++){
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