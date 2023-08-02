#include "MagmaState.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using json = nlohmann::json;


MagmaState::MagmaState(const json& input, Mesh* mesh) : mesh(mesh) {
    auto magma_properties = input["magmaProperties"];
    density_model = std::string(magma_properties["densityModel"]);
    viscosity_model = std::string(magma_properties["viscosityModel"]);

    if (density_model == "constant_density"){
        rho = magma_properties["constantDensity"]["rho"];
    }
    if (viscosity_model == "constant_viscosity"){
        mu = magma_properties["constantViscosity"]["mu"];
    }
}


void MagmaState::updateDensity(DikeData* dike) const{
    if (density_model == "constant_density"){
        dike->density.fill(rho);
    }
}


void MagmaState::updateViscosity(DikeData* dike) const{
    if (viscosity_model == "constant_viscosity"){
        dike->viscosity.fill(mu);
    }
}