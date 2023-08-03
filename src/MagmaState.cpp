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
    else if(viscosity_model == "linear_viscosity"){
        double mu1 = viscosity_properties["much"];
        double mu2 = viscosity_properties["musurf"];
        dike->viscosity = VectorXd::LinSpaced(mesh->size(), mu1, mu2);
    }
}