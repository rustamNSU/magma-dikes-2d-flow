#include "MagmaState.hpp"
#include <exception>
#ifdef USE_OMP
#include <omp.h>
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using json = nlohmann::json;


MagmaState::MagmaState(Mesh* mesh, json&& properties) : 
    mesh(mesh),
    properties(properties)
{
    setDensityModel();
    setViscosityModel();
    thermal_conductivity = this->properties["thermalConductivity"];
    specific_heat = this->properties["specificHeatCapacity"];
    latent_heat = this->properties["latentHeat"];
    saturation_model = SaturationModel::LAVALLEE2015;
}


void MagmaState::setDensityModel(){
    auto model = properties["densityModel"].get<std::string>();
    if (model == DensityModel::constant){
        density_model = DensityModel::CONSTANT;
        density_properties = properties["constantDensity"];
    }
    else if (model == DensityModel::water_saturated){
        density_model = DensityModel::WATER_SATURATED;
        density_properties = properties["waterSaturatedDensity"];
    }
    else{
        throw std::invalid_argument("Density model is incorrect!\n");
    }
    return;
}


void MagmaState::setViscosityModel(){
    auto model = properties["viscosityModel"].get<std::string>();
    if (model == ViscosityModel::constant){
        viscosity_model = ViscosityModel::CONSTANT;
        viscosity_properties = properties["constantViscosity"];
    }
    else if (model == ViscosityModel::vft_const_coeff){
        viscosity_model = ViscosityModel::VFT_CONST_COEFF;
        viscosity_properties = properties["vftConstantViscosity"];
    }
    else if (model == ViscosityModel::vft_const_coeff_avg){
        viscosity_model = ViscosityModel::VFT_CONST_COEFF_AVG;
        viscosity_properties = properties["vftConstantViscosity"];
    }
    else if (model == ViscosityModel::vft_const_coeff_cryst){
        viscosity_model = ViscosityModel::VFT_CONST_COEFF_CRYST;
        viscosity_properties = properties["vftConstantViscosity"];
    }
    else{
        throw std::invalid_argument("Viscosity model is incorrect (ro not realized)!\n");
    }
    return;
}


void MagmaState::updateDensity(DikeData* dike) const{
    if (density_model == DensityModel::CONSTANT){
        double rho = density_properties["rho"].get<double>();
        dike->density.fill(rho);
    }
    else if (density_model == DensityModel::WATER_SATURATED){
        updateGasSaturation(dike);
        updateMeltLiquidDensity(dike);
        if (!INITIAL_DATA){
            updateGasDensity(dike);
        }
        double rhom0 = density_properties["rhom0"].get<double>();
        double rhow0 = density_properties["rhow0"].get<double>();
        double rhoc0 = density_properties["rhoc0"].get<double>();
        dike->rhoc = rhoc0 * (1.0 - dike->alpha) * dike->beta;
        dike->rhom = (1.0 - dike->alpha) * (1.0 - dike->beta) * dike->rhom_liquid;
        dike->density = dike->rhog + dike->rhoc + dike->rhom;

        // @debug
        dike->density.fill(rhom0);
    }
}


void MagmaState::updateGasDensity(DikeData* dike) const{
    int nx = dike->meshX->size();
    int ny = dike->getLayersNumber();
    double rhoc0 = density_properties["rhoc0"].get<double>();
    const auto& hw = dike->hw;
    const auto& p = dike->pressure;
    const auto& T = dike->temperature;
    const auto& beta = dike->beta;
    const auto& gamma = dike->gamma;
    const auto& rhom_liquid = dike->rhom_liquid;
    double R = 1000.0;
    for (int ix = 0; ix < nx; ix++){
        int il = std::max(0, ix-1);
        int ir = std::min(nx-1, ix+1);
        if (std::max({hw[il], hw[ix], hw[ir]}) < 1e-13 && ix > std::min(10, nx)){
            continue;
        }
        else{
            for (int iy = 0; iy < ny; iy++){
                double rhog0 = h2o_vapor_density(p(ix), T(ix, iy), R);
                double t1 = (1-beta(ix, iy))*(Mg0-gamma(ix,iy))*rhom_liquid(ix,iy);
                double t2 = Mg0*beta(ix,iy)*rhoc0;
                double alpha = std::max(0.0, (t1 + t2) / (t1 + t2 + (1.0-Mg0)*rhog0));
                dike->rhog(ix, iy) = rhog0 * dike->alpha(ix, iy);
            }
        }
    }
    return;
}



void MagmaState::updateMeltLiquidDensity(DikeData* dike) const{
    int nx = dike->meshX->size();
    int ny = dike->getLayersNumber();
    const auto& hw = dike->hw;
    const auto& gamma = dike->gamma;
    double rhom0 = density_properties["rhom0"].get<double>();
    double rhow0 = density_properties["rhow0"].get<double>();
    auto func = [=](double gamma){
        return h2o_sat_melt_density(rhom0, rhow0, gamma);
    };
    dike->rhom_liquid = dike->gamma.unaryExpr(func);
    return;
}


void MagmaState::updateViscosity(DikeData* dike) const{
    if (viscosity_model == ViscosityModel::VFT_CONST_COEFF_CRYST){
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
    else if (viscosity_model == ViscosityModel::CONSTANT){
        dike->viscosity.fill(viscosity_properties["mu"]);
    }
    else if (viscosity_model == ViscosityModel::VFT_CONST_COEFF){
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
    else if (viscosity_model == ViscosityModel::VFT_CONST_COEFF_AVG){
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


void MagmaState::updateGasSaturation(DikeData* dike) const{
    int nx = dike->meshX->size();
    int ny = dike->ny;
    const auto& hw = dike->hw;
    const auto& p = dike->pressure;
    const auto& T = dike->temperature;
    const auto& gamma = dike->gamma;
    for (int ix = 0; ix < nx; ix++){
        int il = std::max(0, ix-1);
        int ir = std::min(nx-1, ix+1);
        if (std::max({hw[il], hw[ix], hw[ir]}) < 1e-13 && ix > std::min(10, nx)) continue;
        for (int iy = 0; iy < ny; iy++){
            dike->gamma(ix, iy) = h2o_wt_lavallee2015(p(ix), T(ix, iy), 0.5);
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


void MagmaState::setChamberInitialState(DikeData* dike){
    gamma0 = dike->gamma(0, 0);
    beta0 = dike->betaeq(0, 0);
    rho0 = dike->density(0, 0);
    rhom0 = dike->rhom(0, 0);
    Mg0 = (1.0 - beta0) * gamma0 * rhom0 / rho0;
    INITIAL_DATA = false;
}