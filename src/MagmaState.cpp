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
    saturation_model = SaturationModel::MIXED_H2O_CO2;
}


void MagmaState::setDensityModel(){
    auto model = properties["densityModel"].get<std::string>();
    if (model == DensityModel::constant){
        density_model = DensityModel::CONSTANT;
        density_properties = properties["constant"];
    }
    else if (model == DensityModel::mixed_h2o_co2){
        density_model = DensityModel::MIXED_H2O_CO2;
        density_properties = properties["mixed_h2o_co2"];
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
    else if (model == ViscosityModel::grdmodel08){
        viscosity_model = ViscosityModel::GRDMODEL08;
        viscosity_properties = properties[ViscosityModel::grdmodel08];
        auto composition = viscosity_properties["composition"].get<std::array<double, 11>>();
        sio2 = composition[0];
        grdvisc_model.setComposition(composition);
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
    else if (density_model == DensityModel::MIXED_H2O_CO2){
        updateMeltDensity(dike);
        double rhoc0 = density_properties["rhoc0"].get<double>();
        dike->rhoc = rhoc0 * (1.0 - dike->alpha) * dike->beta;
        dike->rhom = (1.0 - dike->alpha) * (1.0 - dike->beta) * dike->rhom_liquid;
        dike->density = dike->rhog + dike->rhoc + dike->rhom;
    }
}


void MagmaState::updateMeltDensity(DikeData* dike) const{
    int nx = dike->meshX->size();
    int ny = dike->getLayersNumber();
    double rhom0 = density_properties["rhom0"].get<double>();
    double rhoh2o0 = density_properties["rhoh2o0"].get<double>();
    double rhoco20 = density_properties["rhoco20"].get<double>();
    double rhoc0 = density_properties["rhoc0"].get<double>();
    for (int ix = 0; ix <= dike->tip_element; ix++){
        for (int iy = 0; iy < ny; iy++){
            double p = dike->pressure(ix);
            double T = dike->temperature(ix, iy);
            double wth2o = std::min(dissolved_weighted_h2o->getValue(p), chamber.wth2o);
            double wtco2 = std::min(dissolved_weighted_co2->getValue(p), chamber.wtco2);
            double gamma = wth2o + wtco2;
            double xh2og = gas_h2o_co2_ratio->getValue(p);
            dike->wth2o(ix, iy) = wth2o;
            dike->wtco2(ix, iy) = wtco2;
            dike->xh2od(ix, iy) = wth2o / gamma;
            dike->xh2og(ix, iy) = xh2og;
            dike->gamma(ix, iy) = gamma;

            double rhom_liquid = melt_density(rhom0, rhoh2o0, rhoco20, wth2o, wtco2);
            double rhog0 = gas_h2o_co2_density->getValue(p, xh2og, T);
            double t1 = (1-beta(ix, iy))*(chamber.Mg0-gamma(ix,iy))*rhom_liquid(ix,iy);
            double t2 = chamber.Mg0*beta(ix,iy)*rhoc0;
            double alpha1 = (t1 + t2) / (t1 + t2 + (1.0-chamber.Mg0)*rhog0);
            double alpha = std::max(0.0, alpha1);
            dike->alpha(ix, iy) = alpha;
            dike->rhog(ix, iy) = rhog0 * alpha;
        }
    }
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
    else if (viscosity_model == ViscosityModel::GRDMODEL08){
        int nx = mesh->size();
        int ny = dike->getLayersNumber();
        double mu_max = viscosity_properties["muMaxLimit"].get<double>();
        double K = 273.15;
        const auto& T = dike->temperature;
        const auto& beta = dike->beta;
        const auto& gamma = dike->gamma;
        auto& viscosity = dike->viscosity;
        const auto& hw = dike->hw;
        for (int ix = 0; ix <= dike->tip_element; ix++){
            int il = std::max(0, ix-1);
            int ir = std::min(nx-1, ix+1);
            if (std::max({hw[il], hw[ix], hw[ir]}) < 1e-13 && ix > std::min(10, nx)){
                viscosity.row(ix).fill(mu_max);
            }
            else{
                for (int iy = 0; iy < ny; iy++){
                    double theta = theta_coef(beta(ix, iy));
                    // double mu_melt = T(ix, iy) + K > C ? VFT_constant_viscosity(T(ix, iy) + K, A, B, C) : mu_max;
                    double mu_melt = grdvisc_model.calculateViscosity(gamma(ix, iy), T(ix, iy));
                    viscosity(ix, iy) = std::min(theta*mu_melt, mu_max);
                }
            }
        }
    }
}


void MagmaState::updateEquilibriumCrystallization(DikeData* dike) const{
    int nx = dike->meshX->size();
    int ny = dike->ny;
    // #pragma omp parallel for
    const auto& hw = dike->hw;
    for (int ix = 0; ix <= dike->tip_element; ix++){
        int il = std::max(0, ix-1);
        int ir = std::min(nx-1, ix+1);
        for (int iy = 0; iy < ny; iy++){
            double Tl = liquidus_temperature(dike->pressure(ix), sio2);
            double Ts = solidus_temperature(dike->pressure(ix));
            double beta = beta_equilibrium(dike->pressure(ix), dike->temperature(ix, iy), Tl, Ts);
            dike->betaeq(ix, iy) = std::max(beta, chamber.beta);
            dike->Tliquidus(ix, iy) = Tl;
            dike->Tsolidus(ix, iy) = Ts;
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


void MagmaState::setChamberInitialState(
    double pressure_chamber,
    double temperature_chamber
){
    chamber.pressure = pressure_chamber;
    chamber.temperature = temperature_chamber;
    chamber.alpha = 0.0;
    if (density_model == DensityModel::MIXED_H2O_CO2){
        chamber.wth2o = dissolved_weighted_h2o->getValue(chamber.pressure);
        chamber.wtco2 = dissolved_weighted_co2->getValue(chamber.pressure);
        chamber.gamma = chamber.wth2o + chamber.wtco2;
        chamber.xh2od = chamber.wth2o / chamber.gamma;
        chamber.xh2og = gas_h2o_co2_ratio->getValue(chamber.pressure);
        chamber.Tl = liquidus_temperature(chamber.pressure, chamber.xh2od, sio2);
        chamber.Ts = solidus_temperature(chamber.pressure, chamber.xh2od);
        chamber.beta = beta_equilibrium(chamber.pressure, chamber.temperature, chamber.Tl, chamber.Ts);
        double rhom0 = density_properties["rhom0"].get<double>();
        double rhoh2o0 = density_properties["rhoh2o0"].get<double>();
        double rhoco20 = density_properties["rhoco20"].get<double>();
        double rhoc0 = density_properties["rhoc0"].get<double>();
        chamber.rhom_liquid = melt_density(rhom0, rhoh2o0, rhoco20, chamber.wth2o, chamber.wtco2);
        chamber.rhom = (1.0 - chamber.beta) * chamber.rhom_liquid;
        chamber.rhoc = rhoc0 * chamber.beta;
        chamber.density = chamber.rhom + chamber.rhoc;
        chamber.calculateGasRatio();
    }
    return;
}