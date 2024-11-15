#include "MagmaState.hpp"
#include <exception>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using json = nlohmann::json;


MagmaState::MagmaState(Mesh* mesh, json&& properties, InputData* input) :
    mesh(mesh),
    properties(properties),
    input(input)
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
        auto sim_dir = input->getSimDir();
        std::string filepath;
        filepath = density_properties["dissolved_data_path"].get<std::string>();
        auto dissolved_data_path = sim_dir / "dissolved_data.json";
        std::filesystem::copy(filepath, dissolved_data_path, std::filesystem::copy_options::overwrite_existing);
        filepath = density_properties["gas_density_data_path"].get<std::string>();
        auto gas_density_data_path = sim_dir / "gas_density_data.json";
        std::filesystem::copy(filepath, gas_density_data_path, std::filesystem::copy_options::overwrite_existing);

        std::ifstream file;
        json data;
        file.open(dissolved_data_path, std::fstream::in);
        data = json::parse(file);
        auto pressure = data["pressure"].get<std::vector<double>>();
        auto wth2o = data["wth2o"].get<std::vector<double>>();
        auto wtco2 = data["wtco2"].get<std::vector<double>>();
        auto xh2og = data["xh2og"].get<std::vector<double>>();
        dissolved_weighted_h2o = std::make_unique<UniformInterpolation1d>(pressure, wth2o);
        dissolved_weighted_co2 = std::make_unique<UniformInterpolation1d>(pressure, wtco2);
        gas_h2o_co2_ratio = std::make_unique<UniformInterpolation1d>(pressure, xh2og);
        file.close();
        file.open(gas_density_data_path, std::fstream::in);
        data = json::parse(file);
        pressure = data["pressure"].get<std::vector<double>>();
        xh2og = data["xh2og"].get<std::vector<double>>();
        auto temperature = data["temperature"].get<std::vector<double>>();
        auto gas_density = data["density"].get<std::vector<double>>();
        gas_h2o_co2_density = std::make_unique<UniformInterpolation3d>(pressure, xh2og, temperature, gas_density);
        file.close();
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
                double beta = dike->beta(ix, iy);
                double wth2o = std::min(dissolved_weighted_h2o->getValue(p), chamber.wth2o);
                double wtco2 = std::min(dissolved_weighted_co2->getValue(p), chamber.wtco2);
                double gamma = wth2o + wtco2;
                double xh2og = gas_h2o_co2_ratio->getValue(p);
                dike->wth2o(ix, iy) = wth2o;
                dike->wtco2(ix, iy) = wtco2;
                dike->xh2od(ix, iy) = wth2o / gamma;
                dike->xh2og(ix, iy) = xh2og;
                dike->gamma(ix, iy) = gamma;

                double rhom_liquid = melt_density(rhom0, rhoh2o0, rhoco20, wth2o, wtco2);  // Effective melt density
                double rhog0 = gas_h2o_co2_density->getValue(p, xh2og, T);
                dike->rhom_liquid(ix, iy) = rhom_liquid;
                double t1 = (1-beta)*(chamber.Mg0-gamma)*rhom_liquid;
                double t2 = chamber.Mg0*beta*rhoc0;
                double alpha1 = (t1 + t2) / (t1 + t2 + (1.0-chamber.Mg0)*rhog0);
                double alpha = std::max(0.0, alpha1);
                dike->alpha(ix, iy) = alpha;

                double rhog = alpha*rhog0;
                double rhom = (1-alpha)*(1-beta)*rhom_liquid;
                double rhoc = (1.0-alpha)*beta*rhoc0;
                dike->rhog(ix, iy) = rhog;
                dike->rhom(ix, iy) = rhom;
                dike->rhoc(ix, iy) = rhoc;
                dike->density(ix, iy) = rhog + rhom + rhoc;
            }
        }
    }
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
            for (int iy = 0; iy < ny; iy++){
                double theta = theta_coef(beta(ix, iy));
                double mu_melt = grdvisc_model.calculateViscosity(gamma(ix, iy), T(ix, iy));
                viscosity(ix, iy) = std::min(theta*mu_melt, mu_max);
            }
        }
    }
}


void MagmaState::updateEquilibriumCrystallization(DikeData* dike) const{
    int nx = dike->meshX->size();
    int ny = dike->ny;
    for (int ix = 0; ix <= dike->tip_element; ix++){
        for (int iy = 0; iy < ny; iy++){
            double p = dike->pressure(ix);
            double T = dike->temperature(ix, iy);
            double xh2od = dike->xh2od(ix, iy);
            double Tl = liquidus_temperature(p, xh2od, sio2);
            double Ts = solidus_temperature(p, xh2od);
            double beta = beta_equilibrium(p, T, Tl, Ts);
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


std::vector<double> arrayxxd_to_rowwise(const ArrayXXd& arr){
    std::vector<double> result;
    result.reserve(arr.size());
    for (int i = 0; i < arr.rows(); i++){
        for (int j = 0; j < arr.cols(); j++){
            result.push_back(arr(i, j));
        }
    }
    return result;
}


void MagmaState::test() const{
    if (density_model == DensityModel::MIXED_H2O_CO2){
        json data;
        ArrayXd pressure = ArrayXd::LinSpaced(101, 0.0, 1000e6);
        ArrayXd temperature = ArrayXd::LinSpaced(51, 500.0, 1000.0);
        ArrayXd wth2o = pressure;
        ArrayXd wtco2 = pressure;
        ArrayXd xh2og = pressure;
        for (int i = 0; i < pressure.size(); i++){
            wth2o[i] = dissolved_weighted_h2o->getValue(pressure[i]);
            wtco2[i] = dissolved_weighted_co2->getValue(pressure[i]);
            xh2og[i] = gas_h2o_co2_ratio->getValue(pressure[i]);
        }
        data["dissolved"]["pressure"] = pressure;
        data["dissolved"]["wth2o"] = wth2o;
        data["dissolved"]["wtco2"] = wtco2;
        data["dissolved"]["xh2og"] = xh2og;

        ArrayXXd gas_density0 = ArrayXXd::Zero(pressure.size(), temperature.size());
        ArrayXXd gas_density05 = ArrayXXd::Zero(pressure.size(), temperature.size());
        ArrayXXd gas_density1 = ArrayXXd::Zero(pressure.size(), temperature.size());
        for (int i = 0; i < pressure.size(); i++){
            for (int j = 0; j < temperature.size(); j++){
                double p = pressure[i];
                double T = temperature[j];
                gas_density0(i, j) = gas_h2o_co2_density->getValue(p, 0.0, T);
                gas_density05(i, j) = gas_h2o_co2_density->getValue(p, 0.5, T);
                gas_density1(i, j) = gas_h2o_co2_density->getValue(p, 1.0, T);
            }
        }
        data["exsolved"]["pressure"] = pressure;
        data["exsolved"]["temperature"] = temperature;
        data["exsolved"]["case1"]["xh2og"] = 0.0;
        data["exsolved"]["case1"]["density"] = arrayxxd_to_rowwise(gas_density0);
        data["exsolved"]["case2"]["xh2og"] = 0.5;
        data["exsolved"]["case2"]["density"] = arrayxxd_to_rowwise(gas_density05);
        data["exsolved"]["case3"]["xh2og"] = 1.0;
        data["exsolved"]["case3"]["density"] = arrayxxd_to_rowwise(gas_density1);
        auto sim_dir = input->getSimDir();
        auto filepath = sim_dir / "test_magma_properties.json";
        std::ofstream f(filepath);
        f << data.dump(4) << std::endl;
        f.close();
    }
}