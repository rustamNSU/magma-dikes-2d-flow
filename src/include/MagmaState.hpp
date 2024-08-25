#pragma once
#include <string>
#include <nlohmann/json.hpp>
#include "Mesh.hpp"
#include "DikeData.hpp"
#include "ViscosityModels.hpp"
#include "PinatubaCrystallization.hpp"
#include "Degasation.hpp"


class MagmaState{
    private:
        struct DensityModel{
            static constexpr int CONSTANT = 0;
            static constexpr int WATER_SATURATED = 1;
            inline static const std::string constant = "constant";
            inline static const std::string water_saturated = "water_saturated";
        };


        struct ViscosityModel{
            static constexpr int CONSTANT = 0;
            static constexpr int VFT_CONST_COEFF = 1;
            static constexpr int VFT_CONST_COEFF_AVG = 2;
            static constexpr int VFT_CONST_COEFF_CRYST = 3;
            static constexpr int GRDMODEL08 = 4;

            inline static const std::string constant = "constant";
            inline static const std::string vft_const_coeff = "vft_const_coeff";
            inline static const std::string vft_const_coeff_avg = "vft_const_coeff_avg";
            inline static const std::string vft_const_coeff_cryst = "vft_const_coeff_cryst";
            inline static const std::string grdmodel08 = "grdmodel08";
        };

        struct SaturationModel{
            static constexpr int LAVALLEE2015 = 0;
            inline static const std::string lavallee2015 = "lavallee2015";
        };


    private:
        nlohmann::json properties;
        int density_model = 0;
        int viscosity_model = 0;
        int saturation_model = 0;
        nlohmann::json density_properties;
        nlohmann::json viscosity_properties;
        nlohmann::json saturation_properties;
        double thermal_conductivity;
        double specific_heat;
        double latent_heat;
        double gamma0 = 0.0;
        double beta0 = 0.0;
        double rho0 = 0.0;
        double rhom0 = 0.0;
        double Mg0 = 0.0; // m_g+m_d / m_tot in chamber
        Mesh* mesh;
    
    public:
        MagmaState(Mesh* mesh, nlohmann::json&& properties);
        void setDensityModel();
        void setViscosityModel();
        void setCrystallizationModel();
        void setSaturationModel();
        void updateGasDensity(DikeData* dike) const;
        void updateMeltLiquidDensity(DikeData* dike) const;
        void updateDensity(DikeData* dike) const;
        void updateViscosity(DikeData* dike) const;
        void updateEquilibriumCrystallization(DikeData* dike) const;
        void updateGasSaturation(DikeData* dike) const;
        double getThermalConductivity() const;
        double getSpecificHeat() const;
        inline double getLatentHeat() const{
            return latent_heat;
        }
        void setChamberInitialState(DikeData* dike);
};