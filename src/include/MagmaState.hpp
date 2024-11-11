#pragma once
#include <string>
#include <memory>
#include <nlohmann/json.hpp>
#include "Mesh.hpp"
#include "DikeData.hpp"
#include "ViscosityModels.hpp"
#include "PinatubaCrystallization.hpp"
#include "Degasation.hpp"
#include "uniform_interp.hpp"
#include "InputData.hpp"

class DikeModel2d;


class MagmaState{
    private:
        struct DensityModel{
            static constexpr int CONSTANT = 0;
            static constexpr int MIXED_H2O_CO2 = 1;
            inline static const std::string constant = "constant";
            inline static const std::string mixed_h2o_co2 = "mixed_h2o_co2";
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
            static constexpr int MIXED_H2O_CO2 = 0;
            inline static const std::string mixed_h2o_co2 = "mixed_h2o_co2";
        };


        struct Chamber{
            double pressure;
            double temperature;
            double alpha;
            double gamma;
            double wth2o; // dissolved h2o wt.%
            double wtco2; // dissolved co2 wt.%
            double xh2og; // gas h2o / (h2o + co2)
            double xh2od; // dissolved gas h2o / (h2o + co2)
            double beta;
            double density;
            double rhom_liquid;
            double rhom;
            double rhoc;
            double Mg0;
            double Tl;
            double Ts;
            inline void calculateGasRatio(){
                Mg0 = (1.0 - beta) * gamma * rhom_liquid / density;
            }
        };

    private:
        InputData* input;
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
        Chamber chamber;
        Mesh* mesh;
        double sio2 = 63.888;
        GiordanoViscosity grdvisc_model;
        std::unique_ptr<UniformInterpolation1d> dissolved_weighted_h2o;
        std::unique_ptr<UniformInterpolation1d> dissolved_weighted_co2;
        std::unique_ptr<UniformInterpolation1d> gas_h2o_co2_ratio;
        std::unique_ptr<UniformInterpolation3d> gas_h2o_co2_density; // [pressure, xh2od, temperature]
    
    public:
        MagmaState(Mesh* mesh, nlohmann::json&& properties, InputData* input);
        void setDensityModel();
        void setViscosityModel();
        void setCrystallizationModel();
        void setSaturationModel();
        void updateDensity(DikeData* dike) const;
        void updateViscosity(DikeData* dike) const;
        void updateEquilibriumCrystallization(DikeData* dike) const;
        double getThermalConductivity() const;
        double getSpecificHeat() const;
        inline double getLatentHeat() const{
            return latent_heat;
        }
        void setChamberInitialState(
            double pressure_chamber,
            double temperature_chamber
        );
        double getMagmaChamberCrystallization() const{
            return chamber.beta;
        }

        friend class DikeModel2d;
};