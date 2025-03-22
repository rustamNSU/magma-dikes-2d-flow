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
            inline static const std::string constant = "constant_density";
            inline static const std::string mixed_h2o_co2 = "mixed_h2o_co2";
        };

        struct ViscosityModel{
            inline static const std::string constant = "constant_viscosity";
            inline static const std::string vft_const_coeff = "vft_const_coeff";
            inline static const std::string grdmodel08 = "grdmodel08";
        };

        struct CrystallizationModel{
            inline static const std::string equilibrium_crystallization = "equilibrium_crystallization";
            inline static const std::string constant_relaxation_crystallization = "constant_relaxation_crystallization";
            inline static const std::string arrhenius_relaxation_crystallization = "arrhenius_relaxation_crystallization";
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
        nlohmann::json density_model;
        nlohmann::json viscosity_model;
        nlohmann::json crystallization_model;
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
        void updateRelaxationCrystallization(DikeData* dike) const;
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

        void test() const;

        friend class DikeModel2d;
};