#pragma once
#include <string>
#include <nlohmann/json.hpp>
#include "Mesh.hpp"
#include "DikeData.hpp"
#include "ViscosityModels.hpp"
#include "PinatubaCrystallization.hpp"


class MagmaState{
    private:
        nlohmann::json properties;
        std::string density_model = "constant_density";
        std::string viscosity_model = "constant_viscosity";
        double rho;
        double thermal_conductivity;
        double specific_heat;
        double latent_heat;
        nlohmann::json viscosity_properties;
        Mesh* mesh;
    
    public:
        MagmaState(Mesh* mesh, nlohmann::json&& properties);
        void updateDensity(DikeData* dike) const;
        void updateViscosity(DikeData* dike) const;
        void updateEquilibriumCrystallization(DikeData* dike) const;
        double getThermalConductivity() const;
        double getSpecificHeat() const;
        inline double getLatentHeat() const{
            return latent_heat;
        }
};