#pragma once
#include <string>
#include <nlohmann/json.hpp>
#include "Mesh.hpp"
#include "DikeData.hpp"
#include "ViscosityModels.hpp"


class MagmaState{
    private:
        nlohmann::json magma_properties;
        std::string density_model = "constant_density";
        std::string viscosity_model = "constant_viscosity";
        double rho;
        double thermal_conductivity;
        double specific_heat;
        nlohmann::json viscosity_properties;
        Mesh* mesh;
    
    public:
        MagmaState(const nlohmann::json& magma_properties, Mesh* mesh);
        void updateDensity(DikeData* dike) const;
        void updateViscosity(DikeData* dike) const;
        double getThermalConductivity() const;
        double getSpecificHeat() const;
};