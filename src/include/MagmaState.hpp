#pragma once
#include <string>
#include <nlohmann/json.hpp>
#include "Mesh.hpp"
#include "DikeData.hpp"


class MagmaState{
    private:
        std::string density_model = "constant_density";
        std::string viscosity_model = "constant_viscosity";
        double rho;
        double mu;
        Mesh* mesh;
    
    public:
        MagmaState(const nlohmann::json& input, Mesh* mesh);
        void updateDensity(DikeData* dike) const;
        void updateViscosity(DikeData* dike) const;
};