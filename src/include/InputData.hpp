#pragma once
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <string>
#include <tuple>
#include "Mesh.hpp"


class InputData{
    private:
        /* Surrounding rock properties */
        double E;
        double nu;
        double KIc;
        double g = 9.81;
        std::string density_model;
        
        Mesh* mesh;
        Eigen::VectorXd rhoR;
        Eigen::VectorXd plith;
    
    public:
        InputData(const nlohmann::json& input, Mesh* mesh);
        const Eigen::VectorXd& getPlith() const;
        double getg() const;
        std::tuple<double, double, double> getElasticityParameters() const;
        void calculateLithostaticPressure();
};