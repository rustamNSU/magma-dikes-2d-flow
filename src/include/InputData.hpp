#pragma once
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <string>
#include <tuple>
#include "Mesh.hpp"


class InputData{
    private:
        nlohmann::json input;
        /* Surrounding rock properties */
        double E;
        double nu;
        double KIc;
        double g = 9.81;
        std::string density_model;
        std::string reservoir_temperature_model;
        
        Mesh* mesh;
        Eigen::VectorXd reservoir_density;
        Eigen::VectorXd lithostatic_pressure;
        Eigen::VectorXd reservoir_temperature;
    
    public:
        InputData(const nlohmann::json& input, Mesh* mesh);
        const Eigen::VectorXd& getLithostaticPressure() const;
        const Eigen::VectorXd& getReservoirTemperature() const;
        double getGravityAcceleration() const;

        /* return <E, nu, KIc> */
        std::tuple<double, double, double> getElasticityParameters() const;
        void calculateLithostaticPressure();
        void calculateReservoirTemperature();
};