#pragma once
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include "Mesh.hpp"


class MagmaState;
class MassBalance;
class DikeDataWriter;

class ReservoirData{
    public:
        friend class MagmaState;
        friend class MassBalance;
        friend class DikeDataWriter;

    private:
        Mesh* mesh;
        nlohmann::json properties;

        double E; // Young modulus
        double nu; // Poisson coefficient
        double KIc; // Fracture toughness
        double g = 9.81; // Gravity acceleration
        double C; // Reservoir specific heat
        double k; // Reservoir conductivity
        std::string density_model;
        std::string temperature_model;

        int ny; // Number of layers
        double L; // Reservoir width
        Eigen::VectorXd yc; // Center of layers
        Eigen::VectorXd yb; // Coordinates of boundaries between layers
        Eigen::VectorXd dy; // Width of layers
        Eigen::VectorXd initial_temperature; // Initial reservoir temperature
        Eigen::VectorXd density;
        Eigen::VectorXd lithostatic_pressure;
        
        double time = 0;
        Eigen::MatrixXd temperature; // Current reservoir temperature field
    
    private:
        void generateCosineRefinementMesh();
        void calculateDensity();
        void calculateLithostaticPressure();
        void calculateReservoirTemperature();

    public:
        ReservoirData(Mesh* mesh, nlohmann::json properties);
        const Eigen::VectorXd& getInitialTemperature() const;
        const Eigen::MatrixXd& getTemperature() const;
        void setTemperature(Eigen::MatrixXd&& T);
        std::tuple<double, double, double> getElasticityParameters() const;
        double getGravityAcceleration() const;
        const Eigen::VectorXd& getLithostaticPressure() const;
};