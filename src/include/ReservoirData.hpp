#pragma once
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include "Mesh.hpp"


class MagmaState;
class MassBalance;
class DikeDataWriter;
class DikeModel2d;

class ReservoirData{
    public:
        friend class MagmaState;
        friend class MassBalance;
        friend class DikeDataWriter;
        friend class DikeModel2d;

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
        Eigen::ArrayXd yc; // Center of layers
        Eigen::ArrayXd yb; // Coordinates of boundaries between layers
        Eigen::ArrayXd dy; // Width of layers
        Eigen::ArrayXd initial_temperature; // Initial reservoir temperature
        Eigen::ArrayXd density;
        Eigen::ArrayXd lithostatic_pressure;
        
        double time = 0;
        Eigen::ArrayXXd temperature; // Current reservoir temperature field
        Eigen::ArrayXXd temperature_old; // Current reservoir temperature field
    
    private:
        void generateCosineRefinementMesh();
        void calculateDensity();
        void calculateLithostaticPressure();
        void calculateReservoirTemperature();

    public:
        ReservoirData(Mesh* mesh, nlohmann::json properties);
        const Eigen::ArrayXd& getInitialTemperature() const;
        const Eigen::ArrayXXd& getTemperature() const;
        void setTemperature(Eigen::ArrayXXd&& T);
        std::tuple<double, double, double> getElasticityParameters() const;
        double getGravityAcceleration() const;
        const Eigen::ArrayXd& getLithostaticPressure() const;
};