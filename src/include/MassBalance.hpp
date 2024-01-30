#pragma once
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include "Mesh.hpp"
#include "Elasticity.hpp"
#include "DikeData.hpp"
#include "InputData.hpp"
#include "Schedule.hpp"
#include "MagmaState.hpp"
#include "TimestepController.hpp"
#include "ReservoirData.hpp"

class MassBalance{
    private:
        InputData* input;
        ReservoirData* reservoir;
        TimestepController* timestep_controller;
        Elasticity* elasticity;
        Mesh* mesh;
        Schedule* schedule;
        MagmaState* magma_state;
        DikeData* new_dike;
        DikeData* old_dike;

        Eigen::MatrixXd mat;
        Eigen::VectorXd rhs;

        std::string TIMESTEP_SCHEME;
        int MAX_ITERATIONS = 50;
        int MIN_STAB_ITERATIONS = 2;
        double TOLERANCE = 1e-4;
        double MIN_MOBILITY_WIDTH = 1e-10;
        double CUTOFF_VELOCITY = 1e-4;
        double CFL_FACTOR = 0.1;

    public:
        struct explicitSolverOutput{
            bool cfl_condition = true;
            int ratio = 1;
        };
        MassBalance(
            InputData* input,
            TimestepController* timestep_controller,
            Elasticity* elasticity,
            Mesh* mesh,
            Schedule* schedule,
            MagmaState* magma_state,
            ReservoirData* reservoir
        );

        void setNewTimestepData(
            DikeData* new_dike,
            DikeData* old_dike
        );

        void setAlgorithmProperties(const nlohmann::json& properties);
        explicitSolverOutput explicitSolve();
        void updatePressure();
        void updateTemperature();
        // SolverOutput solve();
        void generateMatrix();
        void calculateMobility();
};