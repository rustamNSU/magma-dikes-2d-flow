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

class DikeModel2d{
    private:
        std::string input_path;
        std::shared_ptr<InputData> input;
        std::shared_ptr<ReservoirData> reservoir;
        std::shared_ptr<Elasticity> elasticity;
        std::shared_ptr<TimestepController> timestep_controller;
        std::shared_ptr<Mesh> mesh;
        std::shared_ptr<Schedule> schedule;
        std::shared_ptr<MagmaState> magma_state;
        std::shared_ptr<DikeData> dike;
        std::shared_ptr<DikeData> old_dike;
        nlohmann::json algorithm_properties;

        std::string TIMESTEP_SCHEME;
        int MAX_ITERATIONS = 50;
        int MIN_STAB_ITERATIONS = 2;
        double TOLERANCE = 1e-4;
        double MIN_MOBILITY_WIDTH = 1e-10;
        double CUTOFF_VELOCITY = 1e-4;
        double CFL_FACTOR = 0.1;

    public:
        DikeModel2d(const std::string& input_path);
        void setInitialData();
        void run();
        void explicitSolver();
        void updatePressure();
        void calculateVerticalFlow();
        void solveMassBalance();
};