#pragma once
#include <memory>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include "Mesh.hpp"
#include "DikeData.hpp"
#include "InputData.hpp"
#include "Schedule.hpp"
#include "MagmaState.hpp"
#include "TimestepController.hpp"
#include "ReservoirData.hpp"

class ChannelFlow{
    private:
        std::string input_path;
        std::shared_ptr<InputData> input;
        std::shared_ptr<ReservoirData> reservoir;
        std::shared_ptr<TimestepController> timestep_controller;
        std::shared_ptr<Mesh> mesh;
        std::shared_ptr<Schedule> schedule;
        std::shared_ptr<MagmaState> magma_state;
        std::shared_ptr<DikeData> dike;
        std::shared_ptr<DikeData> old_dike;

        Eigen::MatrixXd mat;
        Eigen::VectorXd rhs;
        nlohmann::json algorithm_properties;

        std::string TIMESTEP_SCHEME;
        int MAX_ITERATIONS = 50;
        int MIN_STAB_ITERATIONS = 2;
        double TOLERANCE = 1e-4;
        double MIN_MOBILITY_WIDTH = 1e-10;
        double CUTOFF_VELOCITY = 1e-4;
        double CFL_FACTOR = 0.1;

    public:
        ChannelFlow(std::string input_path);
        void updateData();

};