#pragma once
#include <Eigen/Dense>
#include "Mesh.hpp"
#include "Elasticity.hpp"
#include "DikeData.hpp"
#include "InputData.hpp"
#include "Schedule.hpp"
#include "MagmaState.hpp"

class MassBalance{
    private:
        InputData* input;
        Elasticity* elasticity;
        Mesh* mesh;
        Schedule* schedule;
        MagmaState* magma_state;
        DikeData* new_dike;
        DikeData* old_dike;

        double old_time;
        double new_time;
        double dt;

        Eigen::MatrixXd mat;
        Eigen::VectorXd rhs;

        int MAX_ITERATIONS = 50;
        int MIN_STAB_ITERATIONS = 2;
        double TOLERANCE = 1e-4;
    
    public:
        MassBalance(
            InputData* input,
            Elasticity* elasticity,
            Mesh* mesh,
            Schedule* schedule,
            MagmaState* magma_state
        );

        void setNewTimestepData(
            DikeData* new_dike,
            DikeData* old_dike
        );

        void solve();

        void generateMatrix();
};