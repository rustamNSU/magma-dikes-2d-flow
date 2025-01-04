#pragma once
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <unordered_map>

#include "Mesh.hpp"
#include "Elasticity.hpp"
#include "CohesiveModel.hpp"
#include "DikeData.hpp"
#include "InputData.hpp"
#include "Schedule.hpp"
#include "MagmaState.hpp"
#include "TimestepController.hpp"
#include "ReservoirData.hpp"
#include "JustTimer.hpp"

class DikeModel2d{
    private:
        struct ExplicitSolverLog{
            bool successful = true;
            bool cfl_ratio = 2;

            inline void setDefault(){
                successful = true;
            }
        };


        struct DikeFront{
            std::vector<double> time;
            std::vector<double> front;
            
            inline void addTimestep(double t, double x){
                time.push_back(t);
                front.push_back(x);
            }
        };


        struct DikeFrontUnique{
            std::vector<double> time;
            std::vector<double> front;
            int last_tip = -1;
            
            inline void addTimestep(double t, double x){
                time.push_back(t);
                front.push_back(x);
            }
        };
        

    private:
        std::string input_path;
        std::shared_ptr<InputData> input;
        std::shared_ptr<ReservoirData> reservoir;
        std::shared_ptr<Elasticity> elasticity;
        std::shared_ptr<CohesiveModel> cohesive_model;
        std::shared_ptr<TimestepController> timestep_controller;
        std::shared_ptr<Mesh> mesh;
        std::shared_ptr<Schedule> schedule;
        std::shared_ptr<MagmaState> magma_state;
        std::shared_ptr<DikeData> dike;
        std::shared_ptr<DikeData> old_dike;
        nlohmann::json algorithm_properties;

        ExplicitSolverLog solver_log;
        std::string VISCOSITY_APPROXIMATION = "harmonic";
        int MAX_ITERATIONS = 50;
        int MIN_STAB_ITERATIONS = 2;
        double TOLERANCE = 1e-6;
        double MIN_MOBILITY_WIDTH = 1e-4;
        double MIN_WIDTH = 1e-3;
        double SHEAR_HEATING = 1.0;
        double LATENT_HEAT = 1.0;
        bool highOrderApproximation = false;

        int DENSITY_MODEL = 0;
        int VISCOSITY_MODEL = 0;
        double total_injected_energy = 0.0;
        double total_injected_mass = 0.0;
        DikeFront dike_front;
        DikeFrontUnique dike_front_unique;

    public:
        DikeModel2d(const std::string& input_path);
        void setInitialData();
        void setAlgorithmProperties();
        void run();
        void implicitSolver();
        // void calculateVerticalFlow();
        // void solveMassBalance();
        void implicitMassBalance();
        void solveEnergyBalance();
        void updateCrystallization();
        void reloadData();
        void updateData();
        void saveData(const std::string &savepath);
        void printTimer() const;
    
    public:
        struct SummaryTable{
            std::unordered_map<std::string, double> error;
            std::unordered_map<std::string, double> table;

            inline void clear(){
                error.clear();
                table.clear();
            }

            void print() const;
        };
    
    private: 
        SummaryTable summary;
};