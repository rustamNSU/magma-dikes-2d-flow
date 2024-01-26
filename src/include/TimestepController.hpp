#pragma once
#include <stack>
#include <nlohmann/json.hpp>
#include <tuple>
#include <initializer_list>
#include <vector>

class TimestepController{
    private:
        nlohmann::json timestep_properties;
        std::vector<double> dt_list;
        std::vector<double> dt_time;
        int base_dt_indx = 0;
        double start_time;
        double end_time;
        double current_time;
        double base_dt;
        double current_dt;
        int base_timesteps; // For numbering filenames of saved data
        int all_timesteps;
        int attempts; // Amount of timesteps into base timestep (dt = base_dt)
        int level;
        int output_save_rate;
        int nonlinear_iteration = 0;

		std::vector<int> relaxation_iterations = {4, 8, 15, 25, 100, 1000};
		std::vector<double> relaxation_parameters = {1.0, 0.7, 0.5, 0.3, 0.1, 0.01};
        std::stack<double> dt_stack;
        std::stack<int> level_stack;
        TimestepController() = delete;

    public:
        TimestepController(const nlohmann::json& timestep_properties);
        void setBaseTimestep(double base_dt);
        double getCurrentTimestep() const;
        double getCurrentTime() const;
        double getNextTime() const;
        int getTimeIteration() const;
        int getTimestepAttempts() const;
        int getLevel() const;
        void update();
        void updateBaseTimestep();
        void divideTimestep();
        bool isFinish() const;
        std::tuple<bool, int> saveTimestepIteration() const;
        double getRelaxationParameter() const;
        void updateNonlinearIteration();
        int getNonlinearIteration() const;
};