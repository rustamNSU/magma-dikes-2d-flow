#include "TimestepController.hpp"
#include <cmath>
#include <cassert>
#include <filesystem>
using json = nlohmann::json;
namespace fs = std::filesystem;


TimestepController::TimestepController(const json& timestep_properties) : timestep_properties(timestep_properties)
{
    start_time = timestep_properties["start_time"];
    current_time = start_time;
    end_time = timestep_properties["end_time"];
    dt_list = timestep_properties["dt_list"].get<std::vector<double>>();
    dt_time = timestep_properties["dt_time"].get<std::vector<double>>();
    saverate_list = timestep_properties["saverate_list"].get<std::vector<int>>();
    
    base_dt = dt_list[0];
    current_dt = base_dt;
    level = 0;
    base_timesteps = 0;
    all_timesteps = 0;
    attempts = 1;
    output_save_rate = saverate_list[0];
}


void TimestepController::setBaseTimestep(double base_dt) {
    assert(level == 0 && "setBaseTimestep used in wrong time, level not equal zero!\n");
    this->base_dt = base_dt;
    current_dt = base_dt;
}


double TimestepController::getCurrentTimestep() const {
    return current_dt;
}


double TimestepController::getCurrentTime() const {
    return current_time;
}


double TimestepController::getNextTime() const {
    return current_time + current_dt;
}


int TimestepController::getTimeIteration() const {
    return base_timesteps;
}


std::tuple<bool, int> TimestepController::saveTimestepIteration() const {
    bool is_save = ((base_timesteps % output_save_rate) == 0) && (level == 0);
    int save_iteration = base_timesteps / output_save_rate;
    return std::make_tuple(is_save, save_iteration);
}


int TimestepController::getTimestepAttempts() const {
    return attempts;
}


int TimestepController::getLevel() const {
    return level;
}


void TimestepController::update() {
    current_time += current_dt;
    attempts += 1;
    nonlinear_iteration = 0;
    if (!dt_stack.empty()){
        dt_stack.pop();
        level_stack.pop();
    }
    if (dt_stack.empty()){
        level = 0;
        updateBaseTimestep();
        current_dt = base_dt;
        all_timesteps++;
        base_timesteps++;
        attempts = 1;
    } else {
        current_dt = dt_stack.top();
        level = level_stack.top();
        all_timesteps++;
    }
}


void TimestepController::divideTimestep(int steps) {
    if (!dt_stack.empty()){
        dt_stack.pop();
        level_stack.pop();
    }
    steps = std::max(2, steps);
    current_dt /= steps;
    level += 1;
    for (int i = 0; i < steps; i++){
        dt_stack.push(current_dt);
        level_stack.push(level);
    }
    attempts += 1;
    nonlinear_iteration = 0;
}


void TimestepController::updateNonlinearIteration() {
    nonlinear_iteration++;
}


bool TimestepController::isFinish() const {
    return !(current_time < end_time - 1e-10);
}


double TimestepController::getRelaxationParameter() const {
    int i = 0;
    for (auto rel_iter : relaxation_iterations) {
        if (nonlinear_iteration < rel_iter){
            return relaxation_parameters[i];
        }
        i++;
    }
    return 1.0;
}


int TimestepController::getNonlinearIteration() const {
    return nonlinear_iteration;
}


void TimestepController::updateBaseTimestep() {
    assert(level == 0 && "setBaseTimestep used in wrong time, level not equal zero!\n");
    if (base_dt_indx >= dt_list.size() - 1) {
        return;
    } else {
        if (std::abs(dt_time[base_dt_indx + 1] - current_time) < 0.1 * base_dt) {
            base_dt_indx++;
            base_dt = dt_list[base_dt_indx];
            output_save_rate = saverate_list[base_dt_indx];
            current_dt = base_dt;
        }
    }
}
