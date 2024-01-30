#include "TimestepController.hpp"
#include <cmath>
using json = nlohmann::json;

TimestepController::TimestepController(const json& timestep_properties) :
    timestep_properties(timestep_properties)
{
    start_time = timestep_properties["startTime"];
    current_time = start_time;
    end_time = timestep_properties["endTime"];
    dt_list = timestep_properties["dtList"].get<std::vector<double>>();
    dt_time = timestep_properties["dtTime"].get<std::vector<double>>();
    base_dt = dt_list[0];
    current_dt = base_dt;
    level = 0;
    base_timesteps = 0;
    all_timesteps = 0;
    attempts = 1;
    output_save_rate = timestep_properties["outputSaveRate"];
}


void TimestepController::setBaseTimestep(double base_dt){
    assert(level == 0 && "setBaseTimestep used in wrong time, level not equal zero!\n");
    this->base_dt = base_dt;
    current_dt = base_dt;
}


double TimestepController::getCurrentTimestep() const{
    return current_dt;
}


double TimestepController::getCurrentTime() const{
    return current_time;
}


double TimestepController::getNextTime() const{
    return current_time + current_dt;
}


int TimestepController::getTimeIteration() const{
    return base_timesteps;
}


std::tuple<bool, int> TimestepController::saveTimestepIteration() const{
    bool is_save = (base_timesteps % output_save_rate) == 0;
    int save_iteration = base_timesteps / output_save_rate;
    return std::make_tuple(is_save, save_iteration);
}


int TimestepController::getTimestepAttempts() const{
    return attempts;
}


int TimestepController::getLevel() const{
    return level;
}


void TimestepController::update(){
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
        base_timesteps++;
        attempts = 1;
        return;
    } else{
        current_dt = dt_stack.top();
        level = level_stack.top();
        all_timesteps++;
        return;
    }
}


void TimestepController::divideTimestep(){
    if (!dt_stack.empty()){
        dt_stack.pop();
        level_stack.pop();
    }
    current_dt /= 2;
    level += 1;
    dt_stack.push(current_dt);
    dt_stack.push(current_dt);
    level_stack.push(level);
    level_stack.push(level);
    attempts += 1;
    nonlinear_iteration = 0;
    return;
}


void TimestepController::updateNonlinearIteration(){
    nonlinear_iteration++;
    return;
}


bool TimestepController::isFinish() const{
    return current_time < end_time - 1e-10;
}


double TimestepController::getRelaxationParameter() const{
    int i = 0;
    for (auto rel_iter : relaxation_iterations){
        if (nonlinear_iteration < rel_iter){
            return relaxation_parameters[i];
        }
        i++;
    }
}


int TimestepController::getNonlinearIteration() const{
    return nonlinear_iteration;
}


void TimestepController::updateBaseTimestep(){
    assert(level == 0 && "setBaseTimestep used in wrong time, level not equal zero!\n");
    if (base_dt_indx >= dt_list.size() - 1){
        return;
    }
    else{
        if (std::abs(dt_time[base_dt_indx+1] - current_time) < 0.1 * base_dt){
            base_dt_indx++;
            base_dt = dt_list[base_dt_indx];
            current_dt = base_dt;
        }
        return;
    }
}