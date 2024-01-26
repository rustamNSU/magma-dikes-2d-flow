#include <iostream>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <iostream>

#include "InputData.hpp"
#include "Mesh.hpp"
#include "ReservoirData.hpp"
#include "Elasticity.hpp"
#include "DikeData.hpp"
#include "MassBalance.hpp"
#include "Writers.hpp"

using json = nlohmann::json;
using Eigen::VectorXd;
using Eigen::MatrixXd;
namespace fs = std::filesystem;



int main(int argc, char ** argv){
    std::string input_path = argv[1];
	std::ifstream f(input_path);
    json input_json = json::parse(f);
	f.close();

    InputData input(input_json);
    TimestepController timestep_controller(input.getTimestepProperties());
    Mesh mesh(input.getMeshProperties());
    ReservoirData reservoir(&mesh, input.getReservoirProperties());

    auto [E, nu, KIc] = reservoir.getElasticityParameters();
    Elasticity elasticity(E, nu, &mesh);

    /* @todo: pls, refactor me */
    auto algorithm_properties = input.getAlgorithmProperties();
    DikeData dike(&mesh, algorithm_properties["numberOfLayers"]);
    Schedule schedule(&mesh, input.getScheduleProperties());
    MagmaState magma_state(&mesh, input.getMagmaProperties());
    magma_state.updateDensity(&dike);
    magma_state.updateViscosity(&dike);
    MassBalance mass_balance(
        &input,
        &timestep_controller,
        &elasticity,
        &mesh,
        &schedule,
        &magma_state,
        &reservoir
    );
    mass_balance.setAlgorithmProperties(input.getAlgorithmProperties());
    DikeDataWriter writer;
    auto [is_save_timestep, save_timestep] = timestep_controller.saveTimestepIteration();
    std::string savepath = (input.getDataDir() / "data_").string() + std::to_string(save_timestep) + ".h5";
    writer.saveData(&dike, savepath);
    DikeData old_dike = dike;
    while (timestep_controller.isFinish()){
        double current_time = timestep_controller.getCurrentTime();
		double dt = timestep_controller.getCurrentTimestep();
		double time_iteration = timestep_controller.getTimeIteration();
        auto [is_save, save_timestep] = timestep_controller.saveTimestepIteration();
        std::cout << time_iteration + 1 << ";  " << save_timestep << ") " << int(current_time / 1000) << "e3 s -> " << int((current_time + dt)/1000) << "e3 s." << std::endl;
        dike.setTime(current_time + dt);
        mass_balance.setNewTimestepData(&dike, &old_dike);
        mass_balance.explicitSolve();
        timestep_controller.update();
        auto level = timestep_controller.getLevel();
        std::tie(is_save, save_timestep) = timestep_controller.saveTimestepIteration();
        if (level == 0 && is_save){
            savepath = (input.getDataDir() / "data_").string() + std::to_string(save_timestep) + ".h5";
            writer.saveData(&dike, savepath);
        }
        old_dike = dike;
    }
    return 0;
}

// int main(int argc, char ** argv){
//     std::string input_path = argv[1];
// 	std::ifstream f(input_path);
//     json input_json = json::parse(f);
// 	f.close();

//     int simID = input_json["simID"];
//     auto work_dir = fs::current_path();
//     auto sim_dir = work_dir / ("simulations/simID" + std::to_string(simID));
//     auto data_dir = sim_dir / "data";
//     fs::create_directories(data_dir);
//     fs::copy(input_path, sim_dir / "input.json", fs::copy_options::overwrite_existing);

//     Mesh mesh(
//         input_json["meshProperties"]["n"],
//         input_json["meshProperties"]["xmin"],
//         input_json["meshProperties"]["xmax"]
//     );
//     InputData input(input_json, &mesh);
//     auto [E, nu, KIc] = input.getElasticityParameters();
//     Elasticity elasticity(E, nu, &mesh);
//     DikeData dike(&mesh);
//     Schedule schedule(
//         &mesh,
//         input_json["scheduleProperties"]["Q"],
//         input_json["scheduleProperties"]["t"],
//         input_json["scheduleProperties"]["rho"]
//     );
//     MagmaState magma_state(input_json, &mesh);
//     magma_state.updateDensity(&dike);
//     magma_state.updateViscosity(&dike);
//     MassBalance mass_balance(
//         &input,
//         &elasticity,
//         &mesh,
//         &schedule,
//         &magma_state
//     );
//     mass_balance.setAlgorithmProperties(input_json["algorithmProperties"]);

//     double start_time = input_json["simulationProperties"]["startTime"];
//     double end_time = input_json["simulationProperties"]["endTime"];
//     double dt = input_json["simulationProperties"]["dt"];
//     double current_time = start_time;
//     int time_iteration = 0;
//     DikeDataWriter writer;
//     std::string savepath = (data_dir / "data_").string() + std::to_string(time_iteration) + ".h5";
//     writer.saveData(&dike, savepath);
//     while (current_time < end_time - 1e-6){
//         std::cout << time_iteration + 1 << ")    " << int(current_time) << " -> " << int(current_time + dt) << std::endl;
//         DikeData old_dike = dike;
//         dike.setTime(current_time + dt);
//         mass_balance.setNewTimestepData(&dike, &old_dike);
//         auto solver_output = mass_balance.solve();

//         std::cout << "  MassBalance:\n"
//                   << "    has converged: " << solver_output.is_converge << "\n"
//                   << "    iters: " << solver_output.iters << "\n"
//                   << "    error: " << solver_output.error << std::endl;
//         current_time += dt;
//         ++time_iteration;
//         savepath = (data_dir / "data_").string() + std::to_string(time_iteration) + ".h5";
//         writer.saveData(&dike, savepath);
//     }
//     return 0;
// }