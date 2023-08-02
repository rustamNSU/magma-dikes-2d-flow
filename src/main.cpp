#include "Mesh.hpp"
#include "Elasticity.hpp"
#include "DikeData.hpp"
#include "InputData.hpp"
#include "MassBalance.hpp"
#include "Writers.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>
#include <filesystem>

using json = nlohmann::json;
using Eigen::VectorXd;
using Eigen::MatrixXd;
namespace fs = std::filesystem;


int main(int argc, char ** argv){
    std::string input_path = argv[1];
	std::ifstream f(input_path);
    json input_json = json::parse(f);
	f.close();

    int simID = input_json["simID"];
    auto work_dir = fs::current_path();
    auto sim_dir = work_dir / ("simulations/simID" + std::to_string(simID));
    auto data_dir = sim_dir / "data";
    fs::create_directories(data_dir);
    fs::copy(input_path, sim_dir / "input.json", fs::copy_options::overwrite_existing);

    Mesh mesh(
        input_json["meshProperties"]["n"],
        input_json["meshProperties"]["xmin"],
        input_json["meshProperties"]["xmax"]
    );
    Elasticity elasticity(1.0, 0.2, &mesh);
    DikeData dike(&mesh);
    InputData input(input_json, &mesh);
    Schedule schedule(
        &mesh,
        input_json["scheduleProperties"]["Q"],
        input_json["scheduleProperties"]["t"],
        input_json["scheduleProperties"]["rho"]
    );
    MagmaState magma_state(input_json, &mesh);
    magma_state.updateDensity(&dike);
    magma_state.updateViscosity(&dike);
    MassBalance mass_balance(
        &input,
        &elasticity,
        &mesh,
        &schedule,
        &magma_state
    );

    double start_time = input_json["simulationProperties"]["startTime"];
    double end_time = input_json["simulationProperties"]["endTime"];
    double dt = input_json["simulationProperties"]["dt"];
    double current_time = start_time;
    int time_iteration = 0;
    DikeDataWriter writer;
    std::string savepath = (data_dir / "data_").string() + std::to_string(time_iteration) + ".h5";
    writer.saveData(&dike, savepath);
    while (current_time < end_time - 1e-6){
        DikeData old_dike = dike;
        dike.setTime(current_time + dt);
        mass_balance.setNewTimestepData(&dike, &old_dike);
        mass_balance.solve();

        current_time += dt;
        ++time_iteration;
        savepath = (data_dir / "data_").string() + std::to_string(time_iteration) + ".h5";
        writer.saveData(&dike, savepath);
    }
    return 0;
}