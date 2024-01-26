#include "InputData.hpp"
#include "fstream"
using json = nlohmann::json;
namespace fs = std::filesystem;


InputData::InputData(const json& input) : input(input)
{
    sim_id = input["simID"];
    work_dir = fs::current_path();
    sim_dir = work_dir / ("simulations/simID" + std::to_string(sim_id));
    data_dir = sim_dir / "data";
    reservoir_dir = sim_dir / "reservoir_data";
    fs::create_directories(data_dir);
    fs::create_directories(reservoir_dir);
    saveInputJson();
}


void InputData::saveInputJson(){
    std::ofstream fout(sim_dir / "input.json", std::ios::trunc);
    fout << std::setw(4) << input << std::endl;
    fout.close();
}


json InputData::getTimestepProperties() const{
    return input["timestepProperties"];
}


json InputData::getMeshProperties() const{
    return input["meshProperties"];
}


json InputData::getReservoirProperties() const{
    return input["reservoirProperties"];
}


json InputData::getAlgorithmProperties() const{
    return input["algorithmProperties"];
}


json InputData::getScheduleProperties() const{
    return input["scheduleProperties"];
}


json InputData::getMagmaProperties() const{
    return input["magmaProperties"];
}


const fs::path& InputData::getSimDir() const{
    return sim_dir;
}


const fs::path& InputData::getDataDir() const{
    return data_dir;
}


const fs::path& InputData::getReservoirDir() const{
    return reservoir_dir;
}