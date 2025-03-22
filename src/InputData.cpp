#include "InputData.hpp"
#include "fstream"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
using json = nlohmann::json;
namespace fs = std::filesystem;


InputData::InputData(const json& input) : input(input)
{
    sim_id = input["sim_id"];
    work_dir = fs::current_path();
    sim_dir = work_dir / ("simulations/simID" + std::to_string(sim_id));
    data_dir = sim_dir / "data";
    reservoir_dir = sim_dir / "reservoir_data";
    for (auto dir_path : {data_dir, reservoir_dir}){
        std::filesystem::remove_all(dir_path);
        if (std::filesystem::exists(dir_path) == false){
            std::filesystem::create_directories(dir_path);
        }
    }
    saveInputJson();
    setMultiSinkLogger();
}


void InputData::saveInputJson(){
    std::ofstream fout(sim_dir / "input.json", std::ios::trunc);
    fout << std::setw(4) << input << std::endl;
}


void InputData::setMultiSinkLogger()
{
    auto log_path = sim_dir / "log.txt";

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    if (input["algorithm_properties"]["is_debug"].get<bool>()) {
        console_sink->set_level(spdlog::level::debug);
    } else {
        console_sink->set_level(spdlog::level::info);
    }
    console_sink->set_pattern("%^[%=8l] %v %$");

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_path.string(), true);
    file_sink->set_level(spdlog::level::trace);
    file_sink->set_pattern("[%=8l] %v");

    auto logger = std::make_shared<spdlog::logger>("multi_sink");
    logger->sinks().push_back(console_sink);
    logger->sinks().push_back(file_sink);
    logger->set_level(spdlog::level::trace);
    logger->flush_on(spdlog::level::info);
    spdlog::set_default_logger(logger);
}


json InputData::getTimestepProperties() const{
    return input["timestep_properties"];
}


json InputData::getMeshProperties() const{
    return input["mesh_properties"];
}


json InputData::getReservoirProperties() const{
    return input["reservoir_properties"];
}


json InputData::getAlgorithmProperties() const{
    return input["algorithm_properties"];
}


json InputData::getScheduleProperties() const{
    return input["schedule_properties"];
}


json InputData::getMagmaProperties() const{
    return input["magma_properties"];
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