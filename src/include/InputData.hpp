#pragma once
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <string>
#include <tuple>
#include <filesystem>
#include <spdlog/spdlog.h>
#include <memory>
#include "Mesh.hpp"


class InputData{
    private:
        int sim_id;
        nlohmann::json input;
        std::filesystem::path work_dir;
        std::filesystem::path sim_dir;
        std::filesystem::path data_dir;
        std::filesystem::path reservoir_dir;
    
    public:
        InputData(const nlohmann::json& input);
    
    private:
        void saveInputJson();
        void setMultiSinkLogger();

    public:
        nlohmann::json getTimestepProperties() const;
        nlohmann::json getMeshProperties() const;
        nlohmann::json getReservoirProperties() const;
        nlohmann::json getAlgorithmProperties() const;
        nlohmann::json getScheduleProperties() const;
        nlohmann::json getMagmaProperties() const;
        const std::filesystem::path& getSimDir() const;
        const std::filesystem::path& getDataDir() const;
        const std::filesystem::path& getReservoirDir() const;

        template<typename T>
        T get(const std::string& key) const{
            return input[key].get<T>();
        }

//         /* Access via json pointer: "/k1/k2/k3"_json_pointer */
//         template<typename T, typename Key>
//         T get(Key&& key) const{
//             return input[key].get<T>();
//         }
};