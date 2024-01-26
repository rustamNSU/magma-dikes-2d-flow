#pragma once
#include <vector>
#include <nlohmann/json.hpp>
#include "Mesh.hpp"

class Schedule{
    private:
        Mesh* mesh;
        nlohmann::json properties;
        std::vector<double> qlist;
        std::vector<double> tlist;
        double rho;
        double T;
        void parseProperties();

    public:
        Schedule(Mesh* mesh, nlohmann::json&& properties);

        double getMassRate(double t1, double t2) const;
        double getMagmaChamberTemperature() const;
        double getMagmaChamberDensity() const;
};