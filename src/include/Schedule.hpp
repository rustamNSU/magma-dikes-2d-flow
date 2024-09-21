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
        double beta = 0.0;
        void parseProperties();

    public:
        Schedule(Mesh* mesh, nlohmann::json&& properties);

        double getMassRate(double t1, double t2) const;
        inline double getMagmaChamberTemperature() const{
            return T;
        }
        
        inline double getMagmaChamberDensity() const{
            return rho;
        }

        inline double getMagmaChamberCrystallization() const{
            return beta;
        }

        inline void setMagmaChamberCrystallization(double beta_chamber){
            beta = beta_chamber;
        }
};