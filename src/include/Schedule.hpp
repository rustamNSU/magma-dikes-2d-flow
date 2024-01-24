#pragma once
#include <vector>
#include "Mesh.hpp"

class MassBalance;

class Schedule{
    private:
        Mesh* mesh;
        std::vector<double> qlist;
        std::vector<double> tlist;
        double rho;
        double T;

    public:
        Schedule(
            Mesh* mesh,
            const std::vector<double>& qlist,
            const std::vector<double>& tlist,
            double rho,
            double T = 0
        );

        double getMassRate(double t1, double t2) const;
        friend class MassBalance;
};