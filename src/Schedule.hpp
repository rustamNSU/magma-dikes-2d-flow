#pragma once
#include <vector>
#include "Mesh.hpp"


class Schedule{
    private:
        Mesh* mesh;
        std::vector<double> qlist;
        std::vector<double> tlist;

    public:
        Schedule(
            Mesh* mesh,
            const std::vector<double>& qlist,
            const std::vector<double>& tlist
        );
};