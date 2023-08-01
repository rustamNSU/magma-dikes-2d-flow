#pragma once
#include <Eigen/Dense>
#include "Mesh.hpp"


class DikeData{
    private:
        Mesh* mesh;
        Eigen::VectorXd width;
        Eigen::VectorXd rho;
        Eigen::VectorXd pressure;
    
    public:
        DikeData(Mesh* mesh);
};