#pragma once
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include "Mesh.hpp"


class InputData{
    private:
        /* Surrounding rock properties */
        double E;
        double nu;
        double KIc;
        double g = 9.81;
        
        Eigen::VectorXd rhoR;
        Eigen::VectorXd plith;
    
    public:
        InputData(const nlohmann::json input);
        const Eigen::VectorXd& getPlith() const;
        double getg() const;
};