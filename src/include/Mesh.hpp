#pragma once
#include <Eigen/Dense>
#include <tuple>
#include <nlohmann/json.hpp>

class Mesh{
    private:
        int n;
        Eigen::VectorXd x;  ///< coordinate of centers of i-cell
        Eigen::VectorXd xl; ///< coordinate of left face of i-cell
        Eigen::VectorXd xr; ///< coordinate of right face of i-cell
        double dx;
        double xmin;
        double xmax;
        
    public:
        Mesh(const nlohmann::json& properties);
        
        const Eigen::VectorXd& getx() const;
        const Eigen::VectorXd& getxl() const;
        const Eigen::VectorXd& getxr() const;
        int size() const;
        double getdx() const;

};