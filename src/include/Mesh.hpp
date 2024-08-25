#pragma once
#include <Eigen/Dense>
#include <tuple>
#include <nlohmann/json.hpp>

class Mesh{
    private:
        int n;
        Eigen::ArrayXd x;  ///< coordinate of centers of i-cell
        Eigen::ArrayXd xl; ///< coordinate of left face of i-cell
        Eigen::ArrayXd xr; ///< coordinate of right face of i-cell
        double dx;
        double xmin;
        double xmax;
        
    public:
        Mesh(const nlohmann::json& properties);
        const Eigen::ArrayXd& getx() const;
        const Eigen::ArrayXd& getxl() const;
        const Eigen::ArrayXd& getxr() const;
        int size() const;
        double getdx() const;

};