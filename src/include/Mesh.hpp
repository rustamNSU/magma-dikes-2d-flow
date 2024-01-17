#pragma once
#include <Eigen/Dense>
#include <tuple>

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
        Mesh(int n, double xmin, double xmax);

        const Eigen::VectorXd& getx() const;
        const Eigen::VectorXd& getxl() const;
        const Eigen::VectorXd& getxr() const;
        int size() const;
        double getdx() const;

};