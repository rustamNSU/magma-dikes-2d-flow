#pragma once
#include <Eigen/Dense>
#include <tuple>

class Mesh{
    private:
        int n;
        Eigen::VectorXd x;
        Eigen::VectorXd xl;
        Eigen::VectorXd xr;
        double dx;
        double xmin;
        double xmax;
        
    public:
        Mesh(int n, double xmin, double xmax);

        const Eigen::VectorXd& getx() const;
        const Eigen::VectorXd& getxl() const;
        const Eigen::VectorXd& getxr() const;
        int size() const;

};