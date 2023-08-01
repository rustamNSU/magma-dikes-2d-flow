#pragma once
#include <Eigen/Dense>
#include "Mesh.hpp"

class Elasticity{
    private:
        double E;
        double nu;
        double Ep;
        Mesh* mesh;
        Eigen::MatrixXd matrix;

        Elasticity() = delete;

    public:
        Elasticity(
            double E,
            double nu,
            Mesh* mesh
        );

        void generateMatrix();
        const Eigen::MatrixXd& get_matrix() const;
};