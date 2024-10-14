#pragma once
#include <Eigen/Dense>
#include "Mesh.hpp"

class Elasticity{
    private:
        double E;
        double nu;
        double Ep;
        Mesh* mesh;
        Eigen::MatrixXd matrix; // C

        int step = 1;
        Eigen::MatrixXd A; // Sparse diagonal from C with (window = 2*step + 1)
        Eigen::MatrixXd B; // C - A

        Elasticity() = delete;

    public:
        Elasticity(
            double E,
            double nu,
            Mesh* mesh
        );

        void generateMatrix();
        const Eigen::MatrixXd& getMatrix() const;


        inline const Eigen::MatrixXd& getA() const{
            return A;
        }


        inline const Eigen::MatrixXd& getB() const{
            return B;
        }
};