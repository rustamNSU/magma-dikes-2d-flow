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

        /* Cohesive stress parameters */
        double dc = 1e-3;
        double dm = 1e-5;
        double Gc = 0.0;
        double sigmac = 0.0;

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


        void setCohesiveParameters(double dc, double dm, double Gc);
        inline double getCohesiveStress(double h) const{
            double w = 2*h;
            if (w < dm) return sigmac * w / dm;
            if (w < dc) return sigmac * (dc - w) / (dc - dm);
            return 0.0;
        }
};