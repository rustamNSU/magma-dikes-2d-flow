#pragma once
#include <Eigen/Dense>
#include "Mesh.hpp"

class CohesiveModel{
    private:
        double E;
        double nu;
        double Ep;
        Mesh* mesh;

        /* Cohesive stress parameters */
        double dc = 0.0;
        double Gc = 0.0;
        double K1c = 0.0;
        double sigmac = 0.0;
        int Ncoh = 6;
        double lc = 0.0;
        CohesiveModel() = delete;

    public:
        CohesiveModel(
            Mesh* mesh,
            double E,
            double nu,
            double K1c,
            int Ncoh = 6
        );

        inline double getCohesiveStress(double h) const{
            double w = 2*h;
            if (w < dc) return sigmac * (1.0 - w / dc);
            return 0.0;
        }


        inline double getCohesiveDerivative(double h) const{
            double w = 2*h;
            if (w < dc) return -2.0 * sigmac / dc;
            return 0.0;
        }

        inline int getCohesiveElementsSize() const{
            return Ncoh;
        }
};