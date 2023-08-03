#pragma once
#include <Eigen/Dense>
#include "Mesh.hpp"


class MagmaState;
class MassBalance;
class DikeDataWriter;
class DikeData{
    private:
        Mesh* mesh;
        Eigen::VectorXd width;
        Eigen::VectorXd density;
        Eigen::VectorXd pressure;
        Eigen::VectorXd overpressure;
        Eigen::VectorXd viscosity;
        Eigen::VectorXd mobility;    // Magma mobility (w^3/12mu)
        double time;
        double MIN_MOBILITY_WIDTH = 1e-10;
    
    public:
        DikeData(Mesh* mesh);
        double getTime() const;
        const Eigen::VectorXd& getMobility() const;
        const Eigen::VectorXd& getDensity() const;
        const Eigen::VectorXd& getWidth() const;
        const Eigen::VectorXd& getPressure() const;

        void setTime(double time);
        void setDensity(const Eigen::VectorXd& vec);
        void setWidth(const Eigen::VectorXd& vec);
        void setViscosity(const Eigen::VectorXd& vec);
        void setPressure(const Eigen::VectorXd& vec);

        void calculateMobility();

        friend class MagmaState;
        friend class MassBalance;
        friend class DikeDataWriter;
};