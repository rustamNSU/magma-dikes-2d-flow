#pragma once
#include <Eigen/Dense>
#include "InputData.hpp"
#include "Mesh.hpp"
#include <nlohmann/json.hpp>


class MagmaState;
class MassBalance;
class DikeDataWriter;


class DikeData{
    private:
        InputData* input;
        Mesh* meshX;
        int ny = 1; // Number of layers into dike
        Eigen::VectorXd yc; // y/h(t, x)
        Eigen::VectorXd yb; // yb/h(t, x)

        Eigen::VectorXd hw; // halfwidth
        Eigen::VectorXd width; // width
        Eigen::VectorXd pressure;
        Eigen::VectorXd overpressure;

        Eigen::MatrixXd density; // Magma density in each layer of mesh element
        Eigen::MatrixXd temperature; // Magma temperature in each layer of mesh element
        Eigen::MatrixXd viscosity; // Magma viscosity in each layer of mesh element
        Eigen::MatrixXd qx; // u d\xi on x_{i+1/2} [nx+1, ny]
        Eigen::MatrixXd qy; // v dx on \xi_{j+1/2} [nx,   ny+1]
        Eigen::MatrixXd A; // A^j on x_{i+1/2}
        Eigen::MatrixXd C; // C^j on x_{i+1/2}
        Eigen::MatrixXd mobility; // Layer mobility into element
        Eigen::VectorXd Qx; // Total flux between elements
        double time;
        nlohmann::json algorithm_properties;
    
    public:
        DikeData(Mesh* mesh, const nlohmann::json& alg_properties);
        
        inline double getTime() const{
            return time;
        }

        inline int getLayersNumber() const{
            return ny;
        }

        // const Eigen::VectorXd& getDensity() const;
        // const Eigen::VectorXd& getWidth() const;
        // const Eigen::VectorXd& getPressure() const;
        // const Eigen::VectorXd& getOverpressure() const;
        // const Eigen::MatrixXd& getYc() const;
        // const Eigen::MatrixXd& getYb() const;
        // const Eigen::MatrixXd& getTemperature() const;
        // const Eigen::MatrixXd& getMobility() const;
        // const Eigen::VectorXd& getTotalMobility() const;

        // void setTime(double time);
        // void setDensity(const Eigen::VectorXd& vec);
        // void setWidth(const Eigen::VectorXd& vec);
        // void setPressure(const Eigen::VectorXd& vec);
        // void setTemperature(const Eigen::MatrixXd& mat);
        // void setViscosity(const Eigen::MatrixXd& mat);

        friend class MagmaState;
        friend class MassBalance;
        friend class DikeDataWriter;
};