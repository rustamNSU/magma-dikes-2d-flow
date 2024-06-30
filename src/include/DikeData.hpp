#pragma once
#include <Eigen/Dense>
#include "InputData.hpp"
#include "Mesh.hpp"
#include <nlohmann/json.hpp>


class MagmaState;
class MassBalance;
class DikeDataWriter;
class DikeModel2d;


class DikeData{
    private:
        InputData* input;
        Mesh* meshX;
        int ny = 1; // Number of layers into dike
        Eigen::VectorXd yc; // y/h(t, x)
        Eigen::VectorXd yb; // yb/h(t, x)

        Eigen::VectorXd hw; // halfwidth
        Eigen::VectorXd pressure;
        Eigen::VectorXd overpressure;

        Eigen::MatrixXd density; // Magma density in each layer of mesh element
        Eigen::MatrixXd temperature; // Magma temperature in each layer of mesh element
        Eigen::MatrixXd viscosity; // Magma viscosity in each layer of mesh element
        Eigen::MatrixXd qx; // u d\xi on x_{i+1/2} [nx+1, ny]
        Eigen::MatrixXd qy; // v dx on \xi_{j+1/2} [nx,   ny+1]
        Eigen::MatrixXd A; // A^j on x_{i+1/2}, without (dp/dx + rho*g)
        Eigen::MatrixXd C; // C^j on x_{i+1/2}, without (dp/dx + rho*g)
        Eigen::MatrixXd mobility; // Layer mobility into element
        Eigen::VectorXd Qx; // Total flux between elements
        Eigen::ArrayXd Twall;
        double time = 0.0;
        nlohmann::json algorithm_properties;
    
    public:
        DikeData(Mesh* mesh, const nlohmann::json& alg_properties);
        
        inline double getTime() const{
            return time;
        }

        inline int getLayersNumber() const{
            return ny;
        }

        inline Eigen::ArrayXd getElementsVolume() const{
            return hw.array() * meshX->getdx();
        }

        inline void setInitialPressure(const Eigen::VectorXd& plith){
            assert(plith.size() == pressure.size());
            pressure = plith;
        }

        inline void setInitialTemperature(const Eigen::VectorXd& T0){
            assert(T0.size() == temperature.rows());
            for (int icol = 0; icol < ny; icol++){
                temperature.col(icol) = T0;
            }
            Twall = T0;
        }

        inline Eigen::ArrayXd getWidth() const{
            return 2.0 * hw;
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
        void save(const std::string path) const;
        friend class MagmaState;
        friend class MassBalance;
        friend class DikeDataWriter;
        friend class DikeModel2d;
};