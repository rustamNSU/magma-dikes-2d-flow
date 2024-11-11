#pragma once
#include <Eigen/Dense>
#include "InputData.hpp"
#include "Mesh.hpp"
#include <nlohmann/json.hpp>
#include <vector>


class MagmaState;
class MassBalance;
class DikeDataWriter;
class DikeModel2d;


class DikeData{
    private:
        InputData* input;
        Mesh* meshX;
        int ny = 1; // Number of layers into dike
        Eigen::ArrayXd yc; // y/h(t, x)
        Eigen::ArrayXd yb; // yb/h(t, x)

        Eigen::ArrayXd hw; // halfwidth
        Eigen::ArrayXd pressure;
        Eigen::ArrayXd overpressure;

        Eigen::ArrayXXd density; // Total magma density in each layer of mesh element (rhom + rhoc + rhog)
        Eigen::ArrayXXd rhom; // Melt density in each layer of mesh element (rhol + rhod)
        Eigen::ArrayXXd rhoc; // Crystal density in each layer of mesh element
        Eigen::ArrayXXd rhog; // Exsolved gas density in each layer of mesh element
        Eigen::ArrayXXd rhom_liquid; // Melt phase density
        Eigen::ArrayXXd temperature; // Magma temperature in each layer of mesh element
        Eigen::ArrayXXd viscosity; // Magma viscosity in each layer of mesh element
        Eigen::ArrayXXd Tliquidus; // Liquidus temperature
        Eigen::ArrayXXd Tsolidus; // Solidus temperature
        Eigen::ArrayXXd betaeq; // Equlibrium crystallization
        Eigen::ArrayXXd beta; // Crystal volume concentration in melt
        Eigen::ArrayXXd alpha; // Exsolved gas volume concentration
        Eigen::ArrayXXd gamma; // Dissolved gas weight concentration in melt
        Eigen::ArrayXXd xh2og; // Gas water weight ratio
        Eigen::ArrayXXd xh2od; // Dissolved water weight ratio
        Eigen::ArrayXXd wth2o; // Dissolved h2o weight concentration in melt
        Eigen::ArrayXXd wtco2; // Dissolved co2 weight concentration in melt
        Eigen::ArrayXXd qx; // u d\xi on x_{i+1/2} [nx+1, ny]
        Eigen::ArrayXXd qy; // v dx on \xi_{j+1/2} [nx,   ny+1]
        Eigen::ArrayXXd mx; // Mass rates for elements (Ox)
        Eigen::ArrayXXd my; // Mass rates for elements (Oy)
        Eigen::ArrayXXd A; // A^j on x_{i+1/2}, without (dp/dx + rho*g)
        Eigen::ArrayXXd C; // C^j on x_{i+1/2}, without (dp/dx + rho*g)
        Eigen::ArrayXXd shear_heat;
        Eigen::ArrayXXd mobility; // Layer mobility into element
        Eigen::ArrayXd Qx; // Total flux between elements
        Eigen::ArrayXd Mx; // Total mass rate between elements
        Eigen::ArrayXd Twall;
        Eigen::ArrayXd G; // Total pressure gradient with buoyancy between elements
        Eigen::ArrayXd magma_to_rock_heat_flux;
        double time = 0.0;
        nlohmann::json algorithm_properties;
        std::vector<bool> open_elements; // Mesh elements which filled with magma
        int tip_element;
        double front;
        double MIN_WIDTH = 1e-10;
    
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

        void updateOpenElements();
        void setMagmaStateAfterTip();
        void setMagmaStateAfterTip(int ntip, int nout);
        double getTotalMass() const;
        Eigen::ArrayXd getElementsHalfMass() const;
        Eigen::ArrayXd getElementsHalfEnergy(double Cm) const;
        double calculateMassBalanceError(const Eigen::ArrayXd &Mold, double dt) const;

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