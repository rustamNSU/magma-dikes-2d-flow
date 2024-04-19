#pragma once
#include <Eigen/Dense>
#include "Mesh.hpp"
#include <nlohmann/json.hpp>


class MagmaState;
class MassBalance;
class DikeDataWriter;
enum class FlowModel{
    dike,
    channel
};


class DikeData{
    public:
        using json = nlohmann::json;

    private:
        FlowModel model = FlowModel::dike;
        Mesh* mesh;
        int ny = 1; // Number of layers into dike
        Eigen::VectorXd width;
        Eigen::VectorXd density;
        Eigen::VectorXd pressure;
        Eigen::VectorXd overpressure;

        Eigen::MatrixXd yc; // y-center of layers (half width)
        Eigen::MatrixXd yb; // y-boundary of layers (half width)
        Eigen::MatrixXd temperature; // Magma temperature in each layer of mesh element
        Eigen::MatrixXd viscosity; // Magma viscosity in each layer of mesh element
        Eigen::MatrixXd mobility; // Layer mobility into element
        Eigen::VectorXd total_mobility; // Total mobility into element
        Eigen::VectorXd total_face_mobility; // Mobility between elements
        Eigen::MatrixXd face_flux;
        Eigen::VectorXd total_face_flux;
        double time;
        json algorithm_properties;
    
    public:
        DikeData(Mesh* mesh, const json& alg_properties);
        double getTime() const;
        int getLayersNumber() const;
        const Eigen::VectorXd& getDensity() const;
        const Eigen::VectorXd& getWidth() const;
        const Eigen::VectorXd& getPressure() const;
        const Eigen::VectorXd& getOverpressure() const;
        const Eigen::MatrixXd& getYc() const;
        const Eigen::MatrixXd& getYb() const;
        const Eigen::MatrixXd& getTemperature() const;
        const Eigen::MatrixXd& getMobility() const;
        const Eigen::VectorXd& getTotalMobility() const;

        void setTime(double time);
        void setDensity(const Eigen::VectorXd& vec);
        void setWidth(const Eigen::VectorXd& vec);
        void setPressure(const Eigen::VectorXd& vec);
        void setTemperature(const Eigen::MatrixXd& mat);
        void setViscosity(const Eigen::MatrixXd& mat);

        friend class MagmaState;
        friend class MassBalance;
        friend class DikeDataWriter;
};