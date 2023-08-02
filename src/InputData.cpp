#include "InputData.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using json = nlohmann::json;


InputData::InputData(const json input, Mesh* mesh) : mesh(mesh){
    int n = mesh->size();
    E = input["reservoirProperties"]["E"];
    nu = input["reservoirProperties"]["nu"];
    g = input["reservoirProperties"]["g"];
    KIc = input["reservoirProperties"]["KIc"];
    plith = VectorXd::Zero(n);
    rhoR = VectorXd::Zero(n);
}


const VectorXd& InputData::getPlith() const{
    return plith;
}


double InputData::getg() const{
    return g;
}