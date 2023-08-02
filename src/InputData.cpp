#include "InputData.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using json = nlohmann::json;


InputData::InputData(const json input){
    E = input["reservoirProperties"]["E"];
    nu = input["reservoirProperties"]["nu"];
    g = input["reservoirProperties"]["g"];
    KIc = input["reservoirProperties"]["KIc"];
}


const VectorXd& InputData::getPlith() const{
    return plith;
}


double InputData::getg() const{
    return g;
}