#include "Schedule.hpp"
#include "Utils.hpp"

using Utils::mean_piecewise;


Schedule::Schedule(Mesh* mesh, nlohmann::json&& properties) : 
    mesh(mesh), 
    properties(properties)
{
    parseProperties();
}


void Schedule::parseProperties(){
    qlist = properties["Q"].get<std::vector<double>>();
    tlist = properties["t"].get<std::vector<double>>();
    rho = properties["rho"].get<double>();
    T = properties["T"].get<double>();
    beta = properties["beta"].get<double>();
}


double Schedule::getMassRate(double t1, double t2) const{
    return rho * mean_piecewise(tlist, qlist, t1, t2);
}