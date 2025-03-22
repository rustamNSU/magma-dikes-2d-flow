#include "Schedule.hpp"
#include "Utils.hpp"

using Utils::mean_piecewise;


Schedule::Schedule(Mesh* mesh, nlohmann::json&& properties) : 
    mesh(mesh), 
    properties(properties)
{
    parseProperties();
}


void Schedule::parseProperties() {
    qlist = properties["volume_rate"].get<std::vector<double>>();
    tlist = properties["time"].get<std::vector<double>>();
    rho = properties["density"].get<double>();
    T = properties["temperature"].get<double>();
    beta = properties["beta"].get<double>();
}


double Schedule::getMassRate(double t1, double t2) const {
    return rho * mean_piecewise(tlist, qlist, t1, t2);
}