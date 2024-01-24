#include "Schedule.hpp"
#include "Utils.hpp"

using Utils::mean_piecewise;


Schedule::Schedule(
    Mesh* mesh,
    const std::vector<double>& qlist,
    const std::vector<double>& tlist,
    double rho,
    double T
) : mesh(mesh), qlist(qlist), tlist(tlist), rho(rho), T(T)
{

}


double Schedule::getMassRate(double t1, double t2) const{
    return rho * mean_piecewise(tlist, qlist, t1, t2);
}