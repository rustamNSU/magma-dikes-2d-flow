#include "Schedule.hpp"


Schedule::Schedule(
    Mesh* mesh,
    const std::vector<double>& qlist,
    const std::vector<double>& tlist
) : mesh(mesh), qlist(qlist), tlist(tlist)
{

}
