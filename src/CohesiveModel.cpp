#include "CohesiveModel.hpp"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

CohesiveModel::CohesiveModel(
    Mesh* mesh,
    double E,
    double nu,
    double K1c,
    int Ncoh
) : E(E), nu(nu), mesh(mesh), K1c(K1c), Ncoh(Ncoh)
{
    double Ep = E / (1.0 - nu*nu);
    double dx = mesh->getdx();
    lc = Ncoh * dx;
    Gc = K1c*K1c / Ep;
    sigmac = std::sqrt((9.0 * M_PI / 32.0 * Ep * Gc / lc));
    dc = 2 * Gc / std::max(1e-10, sigmac);
}
