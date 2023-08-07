#include "MassBalance.hpp"
#include <cmath>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using json = nlohmann::json;

MassBalance::MassBalance(
    InputData* input,
    Elasticity* elasticity,
    Mesh* mesh,
    Schedule* schedule,
    MagmaState* magma_state
) : input(input), elasticity(elasticity), mesh(mesh), schedule(schedule), magma_state(magma_state)
{

}


void MassBalance::setNewTimestepData(
    DikeData* new_dike,
    DikeData* old_dike
){
    this->new_dike = new_dike;
    this->old_dike = old_dike;

    old_time = old_dike->getTime();
    new_time = new_dike->getTime();
    dt = new_time - old_time;
    return;
}


void MassBalance::solve(){
    int iter = 0;
    int min_stab_iter = 0;
    double err_width;
    bool is_convergence = false;
    std::vector<double> err_list;
    err_list.reserve(MAX_ITERATIONS);
    for (; iter < MAX_ITERATIONS; ++iter){
        auto witer = new_dike->getWidth();
        generateMatrix();
        VectorXd sol = mat.fullPivLu().solve(rhs);
        err_width = (witer - sol).norm() / std::max(1.0, witer.norm());
        err_list.push_back(err_width);
        new_dike->setWidth(sol);
        magma_state->updateDensity(new_dike);
        magma_state->updateViscosity(new_dike);
        new_dike->overpressure = elasticity->getMatrix() * new_dike->width;
        new_dike->pressure = new_dike->overpressure + input->getPlith();
        if (err_width < TOLERANCE){
            ++min_stab_iter;
        }
        else{
            min_stab_iter = 0;
        }
        if (min_stab_iter >= MIN_STAB_ITERATIONS){
            is_convergence = true;
            break;
        }
    }
    int a = 1;
    if (!is_convergence){
        a = 2;
    }
    return;
}


void MassBalance::generateMatrix(){
    int n = mesh->size();
    double dx = mesh->getdx();
    double Mch = schedule->getMassRate(old_time, new_time);
    calculateMobility();
    mat.resize(n, n);
    rhs.resize(n);
    mat.fill(0.0);
    rhs.fill(0.0);
    auto C = elasticity->getMatrix();
    auto plith = input->getPlith();

    rhs[0] += Mch / dx;
    for (int i = 1; i < n; i++){
        double lambda = new_dike->getMobility()[i];
        double rho = (new_dike->getDensity()[i] + new_dike->getDensity()[i-1]) / 2.0;
        auto c1 = -rho * lambda / dx / dx;
        auto qrow = c1 * (C.row(i).array() - C.row(i-1).array());
        mat.row(i-1) += qrow.matrix();
        mat.row(i) -= qrow.matrix();

        double qlith = -c1 * (plith[i] - plith[i-1]);
        rhs[i-1] += qlith;
        rhs[i] -= qlith;

        double c2 = rho * lambda / dx * rho * input->getg();
        rhs[i-1] += c2;
        rhs[i] -= c2;
    }

    for (int i = 0; i < n; i++){
        mat(i, i) += new_dike->getDensity()[i] / dt;
        rhs[i] += old_dike->getDensity()[i] * old_dike->getWidth()[i] / dt;
    }
    return;
}


void MassBalance::setAlgorithmProperties(const json& properties){
    MAX_ITERATIONS = properties["massBalanceMaxIterations"];
    MIN_STAB_ITERATIONS = properties["massBalanceMinStabIterations"];
    TOLERANCE = properties["massBalanceTolerance"];
    MIN_MOBILITY_WIDTH = properties["massBalanceMinMobilityWidth"];
    return;
}


void MassBalance::calculateMobility(){
    int n = mesh->size();
    auto& mobility = new_dike->mobility;
    auto& width = new_dike->width;
    auto& viscosity = new_dike->viscosity;

    mobility.fill(0.0);
    for (int i = 1; i < n; ++i){
        double ml = width[i-1] <= MIN_MOBILITY_WIDTH ? 0.0 : std::pow(width[i-1], 3) / 12.0 / viscosity[i-1];
        double mr = width[i] <= MIN_MOBILITY_WIDTH ? 0.0 : std::pow(width[i], 3) / 12.0 / viscosity[i];
        mobility[i] = (ml + mr) / 2.0;
    }
    return;
}