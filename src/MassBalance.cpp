#include "MassBalance.hpp"
#include <cmath>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
        new_dike->pressure = elasticity->getMatrix() * new_dike->width + input->getPlith();
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
    return;
}


void MassBalance::generateMatrix(){
    int n = mesh->size();
    double dx = mesh->getdx();
    double Mch = schedule->getMassRate(old_time, new_time);
    new_dike->calculateMobility();
    mat.resize(n, n);
    rhs.resize(n);
    mat.fill(0.0);
    rhs.fill(0.0);
    auto C = elasticity->getMatrix();
    auto plith = input->getPlith();

    rhs[0] += Mch / dx;
    for (int i = 1; i < n; i++){
        double lambda = new_dike->getMobility()[i];
        auto c1 = -lambda / dx / dx;
        auto qrow = c1 * (C.row(i).array() - C.row(i-1).array());
        mat.row(i-1) += qrow.matrix();
        mat.row(i) -= qrow.matrix();

        double qlith = -c1 * (plith[i] - plith[i-1]);
        rhs[i-1] += qlith;
        rhs[i] -= qlith;

        double rho = (new_dike->getDensity()[i] + new_dike->getDensity()[i-1]) / 2.0;
        double c2 = lambda / dx * rho * input->getg();
        rhs[i-1] += c2;
        rhs[i] -= c2;
    }

    for (int i = 0; i < n; i++){
        mat(i, i) += new_dike->getDensity()[i] / dt;
        rhs[i] += old_dike->getDensity()[i] * old_dike->getWidth()[i] / dt;
    }
    return;
}