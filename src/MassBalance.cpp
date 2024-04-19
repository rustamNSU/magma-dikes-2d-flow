#include "MassBalance.hpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stack>
#include "Utils.hpp"
#include <Eigen/Core>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using json = nlohmann::json;

MassBalance::MassBalance(
    InputData* input,
    TimestepController* timestep_controller,
    Elasticity* elasticity,
    Mesh* mesh,
    Schedule* schedule,
    MagmaState* magma_state,
    ReservoirData* reservoir
) : 
    input(input),
    timestep_controller(timestep_controller),
    elasticity(elasticity),
    mesh(mesh),
    schedule(schedule),
    magma_state(magma_state),
    reservoir(reservoir)
{
    algorithm_properties = input->getAlgorithmProperties();
}


void MassBalance::setNewTimestepData(
    DikeData* new_dike,
    DikeData* old_dike
){
    this->new_dike = new_dike;
    this->old_dike = old_dike;
    return;
}


MassBalance::explicitSolverOutput MassBalance::explicitSolve(){
    explicitSolverOutput output;
    auto nx = mesh->size();
    double current_time = timestep_controller->getCurrentTime();
    double dt = timestep_controller->getCurrentTimestep();
    updatePressure();
    magma_state->updateViscosity(new_dike);
    calculateMobility();
    auto& xc = mesh->getx();
    auto dx = mesh->getdx();
    auto g = reservoir->getGravityAcceleration();
    VectorXd dpdx = VectorXd::Zero(nx+1);
    for (int ix = 1; ix < nx; ix++){
        dpdx[ix] = (new_dike->pressure[ix] - new_dike->pressure[ix-1]) / (xc[ix] - xc[ix-1]) + 
                   (new_dike->density[ix] + new_dike->density[ix-1]) / 2 * g;
        if (dpdx[ix] >= 0){
            dpdx[ix] = 0;
        }
    }
    new_dike->total_face_flux = new_dike->total_face_mobility.cwiseProduct(dpdx);
    for (int ix = 1; ix < nx; ++ix){
        double Qf = new_dike->total_face_flux[ix];
        double Mf = Qf * dt;
        double Min = new_dike->width[ix-1] * dx;
        if (Mf > CFL_FACTOR * Min){
            output.cfl_condition = false;
            output.ratio = std::max(int(std::ceil(Mf / (CFL_FACTOR * Min))), output.ratio);
        }
        if (Mf > 2 * MIN_MOBILITY_WIDTH * dx){
            new_dike->face_flux.row(ix) = Qf * new_dike->mobility.row(ix-1) / new_dike->total_mobility[ix-1];
        }
        else{
            new_dike->total_face_flux[ix] = 0.0;
            new_dike->face_flux.row(ix).fill(0.0);
        }
    }
    if (output.cfl_condition == false){
        return output;
    }
    VectorXd width = VectorXd::Zero(nx);
    VectorXd Qin = VectorXd::Zero(nx);
    Qin[0] = schedule->getMassRate(current_time, current_time+dt) / new_dike->density[0];
    for (int ix = 0; ix < nx; ix++){
        width[ix] = old_dike->width[ix] - dt / dx * (
            new_dike->total_face_flux[ix+1] - new_dike->total_face_flux[ix] - Qin[ix]
        );
    }
    new_dike->setWidth(width);
    updateTemperature();
    return output;
}


void MassBalance::updatePressure(){
    new_dike->overpressure = elasticity->getMatrix() * new_dike->width;
    new_dike->pressure = new_dike->overpressure + reservoir->getLithostaticPressure();
    return;
}
// MassBalance::SolverOutput MassBalance::solve(){
//     int iter = 0;
//     int min_stab_iter = 0;
//     double err_max;
//     double err_l2;
//     SolverOutput solver_output;
//     std::vector<double> err_list;
//     err_list.reserve(MAX_ITERATIONS);
//     for (; iter < MAX_ITERATIONS; ++iter){
//         auto witer = new_dike->getWidth();
//         generateMatrix();
//         VectorXd sol = mat.fullPivLu().solve(rhs);
//         err_l2 = (witer - sol).norm() / std::max(1.0, witer.norm());
//         err_max = (witer - sol).maxCoeff() / std::max(TOLERANCE, witer.maxCoeff());
//         err_list.push_back(err_max);
//         new_dike->setWidth(sol);
//         magma_state->updateViscosity(new_dike);
//         new_dike->overpressure = elasticity->getMatrix() * new_dike->width;
//         new_dike->pressure = new_dike->overpressure + input->getPlith();
//         if (err_l2 < TOLERANCE){
//             ++min_stab_iter;
//         }
//         else{
//             min_stab_iter = 0;
//         }
//         err_list.push_back(err_l2);
//         if (min_stab_iter >= MIN_STAB_ITERATIONS){
//             solver_output.is_converge = true;
//             break;
//         }
//     }
//     solver_output.error = err_max;
//     solver_output.iters = iter;
//     return solver_output;
// }


// void MassBalance::generateMatrix(){
//     int nx = mesh->size();
//     int ny = new_dike->getLayersNumber();
//     double dx = mesh->getdx();
//     double g = input->getGravityAcceleration();
//     double Mch = schedule->getMassRate(old_time, new_time);
//     calculateMobility();
//     mat.resize(nx, nx);
//     rhs.resize(nx);
//     mat.fill(0.0);
//     rhs.fill(0.0);
//     auto C = elasticity->getMatrix();
//     auto plith = input->getLithostaticPressure();

//     rhs[0] += Mch / dx;
//     for (int ix = 1; ix < nx; ix++){
//         double lambda = new_dike->getTotalMobility()[ix];
//         double rho = (new_dike->getDensity()[ix] + new_dike->getDensity()[ix-1]) / 2.0;
//         auto c1 = -rho * lambda / dx / dx;
//         auto qrow = c1 * (C.row(ix).array() - C.row(ix-1).array());
//         mat.row(ix-1) += qrow.matrix();
//         mat.row(ix) -= qrow.matrix();

//         double qlith = -c1 * (plith[ix] - plith[ix-1]);
//         rhs[ix-1] += qlith;
//         rhs[ix] -= qlith;

//         double c2 = rho * lambda / dx * rho * g;
//         rhs[ix-1] += c2;
//         rhs[ix] -= c2;
//     }

//     for (int ix = 0; ix < nx; ix++){
//         mat(ix, ix) += new_dike->getDensity()[ix] / dt;
//         rhs[ix] += old_dike->getDensity()[ix] * old_dike->getWidth()[ix] / dt;
//     }
//     return;
// }


void MassBalance::setAlgorithmProperties(const json& properties){
    TIMESTEP_SCHEME = properties["timestepScheme"];
    if (TIMESTEP_SCHEME == "explicit"){
        MIN_MOBILITY_WIDTH = properties["massBalanceMinMobilityWidth"];
        CUTOFF_VELOCITY = properties["cutoffVelocity"];
        CFL_FACTOR = properties["lubricationCflFactor"];
    }
    // MAX_ITERATIONS = properties["massBalanceMaxIterations"];
    // MIN_STAB_ITERATIONS = properties["massBalanceMinStabIterations"];
    // TOLERANCE = properties["massBalanceTolerance"];
    return;
}


/* We define velocity Ui in i-th layer as Ui = (Ai*y^2 + Ci) * (dp/dx + rho_m*g) */
void MassBalance::calculateMobility(){
    int nx = mesh->size();
    int ny = new_dike->getLayersNumber();
    auto& width = new_dike->width;
    auto& viscosity = new_dike->viscosity;
    auto& mobility = new_dike->mobility;
    auto& total_mobility = new_dike->total_mobility;
    

    mobility.fill(0.0);
    total_mobility.fill(0.0);
    for (int ix = 0; ix < nx; ++ix){
        if (width[ix] < MIN_MOBILITY_WIDTH) continue;
        double dy = width[ix] / ny;
        VectorXd yb = new_dike->yb.row(ix);
        VectorXd yc = new_dike->yc.row(ix);
        VectorXd mu = viscosity.row(ix);
        VectorXd A = 0.5 * mu.cwiseInverse();
        VectorXd C = VectorXd::Zero(ny);
        C(ny-1) = -A(ny-1) * std::pow(yb(ny), 2);
        for (int iy = ny-2; iy >= 0; iy--){
            C(iy) = C(iy+1) + std::pow(yb(iy+1), 2) * (A(iy+1) - A(iy));
        }
        VectorXd ybt = yb(Eigen::seq(1, ny));
        VectorXd ybb = yb(Eigen::seq(0, ny-1));
        mobility.row(ix) = A.array() * (ybt.array().cube() - ybb.array().cube()) / 3 + C.array() * (ybt.array() - ybb.array());
        total_mobility(ix) = 2 * mobility.row(ix).sum();
    }
    new_dike->total_face_mobility.setZero();
    new_dike->total_face_mobility(Eigen::seq(1, nx-1)) = 0.5 * (total_mobility(Eigen::seq(0, nx-2)) + total_mobility(Eigen::seq(1, nx-1)));
    return;
}


void MassBalance::updateTemperature(){
    double dt = timestep_controller->getCurrentTimestep();
    double dx = mesh->getdx();
    int nx = mesh->size();
    int ny = new_dike->ny;

    new_dike->temperature.row(0).fill(schedule->getMagmaChamberTemperature());
    auto Tconv = new_dike->temperature;
    for (int ix = 1; ix < nx; ix++){
        VectorXd Tin = old_dike->temperature.row(ix-1);
        VectorXd Vin = new_dike->face_flux.row(ix) * dt;
        VectorXd Vout = new_dike->face_flux.row(ix+1) * dt;

        VectorXd Told = old_dike->temperature.row(ix);
        VectorXd Vold = (old_dike->yb(ix, Eigen::seq(1, ny)) - old_dike->yb(ix, Eigen::seq(0, ny-1))) * dx;
        VectorXd Tnew = Told;
        VectorXd Vnew = (new_dike->yb(ix, Eigen::seq(1, ny)) - new_dike->yb(ix, Eigen::seq(0, ny-1))) * dx;
        if (Vnew.sum() < MIN_MOBILITY_WIDTH) continue;

        VectorXd Vcfl = Vold - Vout;
        double Vmin = Vcfl.minCoeff();
        if (Vmin < 0){
            std::cout << "-- Vcfl < 0\n";
        }
        Vold = Vold - Vout;

        VectorXd Vall;
        VectorXd Tall;
        // if (Vold.sum() < MIN_MOBILITY_WIDTH){
        //     Vall = Vin;
        //     Tall = Tin;
        // }
        // else{
            Vall.resize(2*ny);
            Tall.resize(2*ny);
            Vall << Vin, Vold;
            Tall << Tin, Told;
        // }

        VectorXi indxs = VectorXi::LinSpaced(Vall.size(), 0, Vall.size()-1);
        std::stable_sort(indxs.begin(), indxs.end(), [&Tall](
            int i1, int i2
        ){
            return Tall[i1] > Tall[i2];
        });

        Vall(indxs);
        Tall(indxs);
        std::stack<double> Vstack(std::deque<double>(Vall.begin(), Vall.end()));
        std::stack<double> Tstack(std::deque<double>(Tall.begin(), Tall.end()));

        int iy = ny-1;
        double V = Vnew[iy];
        double E = 0;
        while (!Vstack.empty()){
            double Vtmp = Vstack.top();
            double Ttmp = Tstack.top();
            Vstack.pop();
            Tstack.pop();
            if (V > Vtmp || iy == 0){
                V = V - Vtmp;
                E += Vtmp * Ttmp;
            }
            else{
                E += V * Ttmp;
                Tnew[iy] = E / Vnew[iy];
                Vstack.push(Vtmp - V);
                Tstack.push(Ttmp);
                iy -= 1;
                V = Vnew[iy];
                E = 0;
            }
        }
        Tnew[0] = E / Vnew[0];
        Tconv.row(ix) = Tnew;
        // Tnew = Tnew.array() - 1000.0;
        // if (std::abs(Tnew.maxCoeff() - Tnew.minCoeff()) > 1e-6){
        //     int aaa = 0;
        // }
    }
    new_dike->temperature = Tconv;

    double rhom = schedule->getMagmaChamberDensity();
    double Cm = magma_state->getSpecificHeat();
    double km = magma_state->getThermalConductivity();
    double Cr = reservoir->C;
    double kr = reservoir->k;
    auto reservoir_boundary_temperature = reservoir->getInitialTemperature();
    auto reservoir_temperature = reservoir->temperature;
    for (int ix = 1; ix < nx; ix++){
        if (new_dike->width[ix] < MIN_MOBILITY_WIDTH) continue;
        VectorXd yb = new_dike->yb.row(ix);
        VectorXd yc = new_dike->yc.row(ix);
        VectorXd dy = yb(Eigen::seq(1, ny)) - yb(Eigen::seq(0, ny-1));
        VectorXd Told = Tconv.row(ix);
        VectorXd Tnew = Told;

        int nr = reservoir->ny;
        double Tinf = reservoir_boundary_temperature[ix];
        VectorXd Trold = reservoir_temperature.row(ix);
        VectorXd Trnew = Trold;
        const auto& ycr = reservoir->yc;
        const auto& ybr = reservoir->yb;
        const auto& dyr = reservoir->dy;
        double rhor = reservoir->density[ix];

        std::vector<double> a(ny+nr), b(ny+nr), c(ny+nr), rhs(ny+nr);
        if (ny == 1){
            double coef = dy[0] * dy[0] * rhom * Cm / dt;
            double ktop = km*kr*dy[0] / (km*ycr[0] + kr*(yb[1] - yc[0]));
            a[0] = 0.0;
            b[0] = coef + ktop;
            c[0] = -ktop;
            rhs[0] = coef * Told[0];
        }
        else{
            double coef = dy[0] * dy[0] * rhom * Cm / dt;
            a[0] = 0.0;
            b[0] = coef + km;
            c[0] = -km;
            rhs[0] = coef * Told[0];
            for (int iy = 1; iy < ny-1; iy++){
                coef = dy[iy] * dy[iy] * rhom * Cm / dt;
                a[iy] = -km;
                b[iy] = coef + 2*km;
                c[iy] = -km;
                rhs[iy] = coef * Told[iy];
            }
            coef = dy[ny-1] * dy[ny-1] * rhom * Cm / dt;
            double coef2 = km*kr*dy[ny-1] / (km*ycr[0] + kr*(yb[ny] - yc[ny-1]));
            a[ny-1] = -km;
            b[ny-1] = coef + km + coef2;
            c[ny-1] = -coef2;
            rhs[ny-1] = coef*Told[ny-1];
        }

        double coef = dyr[0] * dyr[0] * rhor * Cr / dt;
        double kbot = km*kr*dyr[0] / (km*ycr[0] + kr*(yb[ny] - yc[ny-1]));
        double ktop = kr * dyr[0] / (ycr[1] - ycr[0]);
        a[ny] = -kbot;
        b[ny] = coef + ktop + kbot;
        c[ny] = -ktop;
        rhs[ny] = coef * Trold[0];
        for (int ir = 1; ir < nr - 1; ir++){
            coef = dyr[ir] * dyr[ir] * rhor * Cr / dt;
            kbot = kr * dyr[ir] / (ycr[ir] - ycr[ir-1]);
            ktop = kr * dyr[ir] / (ycr[ir+1] - ycr[ir]);
            a[ny+ir] = -kbot;
            b[ny+ir] = coef + ktop + kbot;
            c[ny+ir] = -ktop;
            rhs[ny+ir] = coef * Trold[ir];
        }
        coef = dyr[nr-1] * dyr[nr-1] * rhor * Cr / dt;
        kbot = kr * dyr[nr-1] / (ycr[nr-1] - ycr[nr-2]);
        ktop = kr * dyr[nr-1] / (ybr[nr] - ycr[nr-1]);
        a[ny+nr-1] = -kbot;
        b[ny+nr-1] = coef + ktop + kbot;
        c[ny+nr-1] = 0.0;
        rhs[ny+nr-1] = coef * Trold[nr-1] + ktop * Tinf;
        auto sol = Utils::tridiagonal_solver(a, b, c, rhs);
        Tnew = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sol.data(), ny);
        Trnew = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sol.data() + ny, nr);
        new_dike->temperature.row(ix) = Tnew;
        reservoir->temperature.row(ix) = Trnew;
    }
}


void MassBalance::setInitialData(){
    if (new_dike->model == FlowModel::channel){
        double Tch = schedule->getMagmaChamberTemperature();
        new_dike->temperature.fill(Tch);
        new_dike->pressure = reservoir->getLithostaticPressure();
        magma_state->updateViscosity(new_dike);
    }
}