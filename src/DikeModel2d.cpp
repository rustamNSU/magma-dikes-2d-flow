#include "DikeModel2d.hpp"
#include <fstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stack>
#include "Utils.hpp"
#include <Eigen/Core>
#include <highfive/H5Easy.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::VectorXi;
using json = nlohmann::json;
using H5Easy::File;
using H5Easy::dump;

DikeModel2d::DikeModel2d(const std::string& input_path) :
    input_path(input_path)
{
    std::ifstream f(input_path);
    json input_json = json::parse(f);
	f.close();
    input = std::make_shared<InputData>(input_json);
    timestep_controller = std::make_shared<TimestepController>(input->getTimestepProperties());
    mesh = std::make_shared<Mesh>(input->getMeshProperties());
    reservoir = std::make_shared<ReservoirData>(mesh.get(), input->getReservoirProperties());
    schedule = std::make_shared<Schedule>(mesh.get(), input->getScheduleProperties());
    magma_state = std::make_shared<MagmaState>(mesh.get(), input->getMagmaProperties());
    algorithm_properties = input->getAlgorithmProperties();
    dike = std::make_shared<DikeData>(mesh.get(), algorithm_properties);

    auto [E, nu, KIc] = reservoir->getElasticityParameters();
    elasticity = std::make_shared<Elasticity>(E, nu, mesh.get());
    setAlgorithmProperties();
    setInitialData();
}


void DikeModel2d::setAlgorithmProperties(){
    // TIMESTEP_SCHEME = properties["timestepScheme"];
    MIN_MOBILITY_WIDTH = algorithm_properties["massBalanceMinMobilityWidth"];
    CUTOFF_VELOCITY = algorithm_properties["cutoffVelocity"];
    CFL_FACTOR = algorithm_properties["lubricationCflFactor"];
    // MAX_ITERATIONS = properties["massBalanceMaxIterations"];
    // MIN_STAB_ITERATIONS = properties["massBalanceMinStabIterations"];
    // TOLERANCE = properties["massBalanceTolerance"];
    return;
}


void DikeModel2d::setInitialData(){
    auto plith = reservoir->getLithostaticPressure();
    auto tlith = reservoir->getInitialTemperature();
    dike->setInitialPressure(plith);
    dike->setInitialTemperature(tlith);
    magma_state->updateDensity(dike.get());
    magma_state->updateViscosity(dike.get());
    dike->time = timestep_controller->getCurrentTime();
    old_dike = std::make_shared<DikeData>(*dike);

    auto [is_save, save_timestep] = timestep_controller->saveTimestepIteration();
    std::string savepath = (input->getDataDir() / "data_").string() + std::to_string(save_timestep) + ".h5";
    saveData(savepath);
}


void DikeModel2d::run(){
    while (!timestep_controller->isFinish()){
        auto time = timestep_controller->getCurrentTime();
        auto dt = timestep_controller->getCurrentTimestep();
        solver_log.setDefault();
        explicitSolver();
        if (solver_log.successful){
            timestep_controller->update();
            dike->time = time + dt;
            *old_dike = *dike;
            auto [is_save, save_timestep] = timestep_controller->saveTimestepIteration();
            if (is_save){
                spdlog::info("{:^7} | {:<5} h | dt = {:<5}",
                    save_timestep,
                    time / 3600.0,
                    dt
                );
                std::string savepath = (input->getDataDir() / "data_").string() + std::to_string(save_timestep) + ".h5";
                saveData(savepath);
            }
        }
        else{
            *dike = *old_dike;
            timestep_controller->divideTimestep(solver_log.cfl_ratio);
        }
    }
}


void DikeModel2d::explicitSolver(){
    auto nx = mesh->size();
    auto ny = dike->getLayersNumber();
    auto dt = timestep_controller->getCurrentTimestep();

    magma_state->updateViscosity(dike.get());
    updatePressure();
    calculateVerticalFlow();
    solveMassBalance();
    if (!solver_log.successful){
        return;
    }
    solveEnergyBalance();
}


void DikeModel2d::updatePressure(){
    dike->overpressure = 2 * (elasticity->getMatrix() * old_dike->hw);
    dike->pressure = old_dike->overpressure + reservoir->getLithostaticPressure().matrix();
}


void DikeModel2d::calculateVerticalFlow(){
    int nx = mesh->size();
    int ny = dike->getLayersNumber();
    const auto& x = mesh->getx();
    auto& hw = old_dike->hw;
    auto& viscosity = dike->viscosity;
    auto& qx = dike->qx;
    auto& Qx = dike->Qx;
    const auto& yb = dike->yb;
    const auto& yc = dike->yc;
    double dy = yb[1] - yb[0];
    double g = reservoir->getGravityAcceleration();

    /* qx[ix] - flow beetween (ix-1) and ix elements */
    for (int ix = 1; ix < nx; ++ix){
        if (std::max(hw[ix-1], hw[ix]) < MIN_MOBILITY_WIDTH) continue;
        double h = 0.5 * (hw[ix-1] + hw[ix]);
        double xl = x[ix-1];
        double xr = x[ix];
        double dpdx = (dike->pressure[ix] - dike->pressure[ix-1]) / (xr - xl);
        
        // @todo: add average density
        double rho = dike->density(ix, 0);
        double G = dpdx + rho * g;
        if (G > 0){
            G = 0;
        }
        dike->G[ix] = G;
        ArrayXd mul = viscosity.row(ix-1);
        ArrayXd mur = viscosity.row(ix);
        ArrayXd mu = 2 * (mul * mur) / (mul + mur);
        ArrayXd A = 0.5*h*h * mu.cwiseInverse();
        ArrayXd C = ArrayXd::Zero(ny);
        C(ny-1) = -A(ny-1);
        for (int iy = ny-2; iy >= 0; iy--){
            C(iy) = C(iy+1) + yb(iy+1) * yb(iy+1) * (A(iy+1) - A(iy));
        }
        ArrayXd ybt = yb(Eigen::seq(1, ny));
        ArrayXd ybb = yb(Eigen::seq(0, ny-1));
        qx.row(ix) = (A * (ybt.cube() - ybb.cube()) / 3 + C * (ybt - ybb)) * h * G;
        Qx(ix) = qx.row(ix).sum();
        dike->A.row(ix) = A * G;
        dike->C.row(ix) = C * G;
    }
    return;
}


void DikeModel2d::solveMassBalance(){
    int nx = mesh->size();
    int ny = dike->getLayersNumber();
    auto& h = dike->hw;
    const auto& hold = old_dike->hw;
    double dx = mesh->getdx();
    double dt = timestep_controller->getCurrentTimestep();
    double t1 = timestep_controller->getCurrentTime();

    const auto& Qx = dike->Qx;
    double Q0 = 0.5 * schedule->getMassRate(t1, t1 + dt) / schedule->getMagmaChamberDensity();
    VectorXd Qin = VectorXd::Zero(nx);
    Qin[0] = Q0;
    h = hold + dt/dx*(Qin + Qx(Eigen::seq(0, nx-1)) - Qx(Eigen::seq(1, nx)));
    // ArrayXd cfl_array = CFL_FACTOR * dike->getElementsVolume() - dt * Qx(Eigen::seq(1, nx)).array();
    // if (cfl_array.minCoeff() < 1e-12){
        
    // }
    ArrayXd Vold = old_dike->getElementsVolume();
    for (int ix = 1; ix < nx; ++ix){
        double Vout = Qx[ix] * dt;
        double Vtmp = Vold[ix-1];
        if (Vout > CFL_FACTOR * Vtmp){
            solver_log.successful = false;
            solver_log.cfl_condition = false;
            solver_log.cfl_ratio = std::max(int(std::ceil(Vout / (CFL_FACTOR * Vtmp))), solver_log.cfl_ratio);
        }
    }
    if (!solver_log.successful){
        return;
    }

    auto &qy = dike->qy;
    const auto &qx = dike->qx;
    const auto& yb = dike->yb;
    const auto& yc = dike->yc;
    double dy = yb[1] - yb[0];
    for (int ix = 0; ix < nx; ix++){
        if (h[ix] < MIN_MOBILITY_WIDTH) continue;
        for (int iy = 1; iy < ny; iy++){
            qy(ix, iy) = qy(ix, iy-1) - qx(ix+1, iy-1) + qx(ix, iy-1) - dy*dx/dt*(h[ix] - hold[ix]) + Qin[ix]*dy;
        }
        /* check convergence: must be zero */
        double err = dy*dx/dt*(h[ix] - hold[ix]) - Qin[ix]*dy + qx(ix+1, ny-1) - qx(ix, ny-1) - qy(ix, ny-1);
        double a = 0;
    }
}


void DikeModel2d::solveEnergyBalance(){
    int nx = mesh->size();
    int ny = dike->getLayersNumber();
    const auto& h = dike->hw;
    const auto& hold = old_dike->hw;
    double dx = mesh->getdx();
    double dt = timestep_controller->getCurrentTimestep();
    double t1 = timestep_controller->getCurrentTime();

    const auto& qx = dike->qx;
    const auto& qy = dike->qy;
    const auto& yb = dike->yb;
    const auto& yc = dike->yc;
    double dy = yb[1] - yb[0];

    int nr = reservoir->ny;
    double Cm = magma_state->getSpecificHeat();
    double km = magma_state->getThermalConductivity();

    double Cr = reservoir->C;
    double kr = reservoir->k;
    auto Trinf = reservoir->getInitialTemperature();
    const auto& ycr = reservoir->yc;
    const auto& ybr = reservoir->yb;
    const auto& dyr = reservoir->dy;

    double E0 = 0.5 * schedule->getMassRate(t1, t1 + dt) * schedule->getMagmaChamberTemperature() * Cm;
    int N = ny+nr+1;
    for (int ix = 0; ix < nx; ix++){
        if (h[ix] < MIN_MOBILITY_WIDTH) continue;
        std::vector<double> a(N, 0.0), b(N, 0.0), c(N, 0.0), rhs(N, 0.0);
        // ArrayXd rhom = dike->density.row(ix); // @todo
        double rhom = dike->density(0, 0);
        const ArrayXd Told = old_dike->temperature.row(ix);
        ArrayXd Ein = ArrayXd::Zero(ny);
        if (ix == 0){
            Ein.fill(E0*dy);
        }
        else{
            Ein = rhom * Cm * dike->temperature.row(ix-1).array() * qx.row(ix).array();
        }
        double dyt, dyb, vtp, vtm, vbp, vbm;
        dyt = yc[1] - yc[0];
        vtp = std::max(qy(ix, 1), 0.0); // qy in top
        vtm = -std::max(-qy(ix, 1), 0.0); // 0 if qy > 0 in top
        a[0] = 0;
        b[0] = rhom*Cm*(h[ix]*dx*dy/dt + qx(ix+1, 0) + vtp) + dx*km/h[ix]*(1.0/dyt);
        c[0] = rhom*Cm*vtm - dx*km/h[ix]*(1.0/dyt);
        rhs[0] = rhom*Cm*Told[0]*hold[ix]*dx*dy/dt + Ein[0];

        for (int iy = 1; iy < ny-1; iy++){
            dyt = yc[iy+1] - yc[iy];
            dyb = yc[iy] - yc[iy-1];
            vtp = std::max(qy(ix, iy+1), 0.0); // qy in top
            vtm = -std::max(-qy(ix, iy+1), 0.0); // 0 if qy > 0 in top
            vbp = std::max(qy(ix, iy), 0.0); // qy in bot
            vbm = -std::max(-qy(ix, iy), 0.0); // 0 if qy > 0 in bot
            a[iy] = -rhom*Cm*vbp - dx*km/h[ix]/dyb;
            b[iy] = rhom*Cm*(h[ix]*dx*dy/dt + qx(ix+1, iy) + vtp - vbm) + dx*km/h[ix]*(1.0/dyt + 1.0/dyb);
            c[iy] = rhom*Cm*vtm - dx*km/h[ix]/dyt;
            rhs[iy] = rhom*Cm*Told[iy]*hold[ix]*dx*dy/dt + Ein[iy];
        }
        dyt = yb[ny] - yc[ny-1];
        dyb = yc[ny-1] - yc[ny-2];
        vbp = std::max(qy(ix, ny-1), 0.0); // qy in bot
        vbm = -std::max(-qy(ix, ny-1), 0.0); // 0 if qy > 0 in bot
        a[ny-1] = -rhom*Cm*vbp - dx*km/h[ix]/dyb;
        b[ny-1] = rhom*Cm*(h[ix]*dx*dy/dt + qx(ix+1, ny-1) - vbm) + dx*km/h[ix]*(1.0/dyt + 1.0/dyb);
        c[ny-1] = -dx*km/h[ix]/dyt;
        rhs[ny-1] = rhom*Cm*Told[ny-1]*hold[ix]*dx*dy/dt + Ein[ny-1];
        a[ny] = -dx*km/h[ix]/dyt;
        b[ny] = dx*km/h[ix]/dyt;

        /* Reservoir conduction */
        int ind = ny+1;
        double rhor = reservoir->density[ix];
        auto Tr = reservoir->temperature.row(ix);
        double ktop, kbot;
        dyt = ycr[1] - ycr[0];
        dyb = ycr[0] - ybr[0];
        b[ny] += kr*dx/dyb; // on wall
        c[ny] = -kr*dx/dyb; // on wall
        a[ind] = -kr*dx/dyb;
        b[ind] = rhor*Cr*dx*dyr[0]/dt + kr*dx/dyt + kr*dx/dyb;
        c[ind] = -kr*dx/dyt;
        rhs[ind] = rhor*Cr*Tr[0]*dx*dyr[0]/dt;
        for (int ir = 1; ir < nr-1; ir++){
            ind++;
            dyt = ycr[ir+1] - ycr[ir];
            dyb = ycr[ir] - ycr[ir-1];
            a[ind] = -kr*dx/dyb;
            b[ind] = rhor*Cr*dx*dyr[ir]/dt + kr*dx/dyt + kr*dx/dyb;
            c[ind] = -kr*dx/dyt;
            rhs[ind] = rhor*Cr*Tr[ir]*dx*dyr[ir]/dt;
        }
        ind++;
        dyt = ybr[nr] - ycr[nr-1];
        dyb = ycr[nr-1] - ycr[nr-2];
        a[ind] = -kr*dx/dyb;
        b[ind] = rhor*Cr*dx*dyr[nr-1]/dt + kr*dx/dyt + kr*dx/dyb;
        c[ind] = 0.0;
        rhs[ind] = rhor*Cr*Tr[nr-1]*dx*dyr[nr-1]/dt + kr*dx/dyt*Trinf[ix];

        /* Solve tridiagonal system */
        auto sol = Utils::tridiagonal_solver(a, b, c, rhs);
        dike->temperature.row(ix) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sol.data(), ny);
        dike->Twall[ix] = sol[ny];
        reservoir->temperature.row(ix) =  Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sol.data() + ny+1, nr);
    }
}


void DikeModel2d::saveData(const std::string &savepath){
    File file(savepath, File::Overwrite);
    dump(file, "time", dike->time);
    dump(file, "mesh/xc", dike->meshX->getx());
    dump(file, "mesh/xl", dike->meshX->getxl());
    dump(file, "mesh/xr", dike->meshX->getxr());

    dump(file, "ny", dike->ny);
    dump(file, "yc", dike->yc);
    dump(file, "yb", dike->yb);
    dump(file, "halfwidth", dike->hw);
    dump(file, "width", dike->getWidth());
    dump(file, "pressure", dike->pressure);
    dump(file, "overpressure", dike->overpressure);
    dump(file, "density", dike->density);
    dump(file, "viscosity", dike->viscosity);
    dump(file, "temperature", dike->temperature);
    dump(file, "Twall", dike->Twall);
    dump(file, "qx", dike->qx);
    dump(file, "qy", dike->qy);
    dump(file, "A", dike->A);
    dump(file, "C", dike->C);
    dump(file, "G", dike->G);

    dump(file, "reservoir/yc", reservoir->yc);
    dump(file, "reservoir/yb", reservoir->yb);
    dump(file, "reservoir/temperature", reservoir->temperature);
    dump(file, "reservoir/ny", reservoir->ny);
}