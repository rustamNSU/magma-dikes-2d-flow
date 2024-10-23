#include "DikeModel2d.hpp"
#include <fstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stack>
#include <exception>
#include "Utils.hpp"
#include <Eigen/Core>
#include <highfive/H5Easy.hpp>
#include "Interp.hpp"
#include "LinearSolver.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
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
    MIN_MOBILITY_WIDTH = algorithm_properties["massBalanceMinMobilityWidth"];
    CUTOFF_VELOCITY = algorithm_properties["cutoffVelocity"];
    CFL_FACTOR = algorithm_properties["lubricationCflFactor"];
    VISCOSITY_APPROXIMATION = algorithm_properties["viscosityApproximation"].get<std::string>();
    shearHeating = algorithm_properties["shearHeating"].get<bool>() ? 1.0 : 0.0;
    latentHeatCrystallization = algorithm_properties["latentHeatCrystallization"].get<bool>() ? 1.0 : 0.0;
    highOrderApproximation = algorithm_properties["highOrderApproximation"].get<bool>();


    if (algorithm_properties["isCohesiveStress"].get<bool>()){
        auto [E, nu, K1c] = reservoir->getElasticityParameters();
        cohesive_model = std::make_shared<CohesiveModel>(mesh.get(), E, nu, K1c);
    }
    return;
}


void DikeModel2d::setInitialData(){
    auto plith = reservoir->getLithostaticPressure();
    auto tlith = reservoir->getInitialTemperature();
    magma_state->setChamberInitialState(plith[0], tlith[0]);
    schedule->setMagmaChamberCrystallization(magma_state->chamber.beta);
    dike->setInitialPressure(plith);
    dike->setInitialTemperature(tlith);
    magma_state->updateEquilibriumCrystallization(dike.get());
    dike->beta.setZero();
    dike->beta.row(0) = magma_state->chamber.beta;
    dike->gamma.row(0) = magma_state->chamber.gamma;
    magma_state->updateDensity(dike.get());
    magma_state->updateViscosity(dike.get());
    dike->time = timestep_controller->getCurrentTime();
    old_dike = std::make_shared<DikeData>(*dike);

    auto [is_save, save_timestep] = timestep_controller->saveTimestepIteration();
    std::string savepath = (input->getDataDir() / "data_").string() + std::to_string(save_timestep) + ".h5";
    saveData(savepath);
}


void DikeModel2d::run(){
    JUST_TIMER_START("0) Total run time");
    while (!timestep_controller->isFinish()){
        auto time = timestep_controller->getCurrentTime();
        auto dt = timestep_controller->getCurrentTimestep();
        solver_log.setDefault();
        if (algorithm_properties["timestepScheme"] == "explicit"){
            JUST_TIMER_START("1) Explicit solver time");
            explicitSolver();
            JUST_TIMER_STOP("1) Explicit solver time");
        }
        else if (algorithm_properties["timestepScheme"] == "implicit"){
            JUST_TIMER_START("1) Implicit solver time");
            implicitSolver();
            JUST_TIMER_STOP("1) Implicit solver time");
        }
        else{
            throw std::invalid_argument("bad!\n");
        }

        JUST_TIMER_START("1) Update data time");
        if (solver_log.successful){
            updateData();
        }
        else{
            reloadData();
        }
        JUST_TIMER_STOP("1) Update data time");
    }
    JUST_TIMER_STOP("0) Total run time");
    // JUST_TIMER_PRINT_RATIO(std::cout, "0) Total run time");
    std::stringstream ss;
    std::string token;
    JUST_TIMER_PRINT_RATIO(ss, "0) Total run time");
    while (std::getline(ss, token, '\n')){
        spdlog::info(token);
    }
    JUST_TIMER_CLEAR;
}


void DikeModel2d::explicitSolver(){
    auto nx = mesh->size();
    auto ny = dike->getLayersNumber();
    auto dt = timestep_controller->getCurrentTimestep();
    JUST_TIMER_START("2) Update pressure time");
    updatePressure();
    JUST_TIMER_STOP("2) Update pressure time");
    JUST_TIMER_START("2) Update viscosity time");
    magma_state->updateViscosity(dike.get());
    JUST_TIMER_STOP("2) Update viscosity time");
    JUST_TIMER_START("2) Update density time");
    magma_state->updateDensity(dike.get());
    JUST_TIMER_STOP("2) Update density time");
    dike->setMagmaStateAfterTip();
    calculateVerticalFlow();
    solveMassBalance();
    if (!solver_log.successful){
        return;
    }
    dike->updateOpenElements();
    JUST_TIMER_START("2) Update equilibrium crystallization");
    magma_state->updateEquilibriumCrystallization(dike.get());
    JUST_TIMER_STOP("2) Update equilibrium crystallization");
    updateCrystallization();
    dike->setMagmaStateAfterTip();
    solveEnergyBalance();
}


void DikeModel2d::implicitSolver(){
    auto nx = mesh->size();
    auto ny = dike->getLayersNumber();
    auto dt = timestep_controller->getCurrentTimestep();
    JUST_TIMER_START("2) Update magma state");
    magma_state->updateViscosity(dike.get());
    magma_state->updateDensity(dike.get());
    magma_state->updateEquilibriumCrystallization(dike.get());
    dike->setMagmaStateAfterTip();
    JUST_TIMER_STOP("2) Update magma state");

    JUST_TIMER_START("2) Update magma width and pressure");
    implicitMassBalance();
    JUST_TIMER_STOP("2) Update magma width and pressure");
    updateCrystallization();
    solveEnergyBalance();
}


void DikeModel2d::updatePressure(){
    JUST_TIMER_START("2) Elasticity time");
    dike->overpressure = 2 * (elasticity->getMatrix() * dike->hw.matrix());
    // dike->overpressure = dike->overpressure.unaryExpr([&](double x){
    //     return 0.5*(x + std::abs(x));
    // });
    dike->pressure = dike->overpressure + reservoir->getLithostaticPressure();
    JUST_TIMER_STOP("2) Elasticity time");
}


void DikeModel2d::calculateVerticalFlow(){
    JUST_TIMER_START("2) Vert. flow calc. time");
    int nx = mesh->size();
    int ny = dike->getLayersNumber();
    const auto& x = mesh->getx();
    auto& hw = old_dike->hw;
    auto& viscosity = dike->viscosity;
    const auto& rho = dike->density;
    auto& qx = dike->qx;
    auto& Qx = dike->Qx;
    auto& Mx = dike->Mx;
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
        double dpdx;
        if (ix > 1 && highOrderApproximation){
            dpdx = parab_derivative(0.5*(xl+xr), {x[ix-2], x[ix-1], x[ix]}, {dike->pressure[ix-2], dike->pressure[ix-1], dike->pressure[ix]});
        }
        else{
            dpdx = (dike->pressure[ix] - dike->pressure[ix-1]) / (xr - xl);
        }
        
        // @todo: add average density
        double rho_avg = 0.5 * (rho.row(ix-1).mean() + rho.row(ix).mean());
        double G = dpdx + rho_avg * g;
        if (G > 0){
            G = 0;
        }
        // double GMAX = -2000*g;
        // G = std::max(G, GMAX);
        dike->G[ix] = G;
        ArrayXd mul = viscosity.row(ix-1);
        ArrayXd mur = viscosity.row(ix);
        ArrayXd mu = mul;
        if (VISCOSITY_APPROXIMATION == "harmonic"){
            mu = 2 * (mul * mur) / (mul + mur);
        }
        if (VISCOSITY_APPROXIMATION == "mean"){
            mu = 0.5 * (mul + mur);
        }
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
        Mx(ix) = (qx.row(ix).cwiseProduct(rho.row(ix-1))).sum();
        dike->A.row(ix) = A * G;
        dike->C.row(ix) = C * G;
    }

    double dx = mesh->getdx();
    for (int ix = 1; ix < nx; ++ix){
        if (hw[ix] < MIN_MOBILITY_WIDTH) continue;
        double G = 0.5 * (dike->G[ix] + dike->G[ix+1]);
        double h = hw[ix];
        ArrayXd ybt = yb(Eigen::seq(1, ny));
        ArrayXd ybb = yb(Eigen::seq(0, ny-1));
        dike->shear_heat.row(ix) = dx*(h*h*h*G*G/3) * (ybt.cube() - ybb.cube()).transpose() / dike->viscosity.row(ix).array();
    }
    JUST_TIMER_STOP("2) Vert. flow calc. time");
    return;
}


void DikeModel2d::solveMassBalance(){
    JUST_TIMER_START("2) Mass balance time");
    int nx = mesh->size();
    int ny = dike->getLayersNumber();
    auto& h = dike->hw;
    const auto& hold = old_dike->hw;
    const auto& rho_new = dike->density;
    const auto& rho_old = old_dike->density;
    ArrayXd rho_avg = dike->density.rowwise().mean();
    ArrayXd rho_old_avg = old_dike->density.rowwise().mean();
    double dx = mesh->getdx();
    double dt = timestep_controller->getCurrentTimestep();
    double t1 = timestep_controller->getCurrentTime();
    const auto& Qx = dike->Qx;
    const auto& Mx = dike->Mx;
    double M0 = 0.5 * schedule->getMassRate(t1, t1 + dt);
    ArrayXd Min = ArrayXd::Zero(nx);
    Min[0] = M0;
    h = rho_old_avg / rho_avg * hold + dt/dx*(Min + Mx(Eigen::seq(0, nx-1)) - Mx(Eigen::seq(1, nx))) / rho_avg;
    ArrayXd Vold = old_dike->getElementsVolume();
    for (int ix = 1; ix < nx; ++ix){
        double Vout = Qx[ix] * dt;
        double Vtmp = Vold[ix-1];
        if (Vout > CFL_FACTOR * Vtmp){
            solver_log.successful = false;
            solver_log.cfl_condition = false;
            solver_log.cfl_ratio = std::max(int(std::ceil(Vout / (CFL_FACTOR * Vtmp))), solver_log.cfl_ratio);
            // if (solver_log.cfl_ratio > 1000){
            //     double a = std::sin(solver_log.cfl_ratio);
            //     solver_log.cfl_ratio = std::max(0, int(a));
            // }
        }
    }
    if (!solver_log.successful){
        JUST_TIMER_STOP("2) Mass balance time");
        return;
    }

    auto &qy = dike->qy;
    const auto &qx = dike->qx; 
    const auto& yb = dike->yb;
    const auto& yc = dike->yc;
    double dy = yb[1] - yb[0];
    for (int ix = 0; ix < nx; ix++){
        if (h[ix] < MIN_MOBILITY_WIDTH) continue;
        ArrayXd rho_jb = 0.5 * (rho_new.row(ix)(Eigen::seq(1, ny-1)) + rho_new.row(ix)(Eigen::seq(0, ny-2)));
        ArrayXd rhoj = rho_new.row(ix);
        rhoj(Eigen::seq(1, ny-1)) = rho_jb;
        for (int iy = 1; iy < ny; iy++){
            qy(ix, iy) = (
                rhoj(iy-1) * qy(ix, iy-1) - 
                rho_new(ix, iy-1) * qx(ix+1, iy-1) + 
                rho_new(std::max(0, ix-1), iy-1) * qx(ix, iy-1) - 
                dy*dx/dt*(
                    rho_new(ix, iy-1)*h[ix] - rho_old(ix, iy-1)*hold[ix]
                ) + 
                Min[ix]*dy
            ) / rhoj(iy);
        }
        /* check convergence: must be zero */
        // double err = dy*dx/dt*(h[ix] - hold[ix]) - Qin[ix]*dy + qx(ix+1, ny-1) - qx(ix, ny-1) - qy(ix, ny-1);
        double a = 0;
    }
    JUST_TIMER_STOP("2) Mass balance time");
}


void DikeModel2d::solveEnergyBalance(){
    JUST_TIMER_START("2) Energy balance time");
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

    double rhoc0 = magma_state->density_properties["rhoc0"].get<double>();
    const auto& beta = dike->beta;
    const auto& alpha = dike->alpha;
    const auto& betaeq = dike->betaeq;
    double Lm = magma_state->getLatentHeat();
    double tau = input->getMagmaProperties()["constantRelaxationCrystallization"]["tau"].get<double>();

    const auto& shear_heat = dike->shear_heat;

    double E0 = 0.5 * schedule->getMassRate(t1, t1 + dt) * schedule->getMagmaChamberTemperature() * Cm;
    int N = ny+nr+1;
    for (int ix = 0; ix < nx; ix++){
        if (h[ix] < MIN_MOBILITY_WIDTH) continue;
        std::vector<double> a(N, 0.0), b(N, 0.0), c(N, 0.0), rhs(N, 0.0);
        const ArrayXd Told = old_dike->temperature.row(ix);
        const ArrayXd rho = dike->density.row(ix);
        const ArrayXd rho_old = old_dike->density.row(ix);
        ArrayXd Ein = ArrayXd::Zero(ny);
        if (ix == 0){
            Ein.fill(E0*dy);
        }
        else{
            Ein = rho * Cm * dike->temperature.row(ix-1).transpose() * qx.row(ix).transpose();
        }
        double dyt, dyb, vtp, vtm, vbp, vbm;
        double Vnew = h[ix]*dx*dy;
        double Vold = hold[ix]*dx*dy;
        dyt = yc[1] - yc[0];
        vtp = std::max(qy(ix, 1), 0.0); // qy in top
        vtm = -std::max(-qy(ix, 1), 0.0); // 0 if qy > 0 in top
        a[0] = 0;
        b[0] = rho[0]*Cm*(Vnew/dt + qx(ix+1, 0) + vtp) + dx*km/h[ix]*(1.0/dyt);
        c[0] = rho[1]*Cm*vtm - dx*km/h[ix]*(1.0/dyt);
        rhs[0] = rho_old[0]*Cm*Told[0]*Vold/dt + Ein[0];
        rhs[0] += latentHeatCrystallization*(1.0 - alpha(ix, 0))*rhoc0*Lm*Vnew*(betaeq(ix, 0) - beta(ix, 0))/tau
                + shearHeating * shear_heat(ix, 0);

        for (int iy = 1; iy < ny-1; iy++){
            dyt = yc[iy+1] - yc[iy];
            dyb = yc[iy] - yc[iy-1];
            vtp = std::max(qy(ix, iy+1), 0.0); // qy in top
            vtm = -std::max(-qy(ix, iy+1), 0.0); // 0 if qy > 0 in top
            vbp = std::max(qy(ix, iy), 0.0); // qy in bot
            vbm = -std::max(-qy(ix, iy), 0.0); // 0 if qy > 0 in bot
            a[iy] = -rho[iy-1]*Cm*vbp - dx*km/h[ix]/dyb;
            b[iy] = rho[iy]*Cm*(Vnew/dt + qx(ix+1, iy) + vtp - vbm) + dx*km/h[ix]*(1.0/dyt + 1.0/dyb);
            c[iy] = rho[iy+1]*Cm*vtm - dx*km/h[ix]/dyt;
            rhs[iy] = rho_old[iy]*Cm*Told[iy]*Vold/dt + Ein[iy];
            rhs[iy] += latentHeatCrystallization*(1.0 - alpha(ix, iy))*rhoc0*Lm*Vnew*(betaeq(ix, iy) - beta(ix, iy))/tau + 
                     + shearHeating * shear_heat(ix, iy);
        }
        dyt = yb[ny] - yc[ny-1];
        dyb = yc[ny-1] - yc[ny-2];
        vbp = std::max(qy(ix, ny-1), 0.0); // qy in bot
        vbm = -std::max(-qy(ix, ny-1), 0.0); // 0 if qy > 0 in bot
        a[ny-1] = -rho[ny-2]*Cm*vbp - dx*km/h[ix]/dyb;
        b[ny-1] = rho[ny-1]*Cm*(Vnew/dt + qx(ix+1, ny-1) - vbm) + dx*km/h[ix]*(1.0/dyt + 1.0/dyb);
        c[ny-1] = -dx*km/h[ix]/dyt;
        rhs[ny-1] = rho_old[ny-1]*Cm*Told[ny-1]*Vold/dt + Ein[ny-1];
        rhs[ny-1] += latentHeatCrystallization*(1.0 - alpha(ix, ny-1))*rhoc0*Lm*Vnew*(betaeq(ix, ny-1) - beta(ix, ny-1))/tau + 
                   + shearHeating * shear_heat(ix, ny-1);
        a[ny] = -dx*km/h[ix]/dyt;
        b[ny] = dx*km/h[ix]/dyt;

        /* Reservoir conduction */
        int ind = ny+1;
        double rhor = reservoir->density[ix];
        auto Tr = reservoir->temperature.row(ix);
        double ktop, kbot;
        dyt = ycr[1] - ycr[0];
        dyb = ycr[0] - ybr[0];
        double heat_coef = -kr / dyb;
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
        dike->magma_to_rock_heat_flux[ix] = heat_coef * (sol[ny+1] - sol[ny]);
    }
    JUST_TIMER_STOP("2) Energy balance time");
}


void DikeModel2d::updateCrystallization(){
    JUST_TIMER_START("2) Update crystal time");
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
    double tau = input->getMagmaProperties()["constantRelaxationCrystallization"]["tau"].get<double>();

    // double Q0 = 0.5 * schedule->getMassRate(t1, t1 + dt) / schedule->getMagmaChamberDensity() * schedule->getMagmaChamberCrystallization();
    double Q0 = 0.5 * schedule->getMassRate(t1, t1 + dt) / schedule->getMagmaChamberDensity() * magma_state->getMagmaChamberCrystallization();
    for (int ix = 0; ix < nx; ix++){
        if (h[ix] < MIN_MOBILITY_WIDTH) continue;
        std::vector<double> a(ny, 0.0), b(ny, 0.0), c(ny, 0.0), rhs(ny, 0.0);
        const ArrayXd bold = old_dike->beta.row(ix);
        const ArrayXd beq = dike->betaeq.row(ix);
        ArrayXd Ein = ArrayXd::Zero(ny);
        const ArrayXd alpha = dike->alpha.row(ix);
        const ArrayXd alpha_old = old_dike->alpha.row(ix);
        if (ix == 0){
            Ein.fill(Q0*dy);
        }
        else{
            Ein = (1.0 - dike->alpha.row(ix-1)) * dike->beta.row(ix-1).array() * qx.row(ix).array();
        }
        double dyt, dyb, vtp, vtm, vbp, vbm;
        dyt = yc[1] - yc[0];
        vtp = std::max(qy(ix, 1), 0.0); // qy in top
        vtm = -std::max(-qy(ix, 1), 0.0); // 0 if qy > 0 in top
        a[0] = 0;
        b[0] = (1.0 - alpha[0])*(h[ix]*dx*dy/dt + qx(ix+1, 0) + vtp) + (1.0 - alpha[0])*h[ix]*dx*dy/tau;
        c[0] = (1.0 - alpha[1])*vtm;
        rhs[0] = (1.0 - alpha_old[0])*bold[0]*hold[ix]*dx*dy/dt + Ein[0] + (1.0 - alpha[0])*beq[0]*h[ix]*dx*dy/tau;

        for (int iy = 1; iy < ny-1; iy++){
            vtp = std::max(qy(ix, iy+1), 0.0); // qy in top
            vtm = -std::max(-qy(ix, iy+1), 0.0); // 0 if qy > 0 in top
            vbp = std::max(qy(ix, iy), 0.0); // qy in bot
            vbm = -std::max(-qy(ix, iy), 0.0); // 0 if qy > 0 in bot
            a[iy] = -(1.0 - alpha[iy-1])*vbp;
            b[iy] = (1.0 - alpha[iy])*(h[ix]*dx*dy/dt + qx(ix+1, iy) + vtp - vbm) + (1.0 - alpha[iy])*h[ix]*dx*dy/tau;
            c[iy] = (1.0 - alpha[iy+1])*vtm;
            rhs[iy] = (1.0 - alpha_old[iy])*bold[iy]*hold[ix]*dx*dy/dt + Ein[iy] + (1.0 - alpha[iy])*beq[iy]*h[ix]*dx*dy/tau;
        }
        vbp = std::max(qy(ix, ny-1), 0.0); // qy in bot
        vbm = -std::max(-qy(ix, ny-1), 0.0); // 0 if qy > 0 in bot
        a[ny-1] = -(1.0 - alpha[ny-2])*vbp;
        b[ny-1] = (1.0 - alpha[ny-1])*(h[ix]*dx*dy/dt + qx(ix+1, ny-1) - vbm) + (1.0 - alpha[ny-1])*h[ix]*dx*dy/tau;
        c[ny-1] = 0;
        rhs[ny-1] = (1.0 - alpha_old[ny-1])*bold[ny-1]*hold[ix]*dx*dy/dt + Ein[ny-1] + (1.0 - alpha[ny-1])*beq[ny-1]*h[ix]*dx*dy/tau;

        /* Solve tridiagonal system */
        auto sol = Utils::tridiagonal_solver(a, b, c, rhs);
        dike->beta.row(ix) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sol.data(), ny);
    }
    JUST_TIMER_STOP("2) Update crystal time");
}


void DikeModel2d::implicitMassBalance(){
    using Eigen::seq;
    using Eigen::last;
    int nx = mesh->size();
    int ny = dike->getLayersNumber();
    auto& hold = old_dike->hw;
    VectorXd rho_avg = dike->density.rowwise().mean();
    VectorXd rho_avg_old = old_dike->density.rowwise().mean();
    auto& rho = dike->density;
    auto& rho_old = old_dike->density;
    auto& viscosity = dike->viscosity; 
    double dx = mesh->getdx();
    double dt = timestep_controller->getCurrentTimestep();
    double t1 = timestep_controller->getCurrentTime();
    const auto& yb = dike->yb;
    const auto& yc = dike->yc;
    ArrayXd ybt = yb(seq(1, ny));
    ArrayXd ybb = yb(seq(0, ny-1));
    double dy = yb[1] - yb[0];
    double g = reservoir->getGravityAcceleration();
    
    int tip_old = old_dike->tip_element;
    int N = std::min({tip_old + 2, nx - 1}) + 1;
    ArrayXd W = dike->hw(seq(0, N-1));
    ArrayXd Witer = dike->hw(seq(0, N-1));
    ArrayXd P = dike->pressure(seq(0, N-1));
    VectorXd sol = VectorXd::Zero(2*N); // (W, P)
    VectorXd rhs_base = VectorXd::Zero(2*N);
    MatrixXd mat = MatrixXd::Zero(2*N, 2*N);
    MatrixXd Iww = -2*elasticity->getMatrix()(seq(0, N-1), seq(0, N-1));
    MatrixXd IwwA = -2*elasticity->getA()(seq(0, N-1), seq(0, N-1));
    MatrixXd IwwB = -2*elasticity->getB()(seq(0, N-1), seq(0, N-1));
    VectorXd sigmah = reservoir->getLithostaticPressure()(seq(0, N-1));
    MatrixXd Iwp = MatrixXd::Identity(N, N);
    VectorXd diag_pw = rho_avg(seq(0, N-1)) * dx / dt;
    MatrixXd Ipw = diag_pw.asDiagonal();
    rhs_base(seq(0, N-1)) = sigmah;
    rhs_base(seq(N, last)) = (rho_avg_old(seq(0, N-1)) * dx / dt).array() * hold(seq(0, N-1)).array();
    double M0 = 0.5 * schedule->getMassRate(t1, t1 + dt);
    rhs_base(N) = rhs_base(N) + M0;

    ArrayXd G = ArrayXd::Zero(N+1);
    ArrayXd lambda = ArrayXd::Zero(N+1);
    ArrayXd rhog = ArrayXd::Zero(N+1);
    ArrayXd mobility = ArrayXd::Zero(N+1);
    ArrayXd Qx = ArrayXd::Zero(N+1);
    ArrayXd lambda_total = ArrayXd::Zero(N+1);
    ArrayXd Mx = ArrayXd::Zero(N+1);
    for (int i = 1; i < N; ++i){
        ArrayXd mul = viscosity.row(i-1);
        ArrayXd mur = viscosity.row(i);
        ArrayXd mu = mul;
        if (VISCOSITY_APPROXIMATION == "harmonic"){
            mu = 2 * (mul * mur) / (mul + mur);
        }
        if (VISCOSITY_APPROXIMATION == "mean"){
            mu = 0.5 * (mul + mur);
        }
        ArrayXd acoef = 0.5 * mu.cwiseInverse();
        ArrayXd ccoef = ArrayXd::Zero(ny);
        ccoef(ny-1) = -acoef(ny-1);
        for (int iy = ny-2; iy >= 0; iy--){
            ccoef(iy) = ccoef(iy+1) + yb(iy+1) * yb(iy+1) * (acoef(iy+1) - acoef(iy));
        }
        ArrayXd lambda = acoef * (ybt.cube() - ybb.cube()) / 3 + ccoef * (ybt - ybb);
        lambda_total(i) = (lambda * rho.row(i-1).array().transpose()).sum();
        rhog(i) = 0.5 * (rho_avg(i) + rho_avg(i-1)) * g;
    }
    
    JUST_TIMER_START("2) Solve nonlinear equations");
    int stab_iter = 0;
    double error = 0;
    int iter = 0;
    for (iter = 0; iter < MAX_ITERATIONS; ++iter){
        JUST_TIMER_START("3) Assemble matrix");
        VectorXd rhs = rhs_base;
        VectorXd rhsp = VectorXd::Zero(N);
        MatrixXd Ipp = MatrixXd::Zero(N, N);
        VectorXd sol_iter(2*N);
        sol_iter << W, P;
        Witer = W;
        for (int i = 1; i < N; ++i){
            double h = 0.5 * (W[i] + W[i-1]);
            if (h < MIN_MOBILITY_WIDTH) h = 0;
            mobility(i) = h*h*h*lambda_total(i);
            G(i) = (P(i) - P(i-1)) / dx + rhog(i);
            if (G(i) > 0) mobility(i) = 0;
            rhsp(i-1) += -mobility(i) * rhog(i);
            rhsp(i) += mobility(i) * rhog(i);
            double coef = mobility(i) / dx;
            Ipp(i-1, i-1) += -coef;
            Ipp(i-1, i) += coef;
            Ipp(i, i-1) += coef;
            Ipp(i, i) += -coef;
        }
        rhs(seq(N, last)) += rhsp;
        if (algorithm_properties["isSparseElasticity"].get<bool>() == false){
            mat << Iww, Iwp, Ipw, Ipp;
        }
        else{
            VectorXd rhsw = IwwB * W.matrix();
            rhs(seq(0, N-1)) -= rhsw;
            mat << IwwA, Iwp, Ipw, Ipp;
        }

        if (algorithm_properties["isCohesiveStress"].get<bool>()){
            VectorXd rhs_coh = VectorXd::Zero(N);
            if (N > cohesive_model->getCohesiveElementsSize()){
                for (int i = N - cohesive_model->getCohesiveElementsSize(); i < N; i++){
                    rhs_coh[i] = cohesive_model->getCohesiveStress(W[i]);
                }
            }
            rhs(seq(0, N-1)) += rhs_coh;
        }
        JUST_TIMER_STOP("3) Assemble matrix");

        JUST_TIMER_START("3) Solver time");
        if (algorithm_properties["solverName"] == "denselu"){
            Eigen::FullPivLU<MatrixXd> lu(mat);
            sol = lu.solve(rhs);
        }
        else if(algorithm_properties["solverName"] == "umfpack"){
            Eigen::SparseMatrix<double> mat_sparse = mat.sparseView();
            LinearSolver::solveUmfpack(mat_sparse, rhs, sol);
        }
        else if(algorithm_properties["solverName"] == "pardiso"){
            Eigen::SparseMatrix<double> mat_sparse = mat.sparseView();
            LinearSolver::solvePardiso(mat_sparse, rhs, sol);
        }
        JUST_TIMER_STOP("3) Solver time");
        /* Check for Nan value */
        for (int i = 0; i < sol.size(); ++i){
            if (std::isnan(sol[i])){
                solver_log.successful = false;
                return;
            }
        }
        // error = (sol-sol_iter).lpNorm<Eigen::Infinity>() / std::max(sol.lpNorm<Eigen::Infinity>(), 1e-4);
        error = 0;
        for (int iw = 0; iw < N; iw++){
            double erri = std::abs(sol(iw) - sol_iter(iw)) / std::max(1e-4, std::abs(sol(iw)));
            error = error < erri ? erri : error;
        }
        W = sol(seq(0, N-1));
        P = sol(seq(N, last));
        if (error < TOLERANCE){
            stab_iter++;
        }
        else{
            stab_iter = 0;
        }
        if (stab_iter == 2) break;
    }
    spdlog::trace("-- {:<10}) err = {:<8}, niter. = {:<3}", timestep_controller->getTimeIteration(), error, iter);
    JUST_TIMER_STOP("2) Solve nonlinear equations");
    for (int i = 1; i < N; ++i){
        ArrayXd mul = viscosity.row(i-1);
        ArrayXd mur = viscosity.row(i);
        ArrayXd mu = mul;
        double h = 0.5 * (Witer[i] + Witer[i-1]);
        G(i) = (P(i) - P(i-1)) / dx + rhog(i);
        if (G(i) > 0) G(i) = 0;
        dike->G(i) = G(i);
        if (VISCOSITY_APPROXIMATION == "harmonic"){
            mu = 2 * (mul * mur) / (mul + mur);
        }
        if (VISCOSITY_APPROXIMATION == "mean"){
            mu = 0.5 * (mul + mur);
        }
        ArrayXd acoef = 0.5 * mu.cwiseInverse() * h * h * G(i);
        ArrayXd ccoef = ArrayXd::Zero(ny);
        ccoef(ny-1) = -acoef(ny-1);
        for (int iy = ny-2; iy >= 0; iy--){
            ccoef(iy) = ccoef(iy+1) + yb(iy+1) * yb(iy+1) * (acoef(iy+1) - acoef(iy));
        }
        dike->qx.row(i) = h*(acoef * (ybt.cube() - ybb.cube()) / 3 + ccoef * (ybt - ybb));
        dike->A.row(i) = acoef;
        dike->C.row(i) = ccoef;
        dike->Qx(i) = dike->qx.row(i).sum();
        dike->Mx(i) = mobility(i) * G(i);
    }

    auto &qy = dike->qy;
    const auto &qx = dike->qx;
    ArrayXd Min = ArrayXd::Zero(N);
    Min(0) = M0;
    for (int ix = 0; ix < N; ++ix){
        dike->hw(ix) = W(ix);
        dike->pressure(ix) = P(ix);
        double Gcell;
        if (ix == 0 || ix >= dike->tip_element - 1){
            Gcell = 0.0;
        }
        else{
            Gcell = 0.5 * (G(ix+1) + G(ix));
        }
        Gcell = std::abs(Gcell) < 1000*g ? Gcell : -1000*g;
        dike->shear_heat.row(ix) = dx * (std::pow(W(ix), 3) * Gcell*Gcell / 3) * (ybt.cube() - ybb.cube()).array() / dike->viscosity.row(ix).array().transpose();
        
        ArrayXd rho_jb = 0.5 * (rho.row(ix)(seq(1, ny-1)) + rho.row(ix)(seq(0, ny-2)));
        ArrayXd rhoj = rho.row(ix);
        rhoj(seq(1, ny-1)) = rho_jb;
        for (int iy = 1; iy < ny; iy++){
            qy(ix, iy) = (
                rhoj(iy-1) * qy(ix, iy-1) - 
                rho(ix, iy-1) * qx(ix+1, iy-1) + 
                rho(std::max(0, ix-1), iy-1) * qx(ix, iy-1) - 
                dy*dx/dt*(
                    rho(ix, iy-1)*W[ix] - rho_old(ix, iy-1)*hold[ix]
                ) + 
                Min[ix]*dy
            ) / rhoj(iy);
        }
    }
    dike->overpressure = 2 * (elasticity->getMatrix() * dike->hw.matrix());
    dike->updateOpenElements();
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
    dump(file, "rhoc", dike->rhoc);
    dump(file, "rhog", dike->rhog);
    dump(file, "rhom", dike->rhom);

    dump(file, "temperature", dike->temperature);
    dump(file, "viscosity", dike->viscosity);
    dump(file, "Twall", dike->Twall);
    dump(file, "qx", dike->qx);
    dump(file, "qy", dike->qy);
    dump(file, "A", dike->A);
    dump(file, "C", dike->C);
    dump(file, "G", dike->G);
    dump(file, "alpha", dike->alpha);
    dump(file, "beta", dike->beta);
    dump(file, "gamma", dike->gamma);
    dump(file, "betaeq", dike->betaeq);
    dump(file, "Tliquidus", dike->Tliquidus);
    dump(file, "Tsolidus", dike->Tsolidus);
    dump(file, "TotalFluxElements", dike->Qx);
    dump(file, "TotalMassRateElements", dike->Mx);
    dump(file, "TipElement", dike->tip_element);
    dump(file, "TipFront", dike->front);
    dump(file, "MagmaToRockHeatFlux", dike->magma_to_rock_heat_flux);
    dump(file, "ShearHeat", dike->shear_heat);

    dump(file, "reservoir/yc", reservoir->yc);
    dump(file, "reservoir/yb", reservoir->yb);
    dump(file, "reservoir/temperature", reservoir->temperature);
    dump(file, "reservoir/ny", reservoir->ny);
}


void DikeModel2d::updateData(){
    double dt = timestep_controller->getCurrentTimestep();
    timestep_controller->update();
    dike->time = timestep_controller->getCurrentTime();
    reservoir->time = timestep_controller->getCurrentTime();
    *old_dike = *dike;
    auto [is_save, save_timestep] = timestep_controller->saveTimestepIteration();
    if (is_save){
        JUST_TIMER_START("2) Data saving time");
        spdlog::info("{:>7} | {:>8.4} h | dt = {:>8.4} | tot. iter. = {:>8} |",
            save_timestep,
            timestep_controller->getCurrentTime() / 3600.0,
            dt,
            timestep_controller->getAllIterations()
        );
        std::string savepath = (input->getDataDir() / "data_").string() + std::to_string(save_timestep) + ".h5";
        saveData(savepath);

        std::string reservoir_path = (input->getReservoirDir() / "data_").string() + std::to_string(save_timestep) + ".h5";
        reservoir->saveData(reservoir_path);
        JUST_TIMER_STOP("2) Data saving time");
    }
    // if (save_timestep > 362){
    //     // std::cout << "dt = " << dt << std::endl;
    //     if (dt < 1e-5) throw std::invalid_argument("fuck this all!\n");
    // }
}


void DikeModel2d::reloadData(){
    *dike = *old_dike;
    timestep_controller->divideTimestep(solver_log.cfl_ratio);
}
