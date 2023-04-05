import numpy as np
from typing import Sequence
from copy import deepcopy
import pickle


from src.utils.base_magma_state import mu_melt, beta_T, theta
from src.core.fixed_channel import FixedChannel
from src.utils.adaptive_timestep import AdaptiveTimestep

simID = 21050
TL = 1250 + 273.15 # [K]
TS = 1000 + 273.15 # [K]
alpha = 50          # [J/(m^2*K)]
rho = 2300         # [kg/m^3]
Cp = 1200          # [J/(kg*K)]
k = 20              # [J/()]

T1 = 1300.0 + 273.15
T2 = 800.0 + 273.15
T0 = TL + 30

X1 = 0.0
X2 = 10000.0

W = 1.0
Q = 2.0
P2 = 0.0

t_start = 0.0
t_end = 10000.0
t_cur = 0.0
basic_dt = 25.0
Nx = 50
Ny = 10

tolerance = 1e-6
max_iter = 20

Tr = lambda x: T0 if (x < 3000.0 and x > 6000) else T2
mu_magma = lambda T: mu_melt(T) * theta(beta_T(T, TL, TS))
beta_T1 = lambda T: beta_T(T, TL, TS)

fc_old = FixedChannel(
    xmin=X1,
    xmax=X2,
    width=W,
    nx=Nx,
    ny=Ny)

fc_old.set_initial_state(
    T0=T0,
    Cp=Cp,
    rho=rho,
    mu0=mu_magma(T0),
    beta0=beta_T(T0, TL, TS),
    t0=t_start,
    Tr=np.array([Tr(x) for x in fc_old.xc]),
    alpha=alpha,
    k=k
)

fc_old.update_beta(beta_T1)
fc_old.update_viscosity(mu_magma)

timestep_controller = AdaptiveTimestep(basic_dt)
t_cur = t_start
dt = timestep_controller.get_timestep()
while t_cur + dt < t_end:
    if abs(dt - basic_dt) < 1e-9:
        print("{:10.4f} -> {:10.4f}:".format(t_cur, t_cur + dt))
    print("-- dt = {}".format(dt))
    fc_new = deepcopy(fc_old)
    fc_new.time = t_cur + dt
    err_T = 1.0
    iter = 1
    while iter < max_iter:
        # Solve lubrication
        for ix in range(Nx):
            fc_new.define_mobility(ix)
            fc_new.define_dp(ix, Q / 2)
            
        # Solve energy conservation
        iterT2d = deepcopy(fc_new.T2d)
        fc_new.solve_convect_heat_equation(fc_old)
        fc_new.solve_diffusion_heat_equation(fc_old)
        
        fc_new.update_beta(beta_T1)
        fc_new.update_viscosity(mu_magma)
        err_T = np.linalg.norm(iterT2d - fc_new.T2d) / np.linalg.norm(iterT2d)
        if err_T < tolerance:
            break
        iter += 1
    
    if iter >= max_iter:
        timestep_controller.divide_timestep(dt)
        dt = timestep_controller.get_timestep()
        print("  divide timestep")
        continue
    print("  iter = {}, err = {}".format(iter, err_T))

    fc_old = fc_new
    t_cur += dt
    dt = timestep_controller.get_timestep()
    

filename = "simulations" + "/simID{}.pkl".format(simID)
with open(filename, 'wb') as outp:
    pickle.dump(fc_new, outp, pickle.HIGHEST_PROTOCOL)


a = 10



# T = np.linspace(600, 1400)
# plt.semilogy(T-273.15, mu_magma(T), linewidth=3)
# plt.grid()
# plt.show()

# zz = xc2d/10000 + yc2d
# plt.contourf(xc2d, yc2d, zz)
# plt.show()
