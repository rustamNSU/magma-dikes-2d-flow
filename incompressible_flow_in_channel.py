import numpy as np
from typing import Sequence
from copy import deepcopy
import matplotlib.pyplot as plt


from src.utils.base_magma_state import mu_melt, beta_T, theta
from src.core.fixed_channel import FixedChannel

TL = 1250 + 273.15 # [K]
TS = 1000 + 273.15 # [K]
alpha = 5          # [J/(m^2*K)]
rho = 2300         # [kg/m^3]
Cp = 1200          # [J/(kg*K)]
k = 2              # [J/()]

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
dt = 100.0
Nx = 100
Ny = 10

tolerance = 1e-4
max_iter = 30

Tr = lambda x: (X2 - x) / (X2 - X1) * T1 + (x - X1) / (X2 - X1) * T2
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
    Tr=np.minimum(T0, Tr(fc_old.xc)),
    alpha=alpha,
    k=k
)

fc_old.update_beta(beta_T1)
fc_old.update_viscosity(mu_magma)


t_cur = t_start + dt
while t_cur < t_end:
    fc_new = deepcopy(fc_old)
    fc_new.time = t_cur
    err_T = 1.0
    for iter in range(max_iter):
        # Solve lubrication
        for ix in range(Nx):
            fc_new.define_mobility(ix)
            fc_new.define_dp(ix, Q / 2)
            
        # Solve energy conservation
        iterT2d = deepcopy(fc_new.T2d)
        fc_new.solve_convect_heat_equation(fc_old)
        fc_new.solve_diffusion_heat_equation(fc_old)
        err_T = np.linalg.norm(iterT2d - fc_new.T2d) / np.linalg.norm(iterT2d)
        if err_T < 1e-6:
            break
    print("iter = {}, error = {}".format(iter, err_T))
    fc_old = fc_new



a = 10



# T = np.linspace(600, 1400)
# plt.semilogy(T-273.15, mu_magma(T), linewidth=3)
# plt.grid()
# plt.show()

# zz = xc2d/10000 + yc2d
# plt.contourf(xc2d, yc2d, zz)
# plt.show()
