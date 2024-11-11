import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)


import json
import numpy as np
import matplotlib.pyplot as plt
from py_scripts.utils import set_matplotlib_settings, create_layers_mask
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)


sim_id = 1
datapath = sim_dir + f"/simID{sim_id}/test_magma_properties.json"
data = json.load(open(datapath))


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 3))
P = 1e-6 * np.array(data["dissolved"]["pressure"])
ax1.plot(P, data["dissolved"]["wth2o"], lw=2)
ax2.plot(P, data["dissolved"]["wtco2"], lw=2)
ax3.plot(P, data["dissolved"]["xh2og"], lw=2)
for ax, ylabel in zip([ax1, ax2, ax3], [r"H$_2$O, wt. %", r"CO$_2$, wt. %", r"x$_{H_2O, g}$"]):
    ax.set_xlabel("Pressure, MPa")
    ax.set_ylabel(ylabel)
    ax.grid()
fig.tight_layout()

fig = plt.figure(figsize=(12, 4))
ax1 = fig.add_subplot(1, 3, 1, projection='3d')
ax2 = fig.add_subplot(1, 3, 2, projection='3d')
ax3 = fig.add_subplot(1, 3, 3, projection='3d')

pressure = 1e-6 * np.array(data["exsolved"]["pressure"])
temperature = np.array(data["exsolved"]["temperature"])
P, T = np.meshgrid(pressure, temperature)
xg = []
densities = []
for case in ["case1", "case2", "case3"]:
    xg.append(data["exsolved"][case]["xh2og"])
    density = np.array(data["exsolved"][case]["density"]).reshape((len(pressure), len(temperature))).transpose()
    densities.append(density)

for i, ax in zip(range(3), [ax1, ax2, ax3]):
    D = densities[i]
    ax.plot_surface(P, T, D, cmap=plt.cm.YlGnBu_r)
    ax.set_xlabel(r"Pressure, MPa")
    ax.set_ylabel(r"Temperature, C$^\circ$")
    ax.set_zlabel(r"Gas density, kg/m$^3$")
    ax.set_title(r"x$_{{H_2O, g}}={}$".format(xg[i]))

plt.show()