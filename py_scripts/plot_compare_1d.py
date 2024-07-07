import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from py_scripts.utils import set_matplotlib_settings
from py_scripts.DikeData import DikeData
set_matplotlib_settings()

# simLegends = [str(simID) for simID in simIDs]
simIDs = [14, 13, 12]
simLegends = [r"$L^* = 350000$ + shear heating", r"$L^* = 350000$", r"$L^* = 0$"]
wlim = (0, 3)
plim = (0, 700)
Tlim = (600, 920)
xlim = (-30000, 0)
# simIDs = [1, 3, 2]
# simLegends = [r"$\Delta x = 50$m", r"$\Delta x = 100$m", r"$\Delta x = 200$m"]
# simIDs = [5, 6]
# simLegends = [r"$\mu_{i+1/2} = \frac{2}{1/\mu_i + 1/\mu_{i+1}}$", r"$\mu_{i+1/2} = \mu_i$"]
# wlim = (0, 3)
# plim = (0, 700)
# Tlim = (600, 920)
# xlim = (-30000, -10000.0)

simPaths = [sim_dir + f"/simID{simID}" for simID in simIDs]


dikes = [DikeData(sim_path, step_rate=1) for sim_path in simPaths]
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12))
fig.subplots_adjust(bottom=0.2)
ax1.set_xlabel(r"Depth (m)")
ax1.set_ylabel(r"Width (m)")
ax2.set_xlabel(r"Depth (m)")
ax2.set_ylabel(r"Pressure (MPa)")
ax3.set_xlabel(r"Depth (m)")
ax3.set_ylabel(r"Temperature C$^\circ$")

ax1.set_xlim(xlim)
ax2.set_xlim(xlim)
ax3.set_xlim(xlim)

ax1.set_ylim(wlim)
ax2.set_ylim(plim)
ax3.set_ylim(Tlim)

ax1.grid()
ax2.grid()
ax3.grid()


lws = []
lps = []
lTs = []
colors = ['k', 'b', 'r', 'g']
for i in range(len(dikes)):
    data = dikes[i].data
    x = data[0]["xc"]
    w = data[0]["width"]
    p = 1e-6*data[0]["pressure"]
    T = data[0]["temperature"][:, 0]
    lw, = ax1.plot(x, w, lw=2, ls='-', color=colors[i], label=simLegends[i])
    lp, = ax2.plot(x, p, lw=2, ls='-', color=colors[i], label=simLegends[i])
    lT, = ax3.plot(x, T, lw=2, ls='-', color=colors[i], label=simLegends[i])
    lws.append(lw)
    lps.append(lp)
    lTs.append(lT)

ax1.legend().set_draggable(True)
ax2.legend().set_draggable(True)
ax3.legend().set_draggable(True)

ax_time_slider = fig.add_axes([0.2, 0.05, 0.6, 0.05])
time_slider = Slider(
    ax=ax_time_slider, 
    label='timesteps', 
    valmin=0, 
    valmax=len(dikes[0].data)-1,
    valinit=0, 
    valstep=1)


def time_update(val):
    t = time_slider.val
    fig.suptitle('time = {0:.2f} h'.format(dikes[0].data[t]["time"]/3600))
    for i in range(len(dikes)):
        data = dikes[i].data[t]
        w = data["width"]
        p = 1e-6*data["pressure"]
        T = data["temperature"][:, 0]
        lws[i].set_ydata(w)
        lps[i].set_ydata(p)
        lTs[i].set_ydata(T)

time_slider.on_changed(time_update)
plt.show()