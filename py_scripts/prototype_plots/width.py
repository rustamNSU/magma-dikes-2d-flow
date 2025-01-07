import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_scripts.utils import set_matplotlib_settings, create_layers_mask
from py_scripts.DikeData import DikeData
set_matplotlib_settings(DEFAULT_SIZE=8, LEGEND_SIZE=8, MAC_OS=True)


simIDs = [2, 70, 71]
# simIDs = [20, 25]
simLegends = [str(simID) for simID in simIDs]
wlim = (0, 3)
plim = (0, 900)
simPaths = [sim_dir + f"/simID{simID}" for simID in simIDs]
dikes = [DikeData(sim_path, step_rate=10) for sim_path in simPaths]
xlim = (
    min(data.xmin for data in dikes),
    max(data.xmax for data in dikes)
)


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
ax1.set_xlabel(r"Depth (m)")
ax1.set_ylabel(r"Width (m)")
ax1.set_xlim(xlim)
ax1.set_ylim(wlim)
ax1.grid()

ax2.set_xlabel(r"Depth (m)")
ax2.set_ylabel(r"Pressure (MPa)")
ax2.set_xlim(xlim)
ax2.set_ylim(plim)
ax2.grid()


laxWs   = []
laxPs   = []
colors = ['k', 'b', 'r', 'g']
for i in range(len(dikes)):
    data = dikes[i].data
    x = data[0]["xc"]
    xb = data[0]["xb"]
    w = data[0]["width"]
    p = 1e-6*data[0]["pressure"]
    laxW,   = ax1.plot(x, w, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxP,   = ax2.plot(x, p, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxWs.append(laxW) 
    laxPs.append(laxP) 

ax1.legend().set_draggable(True)
ax2.legend().set_draggable(True)

ax_time_slider = fig.add_axes([0.25, 0.05, 0.65, 0.03])
time_slider = Slider(
    ax=ax_time_slider, 
    label='timesteps', 
    valmin=0, 
    valmax=min(len(dike.data)-1 for dike in dikes),
    valinit=0, 
    valstep=1)


def time_update(val):
    t = time_slider.val
    fig.suptitle('time = {0:.2f} h'.format(dikes[0].data[t]["time"]/3600))
    for i in range(len(dikes)):
        data = dikes[i].data[t]
        w = data["width"]
        p = 1e-6*data["pressure"]
        
        laxWs[i].set_ydata(w) 
        laxPs[i].set_ydata(p) 

time_slider.on_changed(time_update)
plt.show()