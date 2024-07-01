import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_scripts.utils import set_matplotlib_settings
from py_scripts.DikeData import DikeData
set_matplotlib_settings()


simID = 4
sim_path = sim_dir + f"/simID{simID}"
wlim = (0, 3)
plim = (0, 700)
Tmin = 600
Tmax = 920
xlim = (-30000, -10000.0)


dike = DikeData(sim_path, step_rate=10)
data = dike.data
fig, ax = plt.subplots(1, 1, figsize=(16, 4))
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
fig.subplots_adjust(bottom=0.2)
ax.set_xlabel(r"Depth")
ax.set_ylabel(r"$\xi = y/h(x,t)$")

x = data[0]["xc"]
y = data[0]["yb"]
ux = data[0]["ux"][1:-1, :].T
# pax = ax.pcolormesh(x, y, qx, vmin=0, vmax=0.1, shading='flat')
pax = ax.pcolormesh(x, y, ux, shading='flat')
cb = fig.colorbar(pax, cax=cax)

ax_time_slider = fig.add_axes([0.2, 0.05, 0.6, 0.05])
time_slider = Slider(
    ax=ax_time_slider, 
    label='timesteps', 
    valmin=0, 
    valmax=len(data)-1,
    valinit=0, 
    valstep=1)


def time_update(val):
    t = time_slider.val
    fig.suptitle('time = {0:.2f} h'.format(data[t]["time"]/3600))
    ux = data[t]["ux"][1:-1, :].T
    # pax.set_array(ux.ravel())
    ax.clear()
    pax = ax.pcolormesh(x, y, ux, shading='flat')
    fig.colorbar(pax, cax=cax)

time_slider.on_changed(time_update)
plt.show()