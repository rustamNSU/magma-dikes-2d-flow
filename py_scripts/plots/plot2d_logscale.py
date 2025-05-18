import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.widgets import Slider
from py_scripts.utils import set_matplotlib_settings, create_layers_mask
from pysrc import *

set_matplotlib_settings()
simID = 100
sim_path = sim_dir + f"/simID{simID}"
dike = DikeData(sim_path, step_rate=10)
fig, ax = plt.subplots(figsize=(8, 3))
frame = dike.data[0]

# Prepare data: replace values <= 1.0 with 1.0 and apply open_mask
Z = np.where(frame.tau <= 1.0, 1.0, frame.tau)
Z = np.where(frame.open_mask[:, None], Z, np.nan)

# Compute color limits ignoring NaNs
vmin = np.nanmin(Z)
vmax = np.nanmax(Z)
if np.isnan(vmin) or np.isnan(vmax):
    vmin, vmax = 1.0, 10.0

pcm = ax.pcolormesh(frame.xb, frame.yb, Z.T, shading='flat', cmap='jet', norm=LogNorm(vmin=vmin, vmax=vmax))
cb = fig.colorbar(pcm, ax=ax)

ax_time_slider = fig.add_axes([0.15, 0.03, 0.6, 0.015])
time_slider = Slider(
    ax=ax_time_slider, 
    label='timesteps', 
    valmin=0, 
    valmax=len(dike.data)-1,
    valinit=0, 
    valstep=1
)

def time_update(val):
    t = int(time_slider.val)
    frame = dike.data[t]
    Z = np.where(frame.tau <= 1.0, 1.0, frame.tau)
    Z = np.where(frame.open_mask[:, None], Z, np.nan)
    vmin = np.nanmin(Z)
    vmax = np.nanmax(Z)
    if np.isnan(vmin) or np.isnan(vmax):
        vmin, vmax = 1.0, 10.0

    pcm.set_array(Z.T.ravel())
    pcm.norm = LogNorm(vmin=vmin, vmax=vmax)
    cb.update_normal(pcm)
    fig.canvas.draw_idle()

time_slider.on_changed(time_update)
plt.show()