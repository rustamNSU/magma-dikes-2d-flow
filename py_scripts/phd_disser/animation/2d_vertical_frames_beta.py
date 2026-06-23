import sys, os
repository_dir = os.path.abspath(os.getcwd())
sys.path.append(repository_dir)


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_scripts.utils import set_matplotlib_settings, create_layers_mask
from matplotlib.widgets import Button
from pysrc import *
from tqdm import tqdm
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)


simIDs = [117, 110]
sim_dir = repository_dir + "/simulations"
sim_paths = [sim_dir + f"/simID{simID}" for simID in simIDs]
frame_dir = os.path.join(repository_dir, "images", f"beta_frames_{'_'.join(map(str, simIDs))}")
os.makedirs(frame_dir, exist_ok=True)

step_rate = 10
dikes = [DikeData(sim_path, step_rate=step_rate) for sim_path in sim_paths]
timesteps = len(dikes[0].timesteps)
xlim = (-30, -5)
blim = (0, 1.0)
hlim = (0, 1.7)


for it in range(timesteps):
    print(f"Rendering frame {it+1}/{timesteps}", end="\r")
    frame_path = frame_dir + f"/frame_{it:05d}.png"
    fig = plt.figure(figsize=(3*len(simIDs), 5.5), layout="constrained")
    fig.suptitle(f"Time = {(dikes[0].data[it].time / 3600.0):.2f} h")
    col1 = 7
    col2 = 1
    gs = fig.add_gridspec(1, 2*(col1+col2))
    axs = [fig.add_subplot(gs[i*(col1+col2):i*(col1+col2)+col1]) for i in range(len(simIDs))]
    caxs = [fig.add_subplot(gs[i*(col1+col2)+col1:(i+1)*(col1+col2)]) for i in range(len(simIDs))]
    
    for ax, cax, dike in zip(axs, caxs, dikes):
        data = dike.data[it]
        ax.set_title(r"\bf Crystal concentration")
        ax.set_xlabel(r"Halfwidth (m)")
        ax.set_ylabel(r"Depth (km)")
        ax.grid(which='major', linestyle='-', linewidth=0.3, color="#141414")
        ax.grid(which='minor', linestyle='-', linewidth=0.2, color="#141414")
        ax.minorticks_on()
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(xlim)
        ax.set_xlim(hlim)

        xc = data.xc / 1e3
        yc = data.yc
        xb = data.xb / 1e3
        yb = data.yb
        dx = xb[1] - xb[0]
        hw = data.halfwidth
        hwb = (0.5*(hw[0:-1] + hw[1:])).tolist()
        hwb = np.array([hw[0], *hwb, hw[-1]])
        hwb[hwb < 1e-4] = 1e-4

        x1db = xb
        x2db, _ = np.meshgrid(x1db, yb, indexing='ij')
        y2db = np.array([
            a*yb for a in hwb
        ])

        x2dc, _ = np.meshgrid(xc, yb, indexing='ij')
        y2dc = np.array([
            a*yb for a in hw
        ])

        beta = data.beta
        pcmB = ax.pcolormesh(y2db, x2db, beta, shading='flat', cmap='jet', vmin=blim[0], vmax=blim[1])
        cbB = fig.colorbar(pcmB, cax=cax)

        
    plt.savefig(frame_path, dpi=600)
    plt.close(fig)
plt.show()
