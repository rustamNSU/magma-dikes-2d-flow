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
from matplotlib.widgets import Button
from pysrc import *
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)


simID = 110
step_rate = 20
timestep = 100
sim_path = sim_dir + f"/simID{simID}"
dike = DikeData(sim_path, step_rate=step_rate)
data = dike.data[timestep]
xind = [150, 200, 235]
colors = ["k", "k", "k"]
linestyles = ["-", "--", "-."]
markers = ["None", "None", "None"]
nx = len(xind)
xlim = (-30, -5)
Tlim = (200, 900)
blim = (0, 1.0)
T1dlim = (550, 750)

fig = plt.figure(figsize=(10, 5), layout="constrained")
fig.suptitle(f"Time = {(data.time / 3600.0):.2f} h")
col1 = 7
col2 = 1
col3 = 10
gs = fig.add_gridspec(1, 2*(col1+col2)+2*col3)
gs1 = gs[0, 0:2*(col1+col2)].subgridspec(1, 2*(col1+col2))
gs2 = gs[0, 2*(col1+col2):].subgridspec(3, 2)

axT = fig.add_subplot(gs1[0:col1])
caxT = fig.add_subplot(gs1[col1:col1+col2])

nc = col1+col2
axB = fig.add_subplot(gs1[nc:nc+col1])
caxB = fig.add_subplot(gs1[nc+col1:nc+col1+col2])

ax1dT = [fig.add_subplot(gs2[-ix-1, 0]) for ix in range(nx)]
ax1dB = [fig.add_subplot(gs2[-ix-1, 1]) for ix in range(nx)]
ax1dT[-1].set_title(r"\bf Temperature slice")
ax1dB[-1].set_title(r"\bf Crystal concentration slice")


axT.set_title(r"\bf Temperature (C$^\circ$)")
axT.set_xlabel(r"Halfwidth (m)")
axT.set_ylabel(r"Depth (km)")
axT.set_ylim(xlim)
axT.grid(which='major', linestyle='-', linewidth=0.75)
axT.grid(which='minor', linestyle='-', linewidth=0.5)
axT.minorticks_on()
axT.set_axisbelow(True)

axT.spines['top'].set_visible(False)
axT.spines['right'].set_visible(False)

axB.set_title(r"\bf Crystal concentration")
axB.set_xlabel(r"Halfwidth (m)")
axB.set_ylabel(r"Depth (km)")
axT.sharey(axB)
axB.set_ylim(xlim)
axB.grid()
axB.grid(which='major', linestyle='-', linewidth=0.75)
axB.grid(which='minor', linestyle='-', linewidth=0.5)
axB.minorticks_on()
axB.set_axisbelow(True)

axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)

for ax1, ax2 in zip(ax1dT, ax1dB):
    ax1.set_xlabel(r"$y$ (m)")
    ax1.set_ylabel(r"Temperature (C$^\circ$)")
    ax1.grid()
    ax2.set_xlabel(r"$y$ (m)")
    ax2.set_ylabel(r"$\beta$")
    ax2.grid()

xc = data.xc / 1e3
yc = data.yc
xb = data.xb / 1e3
yb = data.yb
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

T = np.where(data.open_mask[:, None], data.temperature, np.nan)
pcmT = axT.pcolormesh(y2db, x2db, T, shading='flat', cmap='jet', vmin=Tlim[0], vmax=Tlim[1])
cbT = fig.colorbar(pcmT, cax=caxT)

beta = data.beta
pcmB = axB.pcolormesh(y2db, x2db, beta, shading='flat', cmap='jet', vmin=blim[0], vmax=blim[1])
cbB = fig.colorbar(pcmB, cax=caxB)

for ix in range(nx):
    xi = xind[ix]
    x = xb[xi+1]
    axT.axhline(x, lw=2, ls=linestyles[ix], color=colors[ix], marker=markers[ix])
    axB.axhline(x, lw=2, ls=linestyles[ix], color=colors[ix], marker=markers[ix])
    Y = yc * hw[xi]
    T1d = T[xi, :]
    B1d = beta[xi, :]
    ax1dT[ix].plot(Y, T1d, lw=2,  color=colors[ix], ls=linestyles[ix], marker=markers[ix])
    ax1dT[ix].set_xlim([0, hw[xi]])
    ax1dT[ix].set_ylim(T1dlim)
    ax1dB[ix].plot(Y, B1d, lw=2,  color=colors[ix], ls=linestyles[ix], marker=markers[ix])
    ax1dB[ix].set_xlim([0, hw[xi]])
    ax1dB[ix].set_ylim(blim)
    
savepath = repository_dir + "/images"
if not os.path.exists(savepath): os.makedirs(savepath)

ax_button = plt.axes([0.7, 0.15, 0.2, 0.075])
def call(event):
    print("button clicked")
    ax_button.set_visible(False)
    plt.savefig(savepath + f"/vertical_T_beta_simID{simID}_timestep{step_rate*timestep}.png", dpi=600)
    
button = Button(ax_button, 'Save image')
button.on_clicked(call)
plt.show()
