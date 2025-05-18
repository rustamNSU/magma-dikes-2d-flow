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
step_rate = 10
timestep = 100
sim_path = sim_dir + f"/simID{simID}"
dike = DikeData(sim_path, step_rate=step_rate)
data = dike.data[timestep]
xind = [50, 150, 225]
colors = ["k", "k", "k"]
linestyles = ["-", "--", "-."]
markers = ["None", "None", "None"]
nx = len(xind)
xlim = (-30, -5)

fig = plt.figure(figsize=(6.5, 5.5), layout="constrained")
row1 = 7
col1 = 7
col2 = 1
col3 = 10
gs = fig.add_gridspec(1, 2*(col1+col2)+col3)
gs1 = gs[0, 0:2*(col1+col2)].subgridspec(1, 2*(col1+col2))
gs2 = gs[0, 2*(col1+col2):].subgridspec(3, 1)

axT = fig.add_subplot(gs1[0:col1])
caxT = fig.add_subplot(gs1[col1:col1+col2])

nc = col1+col2
axB = fig.add_subplot(gs1[nc:nc+col1])
caxB = fig.add_subplot(gs1[nc+col1:nc+col1+col2])

axU1ds = [fig.add_subplot(gs2[-ix-1, 0]) for ix in range(nx)]
axU1ds[-1].set_title(r"\bf Vertical velocity")


axT.set_title(r"\bf Temperature (C$^\circ$)")
axT.set_xlabel(r"Halfwidth (m)")
axT.set_ylabel(r"Depth (km)")
axT.set_ylim(xlim)
axT.grid()
axT.spines['top'].set_visible(False)
axT.spines['right'].set_visible(False)

axB.set_title(r"\bf Crystal concentration")
axB.set_xlabel(r"Halfwidth (m)")
# axB.set_ylabel(r"Depth (km)")
axT.sharey(axB)
axB.set_ylim(xlim)
axB.grid()
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)

for axU1d in axU1ds:
    axU1d.set_xlabel(r"$y$ (m)")
    axU1d.set_ylabel(r"Velocity (m/s)")
    axU1d.grid()

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

T = np.where(data.open_mask[:, None], data.temperature, np.nan)
pcmT = axT.pcolormesh(y2db, x2db, T, shading='flat', cmap='jet')
cbT = fig.colorbar(pcmT, cax=caxT)

beta = data.beta
pcmB = axB.pcolormesh(y2db, x2db, beta, shading='flat', cmap='jet')
cbB = fig.colorbar(pcmB, cax=caxB)

# Velocity field
As = [data.A[xi, :] for xi in xind]
Cs = [data.C[xi, :] for xi in xind]
Ny = 200
Y = np.linspace(0, 1, Ny)
yid = create_layers_mask(yb, Y)

def velocity(A, C):
    Ay = np.array([A[iy] for iy in yid])
    Cy = np.array([C[iy] for iy in yid])
    return Ay*Y*Y + Cy

velocity_arrays = [velocity(As[i], Cs[i]) for i in range(nx)]
common_ylim = [min(np.min(v) for v in velocity_arrays), max(np.max(v) for v in velocity_arrays)]
for ix in range(nx):
    xi = xind[ix]
    x = xb[xi+1]
    axT.axhline(x, lw=2, ls=linestyles[ix], color=colors[ix], marker=markers[ix])
    axB.axhline(x, lw=2, ls=linestyles[ix], color=colors[ix], marker=markers[ix])
    axU1ds[ix].plot(Y * hwb[xi+1], velocity(As[ix], Cs[ix]), lw=2,  color=colors[ix], ls=linestyles[ix], marker=markers[ix])
    axU1ds[ix].set_xlim([0, hwb[xi+1]])
    axU1ds[ix].set_ylim(common_ylim)
    
savepath = repository_dir + "/images/article2024"
if not os.path.exists(savepath): os.makedirs(savepath)

ax_button = plt.axes([0.7, 0.15, 0.2, 0.075])
def call(event):
    print("button clicked")
    ax_button.set_visible(False)
    plt.savefig(savepath + f"/vertical_simID{simID}_timestep{step_rate*timestep}.png", dpi=600)
    
button = Button(ax_button, 'Save image')
button.on_clicked(call)
plt.show()
