import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from py_scripts.utils import set_matplotlib_settings, create_layers_mask
from matplotlib.widgets import Button
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib.patches import Rectangle
from pysrc import *
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)


simID = 110
step_rate = 10
timestep = 100
sim_path = sim_dir + f"/simID{simID}"
dike = DikeData(sim_path, step_rate=step_rate)
data = dike.data[timestep]
xind = [50, 150, 220]
colors = ["k", "k", "k"]
linestyles = ["-", "--", "-."]
markers = ["None", "None", "None"]
nx = len(xind)
xlim = (-30, -5)
Tlim = (400, 900)
blim = (0, 1.0)
Mlim = (3, 10)


# simID = 113
# step_rate = 10
# timestep = 100
# sim_path = sim_dir + f"/simID{simID}"
# dike = DikeData(sim_path, step_rate=step_rate)
# data = dike.data[timestep]
# xind = [50, 150, 220]
# colors = ["k", "k", "k"]
# linestyles = ["-", "--", "-."]
# markers = ["None", "None", "None"]
# nx = len(xind)
# xlim = (-30, -5)
# Tlim = (400, 900)
# blim = (0, 1.0)
# Mlim = (3, 14)


# simID = 154
# step_rate = 10
# timestep = 100
# sim_path = sim_dir + f"/simID{simID}"
# dike = DikeData(sim_path, step_rate=step_rate)
# data = dike.data[timestep]
# xind = [50, 150, 220]
# colors = ["k", "k", "k"]
# linestyles = ["-", "--", "-."]
# markers = ["None", "None", "None"]
# nx = len(xind)
# xlim = (-30, -5)
# Tlim = (400, 900)
# blim = (0, 0.4)
# Mlim = (3, 10)


fig = plt.figure(figsize=(6.5, 5.5), layout="constrained")
row1 = 7
col1 = 7
col2 = 0
col3 = 10
nc = col1+col2
gs = fig.add_gridspec(1, 3*nc+col3)
gs1 = gs[0, 0:3*nc].subgridspec(1, 3*nc)
gs2 = gs[0, 3*nc:].subgridspec(3, 1)

axT = fig.add_subplot(gs1[0:col1])
caxT = inset_axes(
    axT,
    width="10%",        # thickness
    height="40%",      # vertical height
    loc='lower right',
    bbox_to_anchor=(-0.3, 0.02, 1, 1),
    bbox_transform=axT.transAxes,
    borderpad=0)
caxT.yaxis.set_ticks_position('right')
caxT.yaxis.set_label_position('right')
rect = Rectangle(
    (0.55, 0.01),       # lower-left corner
    0.42, 0.42,  # size of rectangle
    transform=axT.transAxes,
    color='white',
    zorder=10,
    clip_on=False
)
axT.add_patch(rect)


axB = fig.add_subplot(gs1[nc:nc+col1])
caxB = inset_axes(
    axB,
    width="10%",        # thickness
    height="40%",      # vertical height
    loc='lower right',
    bbox_to_anchor=(-0.3, 0.02, 1, 1),
    bbox_transform=axB.transAxes,
    borderpad=0)
caxB.yaxis.set_ticks_position('right')
caxB.yaxis.set_label_position('right')
rect = Rectangle(
    (0.55, 0.01),       # lower-left corner
    0.42, 0.42,  # size of rectangle
    transform=axB.transAxes,
    color='white',
    zorder=10,
    clip_on=False
)
axB.add_patch(rect)


axM = fig.add_subplot(gs1[2*nc:2*nc+col1])
caxM = inset_axes(
    axM,
    width="10%",        # thickness
    height="40%",      # vertical height
    loc='lower right',
    bbox_to_anchor=(-0.3, 0.02, 1, 1),
    bbox_transform=axM.transAxes,
    borderpad=0)
caxM.yaxis.set_ticks_position('right')
caxM.yaxis.set_label_position('right')
rect = Rectangle(
    (0.55, 0.01),       # lower-left corner
    0.42, 0.42,  # size of rectangle
    transform=axM.transAxes,
    color='white',
    zorder=10,
    clip_on=False
)
axM.add_patch(rect)

axU1ds = [fig.add_subplot(gs2[-ix-1, 0]) for ix in range(nx)]
axU1ds[-1].set_title(r"Vertical velocity")


axT.set_title(r"Temperature (C$^\circ$)")
axT.set_xlabel(r"Halfwidth (m)")
axT.set_ylabel(r"Depth (km)")
axT.set_ylim(xlim)
axT.grid(which='major', linestyle='-', linewidth=0.75)
axT.grid(which='minor', linestyle='-', linewidth=0.5)
axT.minorticks_on()
axT.set_axisbelow(True)
axT.spines['top'].set_visible(False)
axT.spines['right'].set_visible(False)

axB.set_title(r"Crystal concentration")
axB.set_xlabel(r"Halfwidth (m)")
axB.sharey(axT)
axB.set_ylim(xlim)
axB.grid()
axB.grid(which='major', linestyle='-', linewidth=0.75)
axB.grid(which='minor', linestyle='-', linewidth=0.5)
axB.minorticks_on()
axB.set_axisbelow(True)
axB.tick_params(axis='y', left=False, labelleft=False)
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)


axM.set_title(r"Viscosity, $\log_{10}$ (Pa$\cdot$s)")
axM.set_xlabel(r"Halfwidth (m)")
axM.sharey(axT)
axM.set_ylim(xlim)
axM.grid()
axM.grid(which='major', linestyle='-', linewidth=0.75)
axM.grid(which='minor', linestyle='-', linewidth=0.5)
axM.minorticks_on()
axM.set_axisbelow(True)
axM.tick_params(axis='y', left=False, labelleft=False)
axM.spines['top'].set_visible(False)
axM.spines['right'].set_visible(False)

for ax in axU1ds:
    ax.set_xlabel(r"$y$ (m)")
    ax.set_ylabel(r"Velocity (m/s)")
    ax.grid()
    
    # ax.yaxis.set_major_formatter(FormatStrFormatter("%.1e"))  # full scientific notation
    # formatter = ScalarFormatter(useMathText=False)
    # formatter.set_scientific(True)
    # formatter.set_powerlimits((-3, 3))  # use sci notation for small numbers
    # formatter.set_useOffset(False)
    # formatter.set_useLocale(False)
    # ax.yaxis.set_major_formatter(formatter)
    # ax.ticklabel_format(style='sci', axis='y', scilimits=(-3, 3))

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
pcmT = axT.pcolormesh(y2db, x2db, T, shading='flat', cmap='jet', vmin=Tlim[0], vmax=Tlim[1])
cbT = fig.colorbar(pcmT, cax=caxT)
cbT.ax.set_zorder(11) 

beta = data.beta
pcmB = axB.pcolormesh(y2db, x2db, beta, shading='flat', cmap='jet', vmin=blim[0], vmax=blim[1])
cbB = fig.colorbar(pcmB, cax=caxB)

beta = data.beta
pcmB = axB.pcolormesh(y2db, x2db, beta, shading='flat', cmap='jet', vmin=blim[0], vmax=blim[1])
cbB = fig.colorbar(pcmB, cax=caxB)

viscosity = np.where(data.open_mask[:, None], np.log10(data.viscosity), np.nan)
pcmM = axM.pcolormesh(y2db, x2db, viscosity, shading='flat', cmap='jet', vmin=Mlim[0], vmax=Mlim[1])  # Log scale for viscosity
cbM = fig.colorbar(pcmM, cax=caxM)


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
    axM.axhline(x, lw=2, ls=linestyles[ix], color=colors[ix], marker=markers[ix])
    axU1ds[ix].plot(Y * hwb[xi+1], velocity(As[ix], Cs[ix]), lw=2,  color=colors[ix], ls=linestyles[ix], marker=markers[ix])
    axU1ds[ix].set_xlim([0, hwb[xi+1]])
    axU1ds[ix].set_ylim(common_ylim)
    
savepath = repository_dir + "/images/article2024"
if not os.path.exists(savepath): os.makedirs(savepath)

ax_button = plt.axes([0.7, 0.15, 0.2, 0.075])
def call(event):
    print("button clicked")
    ax_button.set_visible(False)
    plt.savefig(savepath + f"/vertical_review_simID{simID}_timestep{step_rate*timestep}.png", dpi=600)
    
button = Button(ax_button, 'Save image')
button.on_clicked(call)
plt.show()
