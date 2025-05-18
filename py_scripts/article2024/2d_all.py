import sys
from pathlib import Path
repository_dir = Path.cwd()
sim_dir = repository_dir / "simulations"
sys.path.append(str(repository_dir))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_scripts.utils import set_matplotlib_settings, create_layers_mask, create_save_button
from matplotlib.widgets import Button
from pysrc import *
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)

# simID = 110
# step_rate = 10
# timestep = 100
# sim_path = sim_dir / f"simID{simID}"
# dike = DikeData(str(sim_path), step_rate=step_rate)
# data = dike.data[timestep]
# xind = [50, 150, 225]
# colors = ["k", "k", "k"]
# linestyles = ["-", "--", "-."]
# markers = ["None", "None", "None"]
# nx = len(xind)
# xlim = (-30, -5)
simID = 157
step_rate = 10
timestep = 360
sim_path = sim_dir / f"simID{simID}"
dike = DikeData(str(sim_path), step_rate=step_rate)
data = dike.data[timestep]
xind = [50, 150, 215]
# xind = [50, 150, 225]
colors = ["k", "k", "k"]
linestyles = ["-", "--", "-."]
markers = ["None", "None", "None"]
nx = len(xind)
xlim = (-30, -5)

print(f"data.time = {data.time} s.")

# Create figure with new column for viscosity
fig = plt.figure(figsize=(9, 5.5), layout="constrained")
row1 = 7
col1 = 7
col2 = 1
col3 = 10
gs = fig.add_gridspec(1, 3*(col1+col2)+col3)  # Adjusted the grid to include a third column
gs1 = gs[0, 0:3*(col1+col2)].subgridspec(1, 3*(col1+col2))
gs2 = gs[0, 3*(col1+col2):].subgridspec(3, 1)

axs1 = []
caxs1 = []
for i in range(3):
    nc = i*(col1 + col2)
    ax = fig.add_subplot(gs1[nc:nc+col1])
    cax = fig.add_subplot(gs1[nc+col1:nc+col1+col2])
    axs1.append(ax)
    caxs1.append(cax)

axT = axs1[0]
caxT = caxs1[0]
axT.set_title(r"\bf Temperature (C$^\circ$)")
axT.set_xlabel(r"Halfwidth (m)")
axT.set_ylabel(r"Depth (km)")
axT.set_ylim(xlim)
axT.grid()
axT.spines['top'].set_visible(False)
axT.spines['right'].set_visible(False)

axB = axs1[1]
caxB = caxs1[1]
axB.set_title(r"\bf Crystal concentration")
axB.set_xlabel(r"Halfwidth (m)")
axT.sharey(axB)
axB.set_ylim(xlim)
axB.set_xlim([0, 1])
axB.grid()
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)

axV = axs1[2]
caxV = caxs1[2]
axV.set_title(r"\bf Viscosity ($\log_{10}$ Pa$\cdot$s)")
axV.set_xlabel(r"Halfwidth (m)")
axV.set_ylabel(r"Depth (km)")
axV.set_ylim(xlim)
axV.grid()
axV.spines['top'].set_visible(False)
axV.spines['right'].set_visible(False)

axU1ds = [fig.add_subplot(gs2[-ix-1, 0]) for ix in range(nx)]
axU1ds[-1].set_title(r"\bf Vertical velocity")
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

beta = np.where(data.open_mask[:, None], data.beta, np.nan)
pcmB = axB.pcolormesh(y2db, x2db, beta, shading='flat', cmap='jet')
cbB = fig.colorbar(pcmB, cax=caxB)

viscosity = np.where(data.open_mask[:, None], np.log10(data.viscosity), np.nan)
pcmV = axV.pcolormesh(y2db, x2db, viscosity, shading='flat', cmap='jet')  # Log scale for viscosity
cbV = fig.colorbar(pcmV, cax=caxV)

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


# Add time slider at the bottom
ax_slider = fig.add_axes([0.15, 0.02, 0.7, 0.03])
slider = Slider(ax_slider, "Timestep", 0, len(dike.data) - 1, valinit=timestep, valstep=1)

def update(val):
    step = int(slider.val)
    frame = dike.data[step]

    # Update 2D fields
    T = np.where(frame.open_mask[:, None], frame.temperature, np.nan)
    beta = np.where(frame.open_mask[:, None], frame.beta, np.nan)
    viscosity = np.where(frame.open_mask[:, None], np.log10(frame.viscosity), np.nan)

    pcmT.set_array(T.ravel())
    pcmB.set_array(beta.ravel())
    pcmV.set_array(viscosity.ravel())

    pcmT.set_clim(np.nanmin(T), np.nanmax(T))
    pcmB.set_clim(np.nanmin(beta), np.nanmax(beta))
    pcmV.set_clim(np.nanmin(viscosity), np.nanmax(viscosity))

    cbT.update_normal(pcmT)
    cbB.update_normal(pcmB)
    cbV.update_normal(pcmV)

    # # Update velocity plots
    # As = [frame.A[xi, :] for xi in xind]
    # Cs = [frame.C[xi, :] for xi in xind]
    # velocity_arrays = [velocity(As[i], Cs[i]) for i in range(nx)]

    # for ix in range(nx):
    #     axU1ds[ix].lines[0].set_ydata(velocity_arrays[ix])
    
    fig.canvas.draw_idle()

slider.on_changed(update)


savepath = repository_dir / "images/article2024" / f"2d_T_b_mu_simID{simID}_timestep{timestep*step_rate}"
savepath.parent.mkdir(parents=True, exist_ok=True)
bsave = create_save_button(str(savepath), axs=[ax_slider])
plt.show()
