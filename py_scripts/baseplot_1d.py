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
from pysrc import *
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)

simIDs = [110, 113, 115, 116]
simLegends = [f"simID{simID}" for simID in simIDs]
# simIDs = [111, 110, 112]
# simLegends = [
#     r"$\tau_0 = 10^{-5}$",
#     r"$\tau_0 = 10^{-6}$",
#     r"$\tau_0 = 10^{-7}$",
# ]

# Plot limits
wlim = (0, 5)
plim = (0, 900)
Tlim = (300, 900)

# Create list of simulation paths and load DikeData for each simulation
simPaths = [os.path.join(sim_dir, f"simID{simID}") for simID in simIDs]
dikes = [DikeData(sim_path, step_rate=10) for sim_path in simPaths]

# Determine overall x-limits from all dikes
xlim = (
    min(dike.xmin for dike in dikes),
    max(dike.xmax for dike in dikes)
)

# Create figure and a gridspec layout
fig = plt.figure(figsize=(17, 8), constrained_layout=True)
nrows = 4
ncols = 3
nrax = 10
nrslider = 1
gs = fig.add_gridspec(nrows=(nrows*nrax + nrslider), ncols=ncols)

# Create subplots for various quantities
axW = fig.add_subplot(gs[0:nrax, 0])
axW.set_xlabel(r"Depth (m)")
axW.set_ylabel(r"Width (m)")
axW.set_xlim(xlim)
axW.set_ylim(wlim)
axW.grid()

axP = fig.add_subplot(gs[0:nrax, 1])
axP.set_xlabel(r"Depth (m)")
axP.set_ylabel(r"Pressure (MPa)")
axP.set_xlim(xlim)
axP.set_ylim(plim)
axP.grid()

axMu = fig.add_subplot(gs[0:nrax, 2])
axMu.set_xlabel(r"Depth (m)")
axMu.set_ylabel(r"Viscosity (Pa$\cdot$s)")
axMu.set_xlim(xlim)
axMu.grid()

axT = fig.add_subplot(gs[nrax:2*nrax, 0])
axT.set_xlabel(r"Depth (m)")
axT.set_ylabel(r"Averaged temperature (°C)")
axT.set_xlim(xlim)
axT.set_ylim(Tlim)
axT.grid()

axOP = fig.add_subplot(gs[nrax:2*nrax, 1])
axOP.set_xlabel(r"Depth (m)")
axOP.set_ylabel(r"Overpressure (MPa)")
axOP.set_xlim(xlim)
axOP.set_ylim([-10, 20])
axOP.grid()

axQx = fig.add_subplot(gs[2*nrax:3*nrax, 0])
axQx.set_xlabel(r"Depth (m)")
axQx.set_ylabel(r"Total flux rate (m$^2$/s)")
axQx.set_xlim(xlim)
axQx.grid()

axG = fig.add_subplot(gs[2*nrax:3*nrax, 1])
axG.set_xlabel(r"Depth (m)")
axG.set_ylabel(r"Pressure gradient (Pa/m)")
axG.set_xlim(xlim)
axG.grid()

axRho = fig.add_subplot(gs[3*nrax:4*nrax, 0])
axRho.set_xlabel(r"Depth (m)")
axRho.set_ylabel(r"Averaged density (kg/m$^3$)")
axRho.set_xlim(xlim)
axRho.set_ylim([1900, 2600])
axRho.grid()

axBeta= fig.add_subplot(gs[3*nrax:4*nrax, 1])
axBeta.set_xlabel(r"Depth (m)")
axBeta.set_ylabel(r"Averaged crystal content")
axBeta.set_xlim(xlim)
axBeta.set_ylim([0, 1])
axBeta.grid()

axTw = fig.add_subplot(gs[nrax:2*nrax, 2])
axTw.set_xlabel(r"Depth (m)")
axTw.set_ylabel(r"Wall temperature (°C)")
axTw.set_xlim(xlim)
axTw.set_ylim(Tlim)
axTw.grid()

axHeat = fig.add_subplot(gs[2*nrax:3*nrax, 2])
axHeat.set_xlabel(r"Depth (m)")
axHeat.set_ylabel(r"Heat flux (cooling) (W/m$^2$)")
axHeat.set_xlim(xlim)
axHeat.grid()

axShear = fig.add_subplot(gs[3*nrax:4*nrax, 2])
axShear.set_xlabel(r"Depth (m)")
axShear.set_ylabel(r"Width average shear heat rate (W/m$^2$)")
axShear.set_xlim(xlim)
axShear.grid()

# Share x-axis among subplots for alignment
axP.sharex(axW)
axMu.sharex(axW)
axT.sharex(axW)
axOP.sharex(axW)
axQx.sharex(axW)
axG.sharex(axW)
axRho.sharex(axW)
axBeta.sharex(axW)
axTw.sharex(axW)
axHeat.sharex(axW)
axShear.sharex(axW)

# Create lists to hold plot lines for each simulation
laxWs    = []
laxPs    = []
laxMus    = []
laxTs    = []
laxTws   = []
laxOPs   = []
laxQxs   = []
laxGs    = []
laxRhos  = []
laxBetas = []
laxAlphas= []
laxHeats = []
laxShears= []
colors = ['k', 'b', 'r', 'g']

for i, dike in enumerate(dikes):
    # Use attribute access for the DataFrame dataclass
    frame = dike.data[0]
    x = frame.xc
    xb = frame.xb
    w = frame.width
    p = 1e-6 * frame.pressure
    mu = np.mean(frame.viscosity, axis=1)
    T = np.mean(frame.temperature, axis=1)
    Tw = frame.wall_temperature
    rho = np.mean(frame.density, axis=1)
    beta = np.mean(frame.beta, axis=1)
    alpha = np.mean(frame.alpha, axis=1)
    op = 1e-6 * frame.overpressure
    Qx = frame.Qx
    G = frame.G

    laxW, = axW.plot(x, w, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxP, = axP.plot(x, p, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxMu, = axMu.plot(x, mu, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxT, = axT.plot(x, T, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxTw, = axTw.plot(x, Tw, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxOP, = axOP.plot(x, op, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxQx, = axQx.plot(xb, Qx, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxG, = axG.plot(xb, G, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxRho, = axRho.plot(x, rho, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxBeta, = axBeta.plot(x, beta, lw=2, ls='-', color=colors[i], label=simLegends[i] + r", $\beta$")
    laxAlpha, = axBeta.plot(x, alpha, lw=2, ls='--', color=colors[i], label=simLegends[i] + r", $\alpha$")
    laxHeat, = axHeat.plot(x, frame.heatflux, lw=2, ls='-', color=colors[i], label=simLegends[i])
    
    laxWs.append(laxW)
    laxPs.append(laxP)
    laxMus.append(laxMu)
    laxTs.append(laxT)
    laxTws.append(laxTw)
    laxOPs.append(laxOP)
    laxQxs.append(laxQx)
    laxGs.append(laxG)
    laxRhos.append(laxRho)
    laxBetas.append(laxBeta)
    laxAlphas.append(laxAlpha)
    laxHeats.append(laxHeat)

# Make legends draggable
for ax in [axW, axMu, axP, axT, axTw, axOP, axQx, axG, axRho, axBeta, axHeat]:
    ax.legend().set_draggable(True)

# Create a slider for time update
ax_time_slider = fig.add_subplot(gs[nrows*nrax:nrows*nrax + nrslider, 0:2])
time_slider = Slider(
    ax=ax_time_slider,
    label='Timesteps',
    valmin=0,
    valmax=min(len(dike.timesteps) for dike in dikes) - 1,
    valinit=0,
    valstep=1
)

def time_update(val):
    t = int(time_slider.val)
    fig.suptitle('Time = {0:.2f} h'.format(dikes[0].data[t].time / 3600))
    for i, dike in enumerate(dikes):
        frame = dike.data[t]
        mu = np.mean(frame.viscosity, axis=1)
        laxWs[i].set_ydata(frame.width)
        laxPs[i].set_ydata(1e-6 * frame.pressure)
        laxMus[i].set_ydata(mu)
        laxTs[i].set_ydata(np.mean(frame.temperature, axis=1))
        laxTws[i].set_ydata(frame.wall_temperature)
        laxOPs[i].set_ydata(1e-6 * frame.overpressure)
        laxQxs[i].set_ydata(frame.Qx)
        laxGs[i].set_ydata(frame.G)
        laxRhos[i].set_ydata(np.mean(frame.density, axis=1))
        laxBetas[i].set_ydata(np.mean(frame.beta, axis=1))
        laxAlphas[i].set_ydata(np.mean(frame.alpha, axis=1))
        laxHeats[i].set_ydata(frame.heatflux)
        
        # Optionally update axis limits if needed:
        axQx.set_ylim(-1e-2 * max(1e-6, np.nanmax(frame.Qx)), np.nanmax(frame.Qx))
        axG.set_ylim(-1e-2 * max(1e-6, np.nanmax(frame.G)), np.nanmax(frame.G))
        axOP.set_ylim(min(axOP.get_ylim()[0], np.nanmin(1e-6*frame.overpressure)),
                     max(axOP.get_ylim()[1], np.nanmax(1e-6*frame.overpressure)))
        axHeat.set_ylim(min(axHeat.get_ylim()[0], np.nanmin(frame.heatflux)),
                       max(axHeat.get_ylim()[1], np.nanmax(frame.heatflux)))
        axHeat.set_ylim(min(axHeat.get_ylim()[0], np.nanmin(frame.heatflux)),
                       max(axHeat.get_ylim()[1], np.nanmax(frame.heatflux)))
        if i == 0: axMu.set_ylim(0, 0.1)
        axMu.set_ylim(0, max(axMu.get_ylim()[1], np.nanmax(mu)))
        
    fig.canvas.draw_idle()

time_slider.on_changed(time_update)
plt.show()