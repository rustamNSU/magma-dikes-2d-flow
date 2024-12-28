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
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)

# simIDs = [3, 4, 5]
# simLegends = [r"$\Delta x = 100$", r"$\Delta x = 50$", r"$\Delta x = 25$"]
# simIDs = [10, 11, 12]
# simLegends = [r"$N_y = 2$", r"$N_y = 10$", r"$N_y = 30$"]
# simIDs = [12, 13]
# simLegends = [r"$\Delta x = 100$", r"$\Delta x = 50$"]
# simIDs = [12, 14]
# simLegends = [r"$k_m = 2$", r"$k_m = 20000$"]
# simIDs = [10, 11, 12]
# simLegends = [r"$K_{Ic} = 1$ MPa$\cdot$m$^{-1/2}$", r"$K_{Ic} = 100$ MPa$\cdot$m$^{-1/2}$", r"$K_{Ic} = 1000$ MPa$\cdot$m$^{-1/2}$"]
# simIDs = [12, 13]
# simLegends = [r"$K_{Ic} = 1000$ MPa$\cdot$m$^{-1/2}$, $\Delta x = 50$ m", r"$K_{Ic} = 1000$ MPa$\cdot$m$^{-1/2}$, $\Delta x = 25$ m"]
# simIDs = [13, 14]
# simLegends = [r"$K_{Ic} = 1000$ MPa$\cdot$m$^{-1/2}$, $N_{coh} = 6$", r"$K_{Ic} = 1000$ MPa$\cdot$m$^{-1/2}$, $N_{coh} = 10$"]
# simIDs = [100, 101, 102]
# simLegends = [str(simID) for simID in simIDs]

# simIDs = [11, 21]
# simLegends = [r"$6.23$ wt.$\%$, 2d", r"$6.23$ wt.$\%$, 1d"]
simIDs = [102, 101, 100]
simLegends = [r"$3.85$ wt.$\%$", r"$6.18$ wt.$\%$", r"$9.57$ wt.$\%$"]
wlim = (0, 5)
plim = (0, 900)
Tlim = (300, 900)
simPaths = [sim_dir + f"/simID{simID}" for simID in simIDs]
dikes = [DikeData(sim_path, step_rate=1) for sim_path in simPaths]
xlim = (
    min(data.xmin for data in dikes),
    max(data.xmax for data in dikes)
)


fig = plt.figure(figsize=(17, 8), constrained_layout=True)
nrows = 4
ncols = 3
nrax = 10
nrslider = 1
gs = fig.add_gridspec(nrows=(nrows*nrax + nrslider), ncols=ncols)
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

axT = fig.add_subplot(gs[nrax:2*nrax, 0])
axT.set_xlabel(r"Depth (m)")
axT.set_ylabel(r"Averaged temperature (C$^\circ$)")
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
axQx.set_ylabel(r"Total flux rate (m$^2/$s)")
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

axBeta = fig.add_subplot(gs[3*nrax:4*nrax, 1])
axBeta.set_xlabel(r"Depth (m)")
axBeta.set_ylabel(r"Averaged crystal content")
axBeta.set_xlim(xlim)
axBeta.set_ylim([0, 1])
axBeta.grid()

axTw = fig.add_subplot(gs[nrax:2*nrax, 2])
axTw.set_xlabel(r"Depth (m)")
axTw.set_ylabel(r"Wall temperature (C$^\circ$)")
axTw.set_xlim(xlim)
axTw.set_ylim(Tlim)
axTw.grid()

axHeat = fig.add_subplot(gs[2*nrax:3*nrax, 2])
axHeat.set_xlabel(r"Depth (m)")
axHeat.set_ylabel(r"Heat flux (cooling) (W$/$m$^2$)")
axHeat.set_xlim(xlim)
axHeat.grid()

axShear = fig.add_subplot(gs[3*nrax:4*nrax, 2])
axShear.set_xlabel(r"Depth (m)")
axShear.set_ylabel(r"Width average shear heat rate (W$/$m$^2$)")
axShear.set_xlim(xlim)
axShear.grid()

axP.sharex(axW)
axT.sharex(axW)
axOP.sharex(axW)
axQx.sharex(axW)
axG.sharex(axW)
axRho.sharex(axW)
axBeta.sharex(axW)
axTw.sharex(axW)
axHeat.sharex(axW)
axShear.sharex(axW)


laxWs   = []
laxPs   = []
laxTs   = []
laxTws   = []
laxOPs  = []
laxQxs  = []
laxGs = []
laxRhos = []
laxBetas = []
laxAlphas = []
laxHeats = []
laxShears = []
colors = ['k', 'b', 'r', 'g']
for i in range(len(dikes)):
    data = dikes[i].data
    x = data[0]["xc"]
    xb = data[0]["xb"]
    w = data[0]["width"]
    p = 1e-6*data[0]["pressure"]
    T = np.mean(data[0]["temperature"], axis=1)
    Tw = data[0]["Twall"]
    rho = np.mean(data[0]["density"], axis=1)
    beta = np.mean(data[0]["beta"], axis=1)
    alpha = np.mean(data[0]["alpha"], axis=1)
    op = 1e-6*data[0]["overpressure"]
    Qx = data[0]["Qx"]
    G = data[0]["G"]
    laxW,   = axW.plot(x, w, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxP,   = axP.plot(x, p, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxT,   = axT.plot(x, T, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxTw,   = axTw.plot(x, Tw, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxOP,  = axOP.plot(x, op, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxQx,  = axQx.plot(xb, Qx, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxG, = axG.plot(xb, G, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxRho, = axRho.plot(x, rho, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxBeta, = axBeta.plot(x, beta, lw=2, ls='-', color=colors[i], label=simLegends[i] + r" ,$\beta$")
    laxAlpha, = axBeta.plot(x, alpha, lw=2, ls='--', color=colors[i], label=simLegends[i] + r" ,$\alpha$")
    laxHeat, = axHeat.plot(x, data[0]["heatflux"], lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxShear, = axShear.plot(x, data[0]["shear_heat_rate"], lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxWs.append(laxW) 
    laxPs.append(laxP) 
    laxTs.append(laxT)
    laxTws.append(laxTw)
    laxOPs.append(laxOP)
    laxQxs.append(laxQx)
    laxGs.append(laxG)
    laxRhos.append(laxRho)
    laxBetas.append(laxBeta)
    laxAlphas.append(laxAlpha)
    laxHeats.append(laxHeat)
    laxShears.append(laxShear)

axW.legend().set_draggable(True)
axP.legend().set_draggable(True)
axT.legend().set_draggable(True)
axTw.legend().set_draggable(True)
axOP.legend().set_draggable(True)
axQx.legend().set_draggable(True)
axG.legend().set_draggable(True)
axRho.legend().set_draggable(True)
axBeta.legend().set_draggable(True)
axHeat.legend().set_draggable(True)



ax_time_slider = fig.add_subplot(gs[nrows*nrax:nrows*nrax + nrslider, 0:2])
time_slider = Slider(
    ax=ax_time_slider, 
    label='timesteps', 
    valmin=0, 
    valmax=min(max(data.timesteps) for data in dikes),
    valinit=0, 
    valstep=1)


def time_update(val):
    t = time_slider.val
    fig.suptitle('time = {0:.2f} h'.format(dikes[0].data[t]["time"]/3600))
    for i in range(len(dikes)):
        data = dikes[i].data[t]
        w = data["width"]
        p = 1e-6*data["pressure"]
        T = np.mean(data["temperature"], axis=1)
        Tw = data["Twall"]
        rho = np.mean(data["density"], axis=1)
        beta = np.mean(data["beta"], axis=1)
        alpha = np.mean(data["alpha"], axis=1)
        op = 1e-6*data["overpressure"]
        Qx = data["Qx"]
        G = data["G"]
        heatflux = data["heatflux"]
        shear_heat_rate = data["shear_heat_rate"]
        
        laxWs[i].set_ydata(w) 
        laxPs[i].set_ydata(p) 
        laxTs[i].set_ydata(T)
        laxTws[i].set_ydata(Tw)
        laxOPs[i].set_ydata(op)
        laxQxs[i].set_ydata(Qx)
        laxGs[i].set_ydata(G)
        laxRhos[i].set_ydata(rho)
        laxBetas[i].set_ydata(beta)
        laxAlphas[i].set_ydata(alpha)
        laxAlphas[i].set_ydata(alpha)
        laxHeats[i].set_ydata(heatflux)
        laxShears[i].set_ydata(shear_heat_rate)
        
        maxQx = max(1e-6, max(Qx))
        axQx.set_ylim(-1e-2*maxQx, maxQx)
        maxG = max(1e-6, max(G))
        axG.set_ylim(-1e-2*maxG, maxG)
        axOP.set_ylim(min(-1e-6, axOP.get_ylim()[0],  min(op)), max(axOP.get_ylim()[1], max(op)))
        axHeat.set_ylim(min(-1e-6, axHeat.get_ylim()[0], min(heatflux)), max(axHeat.get_ylim()[1], max(heatflux)))
        axShear.set_ylim(min(-1e-6, axShear.get_ylim()[0], min(shear_heat_rate)), max(axShear.get_ylim()[1], max(shear_heat_rate)))

time_slider.on_changed(time_update)
plt.show()