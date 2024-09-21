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

simIDs = [1]
simLegends = [str(simID) for simID in simIDs]
wlim = (0, 3)
plim = (0, 700)
Tlim = (600, 920)
simPaths = [sim_dir + f"/simID{simID}" for simID in simIDs]
dikes = [DikeData(sim_path, step_rate=1) for sim_path in simPaths]
xlim = (
    min(data.xmin for data in dikes),
    max(data.xmax for data in dikes)
)


fig = plt.figure(figsize=(12, 8), constrained_layout=True)
nrows = 4
ncols = 2
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


laxWs   = []
laxPs   = []
laxTs   = []
laxOPs  = []
laxQxs  = []
laxGs = []
laxRhos = []
laxBetas = []
laxAlphas = []
colors = ['k', 'b', 'r', 'g']
for i in range(len(dikes)):
    data = dikes[i].data
    x = data[0]["xc"]
    xb = data[0]["xb"]
    w = data[0]["width"]
    p = 1e-6*data[0]["pressure"]
    T = np.mean(data[0]["temperature"], axis=1)
    rho = np.mean(data[0]["density"], axis=1)
    beta = np.mean(data[0]["beta"], axis=1)
    alpha = np.mean(data[0]["alpha"], axis=1)
    op = 1e-6*data[0]["overpressure"]
    Qx = data[0]["Qx"]
    G = data[0]["G"]
    laxW,   = axW.plot(x, w, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxP,   = axP.plot(x, p, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxT,   = axT.plot(x, T, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxOP,  = axOP.plot(x, op, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxQx,  = axQx.plot(xb, Qx, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxG, = axG.plot(xb, G, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxRho, = axRho.plot(x, rho, lw=2, ls='-', color=colors[i], label=simLegends[i])
    laxBeta, = axBeta.plot(x, beta, lw=2, ls='-', color=colors[i], label=simLegends[i] + r" ,$\beta$")
    laxAlpha, = axBeta.plot(x, alpha, lw=2, ls='--', color=colors[i], label=simLegends[i] + r" ,$\alpha$")
    laxWs.append(laxW) 
    laxPs.append(laxP) 
    laxTs.append(laxT) 
    laxOPs.append(laxOP)
    laxQxs.append(laxQx)
    laxGs.append(laxG)
    laxRhos.append(laxRho)
    laxBetas.append(laxBeta)
    laxAlphas.append(laxAlpha)

axW.legend().set_draggable(True)
axP.legend().set_draggable(True)
axT.legend().set_draggable(True)
axOP.legend().set_draggable(True)
axQx.legend().set_draggable(True)
axG.legend().set_draggable(True)
axRho.legend().set_draggable(True)
axBeta.legend().set_draggable(True)



ax_time_slider = fig.add_subplot(gs[nrows*nrax:nrows*nrax + nrslider, 0:2])
time_slider = Slider(
    ax=ax_time_slider, 
    label='timesteps', 
    valmin=0, 
    valmax=max(max(data.timesteps) for data in dikes),
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
        rho = np.mean(data["density"], axis=1)
        beta = np.mean(data["beta"], axis=1)
        alpha = np.mean(data["alpha"], axis=1)
        op = 1e-6*data["overpressure"]
        Qx = data["Qx"]
        G = data["G"]
        laxWs[i].set_ydata(w) 
        laxPs[i].set_ydata(p) 
        laxTs[i].set_ydata(T) 
        laxOPs[i].set_ydata(op)
        laxQxs[i].set_ydata(Qx)
        laxGs[i].set_ydata(G)
        laxRhos[i].set_ydata(rho)
        laxBetas[i].set_ydata(beta)
        laxAlphas[i].set_ydata(alpha)
        
        maxQx = max(1e-6, max(Qx))
        axQx.set_ylim(-1e-2*maxQx, maxQx)
        maxG = max(1e-6, max(G))
        axG.set_ylim(-1e-2*maxG, maxG)

time_slider.on_changed(time_update)
plt.show()