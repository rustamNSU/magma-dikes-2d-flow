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

FIG_SIZE = [10, 8]
set_matplotlib_settings()
simID = 302
sim_path = sim_dir + f"/simID{simID}"
wlim = (0, 1)
plim = (0, 700)
Tlim = (600, 1100)
xlim = (-30000, -10000.0)
# xlim = (-30000, 0)

dike = DikeData(sim_path, step_rate=10)
data = dike.data
fig = plt.figure(figsize=FIG_SIZE, constrained_layout=True)

n1d = 8
n2d = 30
ncb = 1
gs = fig.add_gridspec(5, n1d+n2d+ncb)
axW = fig.add_subplot(gs[0, n1d:-ncb])
axT = fig.add_subplot(gs[1, n1d:-ncb])
axU = fig.add_subplot(gs[2, n1d:-ncb])
axB = fig.add_subplot(gs[3, n1d:-ncb])
axMu1d = fig.add_subplot(gs[0, 0:n1d])
axT1d = fig.add_subplot(gs[1, 0:n1d])
axU1d = fig.add_subplot(gs[2, 0:n1d])
axB1d = fig.add_subplot(gs[3, 0:n1d])

caxW = fig.add_subplot(gs[0, -ncb:])
caxT = fig.add_subplot(gs[1, -ncb:])
caxU = fig.add_subplot(gs[2, -ncb:])
caxB = fig.add_subplot(gs[3, -ncb:])
caxW.axis('off')



axW.set_xlabel(r"Depth (m)")
axW.set_ylabel(r"Halfwidth (m)")
axMu1d.set_xlabel(r"Viscosity (Pa$\cdot$s)")
axMu1d.set_ylabel(r"$\xi = y/h$")

axT.set_xlabel(r"Depth (m)")
axT.set_ylabel(r"Temperature (C$^\circ$)")
axT1d.set_xlabel(r"Temperature (C$^\circ$)")
axT1d.set_ylabel(r"$\xi = y/h$")

axU.set_xlabel(r"Depth (m)")
axU.set_ylabel(r"Vertical flux (m$^2$/s)")
axU1d.set_xlabel(r"Velocity-$u$ (m/s)")
axU1d.set_ylabel(r"$\xi = y/h$")

axB.set_xlabel(r"Depth (m)")
axB.set_ylabel(r"Crystal concentration")
axB1d.set_xlabel(r"Crystal concentration")
axB1d.set_ylabel(r"$\xi = y/h$")

axW.set_xlim(xlim)
axW.set_ylim(wlim)
axT.set_xlim(xlim)
axU.set_xlim(xlim)
axB.set_xlim(xlim)

axMu1d.set_ylim([0, 1])
axT1d.set_ylim([0, 1])
axU1d.set_ylim([0, 1])
axB1d.set_ylim([0, 1])

axT1d.set_xlim(Tlim)
axU1d.set_xlim([0, None])
axMu1d.set_xlim([0, None])
axB1d.set_xlim([0, 1])

axW.grid()
axT1d.grid()
axU1d.grid()
axB1d.grid()
axMu1d.grid()

xc = data[0]["xc"]
yc = data[0]["yc"]
xb = data[0]["xb"]
yb = data[0]["yb"]
dx = xb[1]-xb[0]

hw = data[0]["halfwidth"]
T = data[0]["Tmask"]
Tl = data[0]["Tliquidus"]
Ts = data[0]["Tsolidus"]
qx = data[0]["qx"][1:-1, :]
beta = data[0]["beta"]
betaeq = data[0]["betaeq"]

# Velocity field
A = data[0]["A"][0, :]
C = data[0]["C"][0, :]
Ny = 200
Y = np.linspace(0, 1, Ny)
yid = create_layers_mask(yb, Y)

def velocity(A, C):
    Ay = np.array([A[iy] for iy in yid])
    Cy = np.array([C[iy] for iy in yid])
    return Ay*Y*Y + Cy

lW, = axW.plot(xc, hw, lw=3)
lMu1d, = axMu1d.plot(data[0]["viscosity"][0, :], yc, lw=2, color="r")

pcmT = axT.pcolormesh(xb, yb, T.T, shading='flat', cmap='jet')
cbT = fig.colorbar(pcmT, cax=caxT)
lT1d, = axT1d.plot(T[0, :], yc, lw=2, color="r", label=r"$T$")
lTl1d, = axT1d.plot(Tl[0, :], yc, lw=2, color="k", ls="--", label=r"$T_{liquidus}$")
lTs1d, = axT1d.plot(Ts[0, :], yc, lw=2, color="b", ls="--", label=r"$T_{solidus}$")
axT1d.legend().set_draggable(True)

pcmU = axU.pcolormesh(xc, yb, qx.T, shading='flat')
cbU = fig.colorbar(pcmU, cax=caxU)
lU1d, = axU1d.plot(velocity(A, C), Y, lw=2, color="r")

pcmB = axB.pcolormesh(xb, yb, beta.T, shading='flat')
cbB = fig.colorbar(pcmB, cax=caxB)
lB1d, = axB1d.plot(beta[0, :], yc, lw=2, color="r", ls='-', label=r"$\beta$")
lBeq1d, = axB1d.plot(betaeq[0, :], yc, lw=2, color="k", ls='--', label=r"$\beta_{eq}$")
axB1d.legend().set_draggable(True)

vlW = axW.axvline(xc[0], lw=2, color="red")
vlT = axT.axvline(xc[0], lw=2, color="red")
vlU = axU.axvline(xc[0], lw=2, color="red")
vlV = axB.axvline(xc[0], lw=2, color="red")

ax_time_slider = fig.add_axes([0.15, 0.03, 0.6, 0.015])
time_slider = Slider(
    ax=ax_time_slider, 
    label='timesteps', 
    valmin=0, 
    valmax=len(data)-1,
    valinit=0, 
    valstep=1)

ax_xlim_slider = plt.axes([0.15, 0.06, 0.6, 0.015])
xlim_slider = Slider(
    ax=ax_xlim_slider, 
    label='xlim', 
    valmin=0, 
    valmax=len(xc)-1, 
    valinit=0, 
    valstep=1)

pcmLists = [pcmT, pcmU, pcmB]
def time_update(val):
    t = time_slider.val
    ix = xlim_slider.val
    fig.suptitle('time = {0:.2f} h'.format(data[t]["time"]/3600))
    hw = data[t]["halfwidth"]
    T = data[t]["Tmask"]
    qx = data[t]["qx"][1:-1, :]
    beta = data[t]["beta"]
    betaeq = data[t]["betaeq"]
    lW.set_ydata(hw)
    
    global pcmLists
    if pcmLists:
        pcmLists[0].remove()
        pcmLists[1].remove()
        pcmLists[2].remove()
        pcmLists = []
    pcmT = axT.pcolormesh(xb, yb, T.T, shading='flat', cmap='jet')
    fig.colorbar(pcmT, cax=caxT)
    
    pcmU = axU.pcolormesh(xc, yb, qx.T, shading='flat')
    fig.colorbar(pcmU, cax=caxU)
    
    pcmB = axB.pcolormesh(xb, yb, beta.T, shading='flat')
    fig.colorbar(pcmB, cax=caxB)
    pcmLists = [pcmT, pcmU, pcmB]
    
    lT1d.set_xdata(data[t]["temperature"][ix, :])
    lTl1d.set_xdata(data[t]["Tliquidus"][ix, :])
    lTs1d.set_xdata(data[t]["Tsolidus"][ix, :])
    
    A = data[t]["A"][ix+1, :]
    C = data[t]["C"][ix+1, :]
    V = velocity(A, C)
    lU1d.set_xdata(V)
    axU1d.set_xlim(0, max(V) + 1e-6)
    
    mu = data[t]["viscosity"][ix, :]
    lMu1d.set_xdata(mu)
    axMu1d.set_xlim(min(mu), max(mu))
    
    lB1d.set_xdata(beta[ix, :])
    lBeq1d.set_xdata(betaeq[ix, :])
    # axB1d.set_xlim(min(V), max(V) + 1e-10)

    
def xlim_update(val):
    t = time_slider.val
    ix = xlim_slider.val
    
    vlW.set_xdata([xc[ix]])
    vlT.set_xdata([xc[ix]])
    vlU.set_xdata([xc[ix]])
    vlV.set_xdata([xc[ix]])
    lT1d.set_xdata(data[t]["temperature"][ix, :])
    lTl1d.set_xdata(data[t]["Tliquidus"][ix, :])
    lTs1d.set_xdata(data[t]["Tsolidus"][ix, :])
    A = data[t]["A"][ix+1, :]
    C = data[t]["C"][ix+1, :]
    V = velocity(A, C)
    lU1d.set_xdata(V)
    axU1d.set_xlim(0, max(V) + 1e-6)
    
    mu = data[t]["viscosity"][ix, :]
    lMu1d.set_xdata(mu)
    axMu1d.set_xlim(min(mu), max(mu))
    
    beta = data[t]["beta"]
    betaeq = data[t]["betaeq"]
    lB1d.set_xdata(beta[ix, :])
    lBeq1d.set_xdata(betaeq[ix, :])
    # axB1d.set_xlim(min(V), max(V) + 1e-10)
    
# # Call update function when slider value is changed
time_slider.on_changed(time_update)
xlim_slider.on_changed(xlim_update)
plt.show()