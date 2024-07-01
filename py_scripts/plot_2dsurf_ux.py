import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_scripts.utils import set_matplotlib_settings
from py_scripts.DikeData import DikeData
set_matplotlib_settings()


simID = 3
sim_path = sim_dir + f"/simID{simID}"
wlim = (0, 3)
plim = (0, 700)
Tmin = 600
Tmax = 920
xlim = (-30000, -10000.0)


dike = DikeData(sim_path, step_rate=10)
data = dike.data
fig, (ax2, ax) = plt.subplots(1, 2, figsize=(16, 4),  width_ratios=[1, 4])
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
fig.subplots_adjust(bottom=0.3)
ax.set_xlabel(r"Depth")
ax.set_ylabel(r"$\xi = y/h(x,t)$")

x = data[0]["xc"]
y = data[0]["yb"]
xb = data[0]["xb"]
ux = data[0]["ux"][1:-1, :].T

# Velocity field
A = data[0]["A"][0, :]
C = data[0]["C"][0, :]
Ny = 200
Y = np.linspace(0, 1, Ny)
yb = data[0]["yb"]
yid = np.zeros(Ny, dtype=int)
ind = 0
for iy in range(Ny):
    if Y[iy] > yb[ind+1] and Y[iy] < 1:
        ind = ind + 1
    yid[iy] = ind

def velocity(A, C):
    Ay = np.array([A[iy] for iy in yid])
    Cy = np.array([C[iy] for iy in yid])
    return Ay*Y*Y + Cy

# pax = ax.pcolormesh(x, y, qx, vmin=0, vmax=0.1, shading='flat')
pax = ax.pcolormesh(x, y, ux, shading='flat')
cb = fig.colorbar(pax, cax=cax)

lu, = ax2.plot(velocity(A, C), Y, lw=2)
ax.axvline(xb[0], lw=2, color="red")
ax2.set_ylim(0, 1)
ax2.set_xlim(0, 1)

ax_time_slider = fig.add_axes([0.25, 0.05, 0.6, 0.05])
time_slider = Slider(
    ax=ax_time_slider, 
    label='timesteps', 
    valmin=0, 
    valmax=len(data)-1,
    valinit=0, 
    valstep=1)

ax_xlim_slider = plt.axes([0.25, 0.15, 0.6, 0.05])
xlim_slider = Slider(
    ax=ax_xlim_slider, 
    label='xlim', 
    valmin=0, 
    valmax=len(xb)-1, 
    valinit=0, 
    valstep=1)


def time_update(val):
    t = time_slider.val
    ix = xlim_slider.val
    fig.suptitle('time = {0:.2f} h'.format(data[t]["time"]/3600))
    ux = data[t]["ux"][1:-1, :].T
    ax.clear()
    pax = ax.pcolormesh(x, y, ux, shading='flat')
    fig.colorbar(pax, cax=cax)
    ax.axvline(xb[ix], lw=2, color="red")
    
    A = data[t]["A"][ix, :]
    C = data[t]["C"][ix, :]
    V = velocity(A, C)
    lu.set_xdata(V)
    ax2.set_xlim(0, max(V) + 1e-3)

    
def xlim_update(val):
    t = time_slider.val
    ix = xlim_slider.val
    ux = data[t]["ux"][1:-1, :].T
    ax.clear()
    pax = ax.pcolormesh(x, y, ux, shading='flat')
    fig.colorbar(pax, cax=cax)
    ax.axvline(xb[ix], lw=2, color="red")
    
    A = data[t]["A"][ix, :]
    C = data[t]["C"][ix, :]
    V = velocity(A, C)
    lu.set_xdata(V)
    ax2.set_xlim(0, max(V) + 1e-3)
    
# Call update function when slider value is changed
time_slider.on_changed(time_update)
xlim_slider.on_changed(xlim_update)

plt.show()