import json
import h5py
import os
import numpy as np
import matplotlib
# matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import interpolate
from matplotlib.widgets import Slider, Button

def createXData(XC, DX):
    x = []
    for xc, dx in zip(XC, DX):
        x.append(xc - 0.5*dx)
        x.append(xc + 0.5*dx)
    return np.array(x)


def createYData(Y):
    yy = []
    for y in Y:
        yy.append(y)
        yy.append(y)
    return np.array(yy)


DEFAULT_SIZE = 10
LEGEND_SIZE = 10
SMALL_SIZE = 10
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
rc('text', usetex=True)
rc('font', size=DEFAULT_SIZE)          # controls default text sizes
rc('axes', titlesize=DEFAULT_SIZE)     # fontsize of the axes title
rc('axes', labelsize=DEFAULT_SIZE)    # fontsize of the x and y labels
rc('xtick', labelsize=DEFAULT_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=DEFAULT_SIZE)    # fontsize of the tick labels
rc('legend', fontsize=LEGEND_SIZE)    # legend fontsize
rc('figure', titlesize=DEFAULT_SIZE)  # fontsize of the figure title
matplotlib.rcParams['text.latex.preamble']= \
    r"\usepackage[utf8]{inputenc} \usepackage[russian]{babel} \usepackage{amsmath} \usepackage{amssymb} \usepackage{bm}"
matplotlib.rcParams['font.family'] = 'serif'


sim_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + "/../simulations")
simIDs = [1]
legends = [str(simID) for simID in simIDs]
timesteps = list(range(0, 2001, 10))
wlim = (0, 3)
plim = (0, 700)
Tlim = (400, 920)
xlim = (-30000, -10000.0)


xc = []
dx = []
width = []
pressure = []
time = []
Tmiddle = []
Twall = []
ny = []
for simID in simIDs:
    xx = []
    dxx = []
    ww = []
    pp = []
    tt = []
    TTm = []
    TTb = []
    for t in timesteps:
        filepath = sim_dir + "/simID{}/data/data_{}.h5".format(simID, t)
        data = h5py.File(filepath, 'r')
        xx.append(np.array(data["mesh"]["xc"]))
        dxx.append(np.array(data["mesh"]["xr"]) - np.array(data["mesh"]["xl"]))
        ww.append(np.array(data["width"]))
        pp.append(1e-6*np.array(data["pressure"]))
        tt.append(float(np.array(data["time"])))
        nyy = int(np.array(data["ny"]))
        TT = np.array(data["temperature"])
        TTm.append(TT[:, 0])
        TTb.append(np.array(data["Twall"]))
    xc.append(xx)
    dx.append(dxx)
    width.append(ww)
    pressure.append(pp)
    time.append(tt)
    Tmiddle.append(TTm)
    Twall.append(TTb)
    ny.append(nyy)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6))
fig.subplots_adjust(bottom=0.2)
ax1.set_xlabel("x")
ax1.set_ylabel("width [m]")
ax2.set_xlabel("x")
ax2.set_ylabel("pressure [MPa]")
ax3.set_xlabel("x")
ax3.set_ylabel("temperature [C]")

lws = []
lps = []
lTms = []
lTbs = []
colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628']
for i in range(len(simIDs)):
    simID = simIDs[i]
    X = createXData(xc[i][0], dx[i][0])
    W = createYData(width[i][0])
    P = createYData(pressure[i][0])
    Tm = createYData(Tmiddle[i][0])
    Tb = createYData(Twall[i][0])
    Ny = ny[i]
    
    lw, = ax1.plot(X, W, linewidth=2, label=legends[i], linestyle='-', color=colors[i])
    lp, = ax2.plot(X, P, linewidth=2, label=legends[i], linestyle='-', color=colors[i])
    lTm, = ax3.plot(X, Tm, linewidth=2, label=legends[i] + ", mid.", linestyle='-', color=colors[i])
    lTb, = ax3.plot(X, Tb, linewidth=2, label=legends[i] + ", bound.", linestyle='--', color=colors[i])
    lws.append(lw)
    lps.append(lp)
    lTms.append(lTm)
    lTbs.append(lTb)

ax1.set_xlim(xlim)
ax1.set_ylim(wlim)
ax1.grid()
ax1.legend().set_draggable(True)

ax2.set_xlim(xlim)
ax2.set_ylim(plim)
ax2.grid()
ax2.legend().set_draggable(True)

ax3.set_xlim(xlim)
ax3.set_ylim(Tlim)
ax3.grid()
ax3.legend().set_draggable(True)

ax_iter_slider = fig.add_axes([0.25, 0.05, 0.65, 0.03])
iter_slider = Slider(
    ax=ax_iter_slider, 
    label='timesteps', 
    valmin=timesteps[0], 
    valmax=timesteps[-1], 
    valinit=timesteps[0], 
    valstep=timesteps)

# xlim_slider = Slider(
#     ax=ax_xlim_slider, 
#     label='xlim', 
#     valmin=1, 
#     valmax=max_X, 
#     valinit=xlim, 
#     valstep=xstep)


def iter_update(val):
    t = iter_slider.val
    t_index = timesteps.index(t)
    fig.suptitle(f'time = {tt[t_index]} s.')
    for i in range(len(simIDs)):
        X = xc[i][t_index]
        W = createYData(width[i][t_index])
        P = createYData(pressure[i][t_index])
        Tm = createYData(Tmiddle[i][t_index])
        Tb = createYData(Twall[i][t_index])
        
        lws[i].set_ydata(W)
        lps[i].set_ydata(P)
        lTms[i].set_ydata(Tm)
        lTbs[i].set_ydata(Tb)
    
# def xlim_update(val):
#     x_slider = xlim_slider.val
#     ax1.set_xlim(-x_slider, x_slider)
#     ax2.set_xlim(-x_slider, x_slider)
#     ax3.set_xlim(-x_slider, x_slider)
    
# Call update function when slider value is changed
iter_slider.on_changed(iter_update)
# xlim_slider.on_changed(xlim_update)
# fig.tight_layout()
# fig.subplots_adjust(left=0.1, bottom=0.25, right=0.9, top=0.98, wspace=0.25, hspace=0.3)
plt.show()