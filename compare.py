import pickle
import numpy as np
from src.core.fixed_channel import FixedChannel

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
DEFAULT_SIZE = 10
LEGEND_SIZE = 10
SMALL_SIZE = 8
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

sim_dir = "simulations"
# simIDs = [11050, 12050, 14050, 110100, 120100, 140100]
simIDs = [105100, 110100, 120100, 140100]
data = []
for simID in simIDs:
    with open(sim_dir + "/simID{}.pkl".format(simID), 'rb') as inp:
        data.append(pickle.load(inp))

fig, ax = plt.subplots(3, 1)
ax[0].set_xlabel("$x$, km")
ax[0].set_ylabel("$p$, Pa")
ax[0].grid()

ax[1].set_xlabel("$x$, km")
ax[1].set_ylabel("$p$, Pa")
ax[1].grid()
ax[1].set_xlim((0, 1))
ax[1].set_ylim((3.75e5, 4.2e5))

for fc, simID in zip(data, simIDs):
    xc = fc.xc / 1e3
    p = fc.p
    ax[0].plot(xc, p, label=str(simID), linewidth = 2)
    ax[1].plot(xc, p, linewidth = 2)
    
ax[2].set_xlabel("$x$, km")
ax[2].set_ylabel("$T$, $^\circ$C")
ax[2].grid()

for fc, simID in zip(data, simIDs):
    xc = fc.xc / 1e3
    T = np.mean(fc.T2d, axis=1)
    ax[2].plot(xc, T, linewidth = 2)

fig.tight_layout()
fig.legend().set_draggable(True)
plt.show()