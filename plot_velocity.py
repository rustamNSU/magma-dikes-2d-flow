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
simID = 140100
# simID = 102100
with open(sim_dir + "/simID{}.pkl".format(simID), 'rb') as inp:
    fc = pickle.load(inp)

fig, ax = plt.subplots(2, 1)
ax[0].set_xlabel("$x$, km")
ax[0].set_ylabel("$y$, m")
ax[0].set_ylim((0, 1))
ax[0].grid()

ax[1].set_xlabel("$x$, km")
ax[1].set_ylabel("$y$, m")
ax[1].set_ylim((0, 1))
ax[1].grid()

X = fc.xc2d / 1e3
Y = fc.yc2d
T = fc.T2d - 273.15
U = np.multiply(fc.q2d, fc.dp[:, np.newaxis])

ch = ax[0].contourf(X, Y, T)
clb = fig.colorbar(ch, ax=ax[0])
# vf = ax[0].quiver(X, Y, U, 0)

ch = ax[1].contourf(X, Y, -fc.q2d)
clb = fig.colorbar(ch, ax=ax[1])

fig.tight_layout()
fig.legend().set_draggable(True)
plt.show()