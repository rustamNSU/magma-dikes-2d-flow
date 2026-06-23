import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button
from matplotlib.ticker import ScalarFormatter, FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_scripts.utils import set_matplotlib_settings, create_layers_mask
from py_scripts.DikeData import DikeData
from scipy.interpolate import splrep, splev, UnivariateSpline, PchipInterpolator
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)

simIDs = [20, 23]
timesteps = [147, 400]
# simIDs = [20, 24]
# timesteps = [147, 202]
simPaths = [sim_dir + f"/simID{simID}" for simID in simIDs]
dikes = [DikeData(simPath, step_rate=10).data for simPath in simPaths]
# simLegends = ["this study", "Roper, Lister (2007)"]
# datas = np.genfromtxt(simPath + "/front.txt", delimiter=";")
# time = data[:, 0]
# front = data[:, 1]

# data = np.genfromtxt(simPath + "/front_unique.txt", delimiter=";")
# time = data[:-1:10, 0]
# front = data[:-1:10, 1]
# dt = np.diff(time)
# v = np.diff(front) / np.diff(time)
# tv = time[:-1]




fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.set_xlabel(r"$x$ (km)")
ax1.set_ylabel(r"Width (m)")
ax1.grid()
ax2.set_xlabel(r"Time (s)")
ax2.set_ylabel(r"Ascend velocity (m$/$s)")
ax2.grid()


for i in range(len(simIDs)):
    timestep = timesteps[i]
    xc = dikes[i][timestep]["xc"] / 1000.0
    w = dikes[i][timestep]["width"]
    w[w <= 1e-3] = np.nan
    ax1.plot(xc, w, lw=3, label=simIDs[i], ls='-', color='b')
# ax2.plot(tv, v, lw=3, label="this study", ls='-', color='b')

# Roper, Lister (2007)
data_list = [np.genfromtxt(repository_dir + path, delimiter=";") for path in [
    "/data/roper_lister/roper_lister_K1.csv",
    "/data/roper_lister/roper_lister_K8.csv"
]]
xScales = [1825.74, 1825.74/4]
for i in range(2):
    data = data_list[i]
    xScale = xScales[i]
    X = -data[:, 0] * xScale / 1e3 - 15.2
    Y = 2*data[:, 1]
    X = np.append(X, [-40])
    Y = np.append(Y, [2])
    ax1.plot(X, Y, lw=3, label="Roper, Lister (2007)", ls='--', color='r')
# ax2.axhline(y=10.0, lw=3, label="Roper, Lister (2007)", ls='--', color='r')

# ax1.legend().set_draggable(True)
# ax2.legend().set_draggable(True)
# ax1.set_xlim([-30, -10])
# ax1.set_ylim([0, 4])
# ax2.set_xlim([50, None])
# ax2.set_ylim([8, 15])
# formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
# ax2.yaxis.set_major_formatter(formatter)
plt.show()