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

simIDs = [1, 2, 3]
simIDs = [2, 4]
simIDs = [2, 12]
simIDs = [20]
simLegends = [str(simID) for simID in simIDs]
simColors = ['r', 'b', 'k', 'r', 'b', 'k']
simLinesteyle = ['-', '-', '-', '--', '--', '--']

simPaths = [sim_dir + f"/simID{simID}" for simID in simIDs]
timeList = []
frontList = []
vtimeList = []
vList = []
for simPath in simPaths:
    data = np.genfromtxt(simPath + "/front.txt", delimiter=";")
    time = data[:, 0]
    front = data[:, 1]
    timeList.append(time)
    frontList.append(front)
    
    data = np.genfromtxt(simPath + "/front_unique.txt", delimiter=";")
    time = data[:, 0]
    front = data[:, 1]
    v = np.diff(front) / np.diff(time)
    tv = time[:-1]
    
    vtimeList.append(tv)
    vList.append(v)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.set_xlabel(r"Time (hour)")
ax1.set_ylabel(r"Front (km)")
ax1.grid()
ax2.set_xlabel(r"Time")
ax2.set_ylabel(r"Ascend velocity (m$/$s)")
ax2.grid()
for sid in range(len(simIDs)):
    time = timeList[sid]
    front = frontList[sid]
    tv = vtimeList[sid]
    v = vList[sid]
    ax1.plot(time / 3600.0, front / 1000.0, lw=3, label=simLegends[sid], ls=simLinesteyle[sid], color=simColors[sid])
    ax2.loglog(tv, v, lw=3, label=simLegends[sid], ls=simLinesteyle[sid], color=simColors[sid])
    
ax1.legend().set_draggable(True)
ax2.legend().set_draggable(True)
ax2.set_xlim([100, None])
ax2.set_xticks([3600, 10*3600, 24*3600, 7*24*3600])
ax2.set_xticklabels(["1h", "10h", "1d", "1 w"])
formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
ax2.yaxis.set_major_formatter(formatter)
plt.show()