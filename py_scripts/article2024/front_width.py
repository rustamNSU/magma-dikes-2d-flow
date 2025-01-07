import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button
from matplotlib.ticker import ScalarFormatter, FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_scripts.utils import set_matplotlib_settings, create_save_button
from py_scripts.DikeData import DikeData
from scipy.interpolate import splrep, splev, UnivariateSpline, PchipInterpolator
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=8)

# # h2o wt.%
# simIDs = [1, 2, 3]
# savename = repository_dir + "/images/article2024/FV_h2owt"
# simLegends = [r"$3.85$ wt.$\%$", r"$6.18$ wt.$\%$", r"$9.57$ wt.$\%$"]

# # dx: 100m vs 50m
# simIDs = [2, 4]
# savename = repository_dir + "/images/article2024/FV_dx_size"
# simLegends = [r"$\Delta x = 100$ m", r"$\Delta x = 50$ m"]

# # quasi 2d vs 1d
# simIDs = [2, 12]
# savename = repository_dir + "/images/article2024/FV_2dvs1d"
# simLegends = [r"quasi-2d model", r"1d model"]

# # 900 vs 850
# simIDs = [2, 32]
# savename = repository_dir + "/images/article2024/FV_chamber_temperature"
# simLegends = [r"$T_{ch}=900$ $^\circ$C", r"$T_{ch}=850$ $^\circ$C"]

# # tau: 1d vs 3d vs 0d
# simIDs = [52, 2, 42]
# savename = repository_dir + "/images/article2024/FV_cryst_relaxation"
# simLegends = [r"$\beta = \beta_{eq}$", r"$\tau = 1$ d", r"$\tau = 3$ d"]

# # KIc: 1 vs 100 vs 500 vs 1000
# simIDs = [2, 61, 62, 63]
# timesteps = [110, 110, 110, 110]
# savename = repository_dir + "/images/article2024/front_width_toughness"
# simLegends = [r"$K_{1c} = 1$ MPa$\cdot$m$^{1/2}$", r"$K_{1c} = 100$ MPa$\cdot$m$^{1/2}$", r"$K_{1c} = 500$ MPa$\cdot$m$^{1/2}$", r"$K_{1c} = 1000$ MPa$\cdot$m$^{1/2}$"]

# KIc: 1 vs 100 vs 500 vs 1000
simIDs = [70, 2, 71]
timesteps = [110, 110, 110]
savename = repository_dir + "/images/article2024/front_width_elasticity"
simLegends = [r"$E = 5$ GPa", r"$E = 20$ GPa", r"$E = 80$ GPa"]

# savename = repository_dir + "/images/front_width_" + + "_".join([str(sid) for sid in simIDs])
# simLegends = [str(simID) for simID in simIDs]
# simColors = ['r', 'b', 'k', 'r', 'b', 'k']
# simLinesteyle = ['-', '-', '-', '--', '--', '--']
simColors = ['r', 'k', 'b', 'darkgreen', 'b', 'k']
simLinesteyle = ['-', '-', '-', '-', '--', '--']

simPaths = [sim_dir + f"/simID{simID}" for simID in simIDs]
timeList = []
frontList = []
for simPath in simPaths:
    data = np.genfromtxt(simPath + "/front.txt", delimiter=";")
    time = data[:, 0]
    front = data[:, 1]
    timeList.append(time)
    frontList.append(front)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.set_xlabel(r"Time (hour)")
ax1.set_ylabel(r"Front (km)")
ax1.grid()
ax2.set_xlabel(r"Depth (m)")
ax2.set_ylabel(r"Width (m)")
ax2.grid()
for sid in range(len(simIDs)):
    time = timeList[sid]
    front = frontList[sid]
    ax1.plot(time / 3600.0, front / 1000.0, lw=3, label=simLegends[sid], ls=simLinesteyle[sid], color=simColors[sid])
    
    dike = DikeData(simPaths[sid], step_rate=10)
    data = dike.data[timesteps[sid]]
    xc = dike.data[timesteps[sid]]["xc"] / 1000.0
    w = dike.data[timesteps[sid]]["width"]
    mask = np.ma.where(w > 1e-4, True, False)
    mask = np.convolve(mask, np.array([True, True, True]), 'same')
    w[np.logical_not(mask)] = np.nan
    ax2.plot(xc, w, lw=3, label=simLegends[sid], ls=simLinesteyle[sid], color=simColors[sid])

ax2.set_ylim(bottom=0.0)
ax1.legend().set_draggable(True)
ax2.legend().set_draggable(True)


savepath = Path(savename)
os.makedirs(savepath.parent.absolute(), exist_ok=True)
bsave = create_save_button(savename)

plt.show()