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
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10, MAC_OS=True)

savename = repository_dir + "/images/article2024/roper_lister_verification"
simID = 20
simPath = sim_dir + f"/simID{simID}"
dike_data = DikeData(simPath, step_rate=10).data
simLegends = ["this study", "Roper, Lister (2007)"]
data = np.genfromtxt(simPath + "/front.txt", delimiter=";")
time = data[:, 0]
front = data[:, 1]

data = np.genfromtxt(simPath + "/front_unique.txt", delimiter=";")
time = data[:-1:10, 0]
front = data[:-1:10, 1]
dt = np.diff(time)
v = np.diff(front) / np.diff(time)
tv = time[:-1]




fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.set_xlabel(r"$x$ (km)")
ax1.set_ylabel(r"Width (m)")
ax1.grid()
ax2.set_xlabel(r"Time (s)")
ax2.set_ylabel(r"Ascend velocity (m$/$s)")
ax2.grid()


for i in [45, 90, 180]:
    label_text = "this study" if i == 45 else '_nolegend_'
    w = np.copy(dike_data[i]["width"])
    w[w <= 1e-3] = np.nan
    ax1.plot(dike_data[i]["xc"] / 1000.0, w, lw=3, label=label_text, ls='-', color='b')
    
    xc = dike_data[i]["xc"] / 1000.0
    w = dike_data[i]["width"]
    xc[w < 1e-3] = -1e6
    X = max(xc)
    ax1.text(X, 1.0, r"{0:.2f}h".format(dike_data[i]["time"]/3600), horizontalalignment='center', verticalalignment='center',
                    fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=1, pad=0.2),
                    zorder=30)
ax2.plot(tv, v, lw=3, label="this study", ls='-', color='b')

# Roper, Lister (2007)
data = np.genfromtxt(repository_dir + "/data/roper_lister/roper_lister_K1.csv", delimiter=";")
xScale = 1825.74
X = -data[:, 0] * xScale / 1e3 - 11.97
Y = 2*data[:, 1]
X = np.append(X, [-40])
Y = np.append(Y, [2])
ax1.plot(X, Y, lw=3, label="Roper, Lister (2007)", ls='--', color='r')
ax2.axhline(y=10.0, lw=3, label="Roper, Lister (2007)", ls='--', color='r')

ax1.legend().set_draggable(True)
ax2.legend().set_draggable(True)
ax1.set_xlim([-30, -10])
ax1.set_ylim([0, 4])
ax2.set_xlim([50, None])
ax2.set_ylim([8, 15])
formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
ax2.yaxis.set_major_formatter(formatter)
savepath = Path(savename)
os.makedirs(savepath.parent.absolute(), exist_ok=True)
bsave = create_save_button(savename)

plt.show()