import sys
from pathlib import Path
repository_dir = Path.cwd()
sim_dir = repository_dir / "simulations"
sys.path.append(str(repository_dir))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from itertools import cycle
from pysrc import *
from py_scripts.utils import set_matplotlib_settings
set_matplotlib_settings(DEFAULT_SIZE=12, LEGEND_SIZE=12)

simIDs = [160, 161]
timesteps = [360, 360]
simLegends = [
    r"quasi-2d",
    r"1d",
]
colors = cycle(['k', 'r', 'b', 'b'])
linestyles = cycle(['-', '--', '-.'])
markers = cycle(['o', 's', 'D'])  # circle, square, diamond

hour = 3600
day = hour * 24
time_shift = 100  # avoid log(0)
timeList, frontList, vtimeList, vList = [], [], [], []
wtimeList, wList = [], []
dikes = []
for simID in simIDs:
    simPath = sim_dir / f"simID{simID}"
    time, front = np.genfromtxt(simPath / "front.txt", delimiter=";").T
    time += time_shift
    timeList.append(time / day)
    frontList.append(front)

    time_u, front_u = np.genfromtxt(simPath / "front_unique.txt", delimiter=";").T
    time_u += time_shift
    v = np.diff(front_u)[1:] / np.diff(time_u)[1:]
    vtime = time_u[2:]
    vtimeList.append(vtime / day)
    vList.append(v)
    
    dike = DikeData(str(simPath), step_rate=10)
    dikes.append(dike)
    wmax = []
    time = []
    for frame in dike.data:
        wmax.append(max(frame.width))
        time.append((frame.time + time_shift) / day)
    wtimeList.append(time)
    wList.append(wmax)


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))
fig.subplots_adjust(hspace=0.4)

xmin = min(min(t) for t in timeList)
xmax = max(max(t) for t in timeList)
# xticks = [1, 10, 100]
# xticklabels = ["1", "10", "100"]

# ax1.set_xscale("log")
ax1.set_xlim([0, xmax])
ax1.grid()
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("Front (km)")
ax1.set_ylim([-30.5, -2])


ax2.set_xlabel(r"$x$ (km)")
ax2.set_ylabel(r"Width (m)")
ax2.set_ylim([0, 4])
ax2.grid()

# Plot curves
for i, dike in enumerate(dikes):
    color = next(colors)
    linestyle = next(linestyles)
    marker = next(markers)

    ax1.plot(timeList[i], frontList[i] / 1000,
             lw=3, ls=linestyle, color=color, marker=marker,
             label=simLegends[i], markevery=0.1)
    
    frame = dike.data[timesteps[i]]
    x = frame.xc
    w = frame.width
    mask = np.ma.where(w > 1e-4, True, False)
    mask = np.convolve(mask, np.array([True, True, True]), 'same')
    w[np.logical_not(mask)] = np.nan
    ax2.plot(x / 1000.0, w, lw=3, ls=linestyle, color=color, marker=marker,
             label=simLegends[i], markevery=0.1)

# Add legends
for ax in (ax1, ax2):
    ax.legend(loc="best").set_draggable(True)

props = dict(ha='center', va='top', fontsize=12)
# ax1.text(0.5, 1.11, r"\textbf{(a)}", transform=ax1.transAxes, **props)
# ax2.text(0.5, 1.11, r"\textbf{(c)}", transform=ax3.transAxes, **props)

ax_button = plt.axes([0.7, 0.05, 0.2, 0.075])  # Position of the button
def save_image(event):
    savepath = repository_dir / "images/article2024" / f"FW_{'_'.join(map(str, simIDs))}"
    savepath.parent.mkdir(parents=True, exist_ok=True)
    ax_button.set_visible(False)
    fig.savefig(str(savepath) + ".png", dpi=600, bbox_inches='tight', pad_inches=0)
    fig.savefig(str(savepath) + ".pdf", bbox_inches='tight', pad_inches=0)
    print(f"Image saved at {savepath}.png")

# Create save button
button = Button(ax_button, 'Save image')
button.on_clicked(save_image)

plt.show()