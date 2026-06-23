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

simIDs = [110, 113, 154]
simLegends = [
    r"$\tau_c(T) = \tau_0 \exp\left(\frac{E_a}{R T}\right)$",
    r"$\tau_c = 20$ s",
    r"$\tau_c = 1$ week",
]
colors = cycle(['k', 'r', 'b', 'b'])
linestyles = cycle(['-', '--', '-.'])
markers = cycle(['o', 's', 'D'])  # circle, square, diamond

hour = 3600
time_shift = 100  # avoid log(0)
timeList, frontList, vtimeList, vList = [], [], [], []
wtimeList, wList = [], []
for simID in simIDs:
    simPath = sim_dir / f"simID{simID}"
    time, front = np.genfromtxt(simPath / "front.txt", delimiter=";").T
    time += time_shift
    timeList.append(time / hour)
    frontList.append(front)

    time_u, front_u = np.genfromtxt(simPath / "front_unique.txt", delimiter=";").T
    time_u += time_shift
    v = np.diff(front_u)[1:] / np.diff(time_u)[1:]
    vtime = time_u[2:]
    vtimeList.append(vtime / hour)
    vList.append(v)
    
    dike = DikeData(str(simPath), step_rate=10)
    wmax = []
    time = []
    for frame in dike.data:
        wmax.append(max(frame.width))
        time.append((frame.time + time_shift) / hour)
    wtimeList.append(time)
    wList.append(wmax)


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 8), sharex=True)
fig.subplots_adjust(hspace=0.4)

axin1 = inset_axes(
    ax1,
    width="25%",       # ширина относительно ax1
    height="45%",      # высота
    bbox_to_anchor=(0.33, -0.1, 1, 1),  # (x, y, width, height) в координатах ax1
    bbox_transform=ax1.transAxes,
    loc='center'
)
axin1.set_xscale("log")
axin1.set_xlim(60, 100)
axin1.set_ylim(-8, -4)
axin1.tick_params(
    axis='both',
    which='both',
    bottom=False,
    left=True,
    labelbottom=False,
    labelleft=True
)
axin1.grid(True, which='major', linestyle='--', linewidth=0.8, alpha=0.5)
axin1.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.3)
axin1.minorticks_on()
axin1.tick_params(labelsize=8)
mark_inset(ax1, axin1, loc1=2, loc2=1, fc="none", ec="0.5")


xmin = min(min(t) for t in timeList)
xmax = max(max(t) for t in timeList)
xticks = [1, 10, 100]
xticklabels = ["1", "10", "100"]

for ax in (ax1, ax2, ax3):
    ax.set_xscale("log")
    ax.set_xlim([xmin, 2 * xmax])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.grid(True, which='major', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.3)
    
ax1.tick_params(labelbottom=True)
ax1.set_xlabel("Time (h)")
ax1.set_ylabel("Front (km)")
ax1.set_ylim([-30.5, -2])

ax2.set_xlabel("Time (h)")
ax2.set_ylabel("Velocity (m/s)")
ax2.set_yscale("log")
ax2.set_yticks([1, 0.1, 0.01])
ax2.set_yticklabels(["1", "0.1", "0.01"])

ax3.set_xlabel("Time (h)")
ax3.set_ylabel("Maximum width (m)")
ax3.set_ylim([0, 2])

# Plot curves
for i in range(len(simIDs)):
    color = next(colors)
    linestyle = next(linestyles)
    marker = next(markers)

    ax1.plot(timeList[i], frontList[i] / 1000,
             lw=3, ls=linestyle, color=color, marker=marker,
             label=simLegends[i], markevery=0.1)
    
    axin1.plot(timeList[i], frontList[i] / 1000,
            lw=2, ls=linestyle, color=color, marker=marker,
            markevery=0.3)

    ax2.plot(vtimeList[i], vList[i],
             lw=3, ls=linestyle, color=color, marker=marker,
             label=simLegends[i], markevery=0.1)
    
    ax3.plot(wtimeList[i], wList[i], lw=3, ls=linestyle, color=color, marker=marker,
             label=simLegends[i], markevery=0.1)

# Add legends
for ax in (ax1,):
    ax.legend(loc="best").set_draggable(True)

props = dict(ha='center', va='top', fontsize=12)
ax1.text(0.5, 1.11, r"\textbf{(a)}", transform=ax1.transAxes, **props)
ax2.text(0.5, 1.11, r"\textbf{(b)}", transform=ax2.transAxes, **props)
ax3.text(0.5, 1.11, r"\textbf{(c)}", transform=ax3.transAxes, **props)

ax_button = plt.axes([0.7, 0.05, 0.2, 0.075])  # Position of the button
def save_image(event):
    savepath = repository_dir / "images/article2024" / f"FVW_{'_'.join(map(str, simIDs))}"
    savepath.parent.mkdir(parents=True, exist_ok=True)
    ax_button.set_visible(False)
    fig.savefig(str(savepath) + ".png", dpi=600, bbox_inches='tight', pad_inches=0)
    fig.savefig(str(savepath) + ".pdf", bbox_inches='tight', pad_inches=0)
    print(f"Image saved at {savepath}.png")

# Create save button
button = Button(ax_button, 'Save image')
button.on_clicked(save_image)

plt.show()