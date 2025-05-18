import sys
from pathlib import Path
repository_dir = Path.cwd()
sim_dir = repository_dir / "simulations"
sys.path.append(str(repository_dir))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from itertools import cycle
from pysrc import *
from py_scripts.utils import set_matplotlib_settings
set_matplotlib_settings(DEFAULT_SIZE=14, LEGEND_SIZE=14)

simIDs = [110, 130, 131]
simLegends = [
    r"$\Delta x = 100$ m, $N_y = 30$",
    r"$\Delta x = 100$ m, $N_y = 50$",
    r"$\Delta x = 50$ m, $N_y = 30$",
]
colors = cycle(['k', 'r', 'g', 'b'])
linestyles = cycle(['-', '--', '-.'])
markers = cycle(['o', 's', 'D'])  # circle, square, diamond

time_shift = 100  # avoid log(0)
timeList, frontList, vtimeList, vList = [], [], [], []

for simID in simIDs:
    simPath = sim_dir / f"simID{simID}"

    time, front = np.genfromtxt(simPath / "front.txt", delimiter=";").T
    time += time_shift
    timeList.append(time)
    frontList.append(front)

    time_u, front_u = np.genfromtxt(simPath / "front_unique.txt", delimiter=";").T
    time_u += time_shift
    v = np.diff(front_u)[1:] / np.diff(time_u)[1:]
    vtime = time_u[2:]
    vtimeList.append(vtime)
    vList.append(v)

# --- Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
fig.subplots_adjust(hspace=0.4)

# Global X limits and ticks
xmin = min(min(t) for t in timeList)
xmax = max(max(t) for t in timeList)
xticks = [3600, 36000, 360000]
xticklabels = ["1h", "10h", "100h"]

for ax in (ax1, ax2):
    ax.set_xscale("log")
    ax.set_xlim([xmin, 2 * xmax])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.grid(True, which='major', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.3)
    
ax1.tick_params(labelbottom=True)
ax1.set_xlabel("Time (h)")
ax1.set_ylabel("Front (km)")
ax2.set_xlabel("Time (h)")
ax2.set_ylabel("Velocity (m/s)")
ax2.set_yscale("log")
ax2.set_yticks([1, 0.1, 0.01])
ax2.set_yticklabels(["1", "0.1", "0.01"])

# Plot curves
for i in range(len(simIDs)):
    color = next(colors)
    linestyle = next(linestyles)
    marker = next(markers)

    ax1.plot(timeList[i], frontList[i] / 1000,
             lw=3, ls=linestyle, color=color, marker=marker,
             label=simLegends[i], markevery=0.1)

    ax2.plot(vtimeList[i], vList[i],
             lw=3, ls=linestyle, color=color, marker=marker,
             label=simLegends[i], markevery=0.1)

# Add legends
for ax in (ax1, ax2):
    ax.legend(fontsize=12, loc="best").set_draggable(True)

props = dict(ha='center', va='top', fontsize=14)
ax1.text(0.5, 1.11, r"\textbf{(a)}", transform=ax1.transAxes, **props)
ax2.text(0.5, 1.11, r"\textbf{(b)}", transform=ax2.transAxes, **props)

ax_button = plt.axes([0.7, 0.05, 0.2, 0.075])  # Position of the button
def save_image(event):
    savepath = repository_dir / "images/article2024" / f"FV_{'_'.join(map(str, simIDs))}"
    savepath.parent.mkdir(parents=True, exist_ok=True)
    ax_button.set_visible(False)
    fig.savefig(str(savepath) + ".png", dpi=600, bbox_inches='tight', pad_inches=0)
    fig.savefig(str(savepath) + ".pdf", bbox_inches='tight', pad_inches=0)
    print(f"Image saved at {savepath}.png")

# Create save button
button = Button(ax_button, 'Save image')
button.on_clicked(save_image)

plt.show()