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
set_matplotlib_settings(DEFAULT_SIZE=12, LEGEND_SIZE=12)

timesteps = [20, 20, 36, 36]
simIDs = [110, 117, 160, 161]
simLegends = [
    r"quasi-2d",
    r"1d",
    r"quasi-2d, episodic",
    r"1d, episodic",
]
colors = cycle(['k', 'r', 'b', 'g'])
linestyles = cycle(['-', '--'])
markers = cycle(['o', 's'])  # circle, square, diamond

hour = 3600
time_shift = 100  # avoid log(0)
timeList, frontList, vtimeList, vList = [], [], [], []
simPaths = [sim_dir / f"simID{simID}" for simID in simIDs]
dikes = [DikeData(str(sim_path), step_rate=100) for sim_path in simPaths]
xmin = -30000
xmax = -30000
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

# --- Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6))
fig.subplots_adjust(hspace=0.4)

# Global X limits and ticks
xmin = min(min(t) for t in timeList)
xmax = max(max(t) for t in timeList)
xticks = [1, 10, 100]
xticklabels = ["1", "10", "100"]


ax1.set_xscale("log")
ax1.set_xlim([xmin, 2 * xmax])
ax1.set_xticks(xticks)
ax1.set_xticklabels(xticklabels)
ax1.grid(True, which='major', linestyle='--', linewidth=0.8, alpha=0.5)
ax1.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.3)
ax1.tick_params(labelbottom=True)
ax1.set_xlabel("Time (h)")
ax1.set_ylabel("Front (km)")


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
    
ax2.set_xlabel(r"$x$ (km)")
ax2.set_ylabel(r"Width (m)")
ax2.set_ylim([0, 4])
ax2.grid()

# Add legends
for ax in (ax1, ax2):
    ax.legend(loc="best").set_draggable(True)


ax_button = plt.axes([0.7, 0.05, 0.2, 0.075])  # Position of the button
def save_image(event):
    savepath = repository_dir / "images/article2024" / f"front_width_{'_'.join(map(str, simIDs))}"
    savepath.parent.mkdir(parents=True, exist_ok=True)
    ax_button.set_visible(False)
    fig.savefig(str(savepath) + ".png", dpi=600, bbox_inches='tight', pad_inches=0)
    fig.savefig(str(savepath) + ".pdf", bbox_inches='tight', pad_inches=0)
    print(f"Image saved at {savepath}.png")

# Create save button
button = Button(ax_button, 'Save image')
button.on_clicked(save_image)
 
plt.show()