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
set_matplotlib_settings(DEFAULT_SIZE=14, LEGEND_SIZE=12)

simIDs = [110, 113, 154]
simLegends = [
    r"$\tau_c(T) = \tau_0 \exp\left(\frac{E_a}{R T}\right)$",
    r"$\tau_c = 20$ s",
    r"$\tau_c = 1$ week",
]
colors = cycle(['k', 'r', 'b', 'b'])
linestyles = cycle(['-', '--', '-.'])
markers = cycle(['o', 's', 'D'])  # circle, square, diamond


fig, ax = plt.subplots(1, 1, figsize=(7, 3))
fig.tight_layout()
simPaths = [sim_dir / f"simID{simID}" for simID in simIDs]
dikes = [DikeData(str(sim_path), step_rate=10) for sim_path in simPaths]
hour = 3600
time_shift = 100
timeList = []
for i, dike in enumerate(dikes):
    color = next(colors)
    linestyle = next(linestyles)
    marker = next(markers)
    wmax = []
    time = []
    for frame in dike.data:
        wmax.append(max(frame.width))
        time.append((frame.time + time_shift) / hour)
    timeList.append(time)
    
    ax.plot(time, wmax, lw=3, ls=linestyle, color=color, marker=marker,
             label=simLegends[i], markevery=0.1)

xmin = min(min(t) for t in timeList)
xmax = max(max(t) for t in timeList)
xticks = [1, 10, 100]
xticklabels = ["1", "10", "100"]
ax.set_xscale("log")
ax.set_xlim([xmin, 2 * xmax])
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
ax.grid(True, which='major', linestyle='--', linewidth=0.8, alpha=0.5)
ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.3)  
ax.set_xlabel("Time (h)")
ax.set_ylabel("Maximum width (m)")
ax.set_ylim([0, 2])
ax.legend(loc="best").set_draggable(True)

ax_button = plt.axes([0.7, 0.05, 0.2, 0.075])  # Position of the button
def save_image(event):
    savepath = repository_dir / "images/article2024" / f"MAXWIDTH_{'_'.join(map(str, simIDs))}"
    savepath.parent.mkdir(parents=True, exist_ok=True)
    ax_button.set_visible(False)
    fig.savefig(str(savepath) + ".png", dpi=600, bbox_inches='tight', pad_inches=0)
    fig.savefig(str(savepath) + ".pdf", bbox_inches='tight', pad_inches=0)
    print(f"Image saved at {savepath}.png")

# Create save button
button = Button(ax_button, 'Save image')
button.on_clicked(save_image)
plt.show()