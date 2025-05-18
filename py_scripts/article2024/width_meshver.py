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

timesteps = [20, 20, 20]
simIDs = [110, 130, 131]
simLegends = [
    r"$\Delta x = 100$ m, $N_y = 30$",
    r"$\Delta x = 100$ m, $N_y = 50$",
    r"$\Delta x = 50$ m, $N_y = 30$",
]
colors = cycle(['k', 'r', 'g', 'b'])
linestyles = cycle(['-', '--', '-.'])
markers = cycle(['o', 's', 'D'])  # circle, square, diamond


fig, ax = plt.subplots(1, 1, figsize=(7, 3))
fig.tight_layout()
simPaths = [sim_dir / f"simID{simID}" for simID in simIDs]
dikes = [DikeData(str(sim_path), step_rate=100) for sim_path in simPaths]
xmin = -30000
xmax = -30000
for i, dike in enumerate(dikes):
    color = next(colors)
    linestyle = next(linestyles)
    marker = next(markers)
    frame = dike.data[timesteps[i]]
    x = frame.xc
    w = frame.width
    
    mask = np.ma.where(w > 1e-4, True, False)
    mask = np.convolve(mask, np.array([True, True, True]), 'same')
    w[np.logical_not(mask)] = np.nan
    ax.plot(x / 1000.0, w, lw=3, ls=linestyle, color=color, marker=marker,
             label=simLegends[i], markevery=0.1)
    
ax.set_xlabel(r"$x$ (km)")
ax.set_ylabel(r"Width (m)")
ax.set_ylim([0, 1.5])
ax.grid()
ax.legend(loc="best").set_draggable(True)

ax_button = plt.axes([0.7, 0.05, 0.2, 0.075])  # Position of the button
def save_image(event):
    savepath = repository_dir / "images/article2024" / f"WIDTH_{'_'.join(map(str, simIDs))}"
    savepath.parent.mkdir(parents=True, exist_ok=True)
    ax_button.set_visible(False)
    fig.savefig(str(savepath) + ".png", dpi=600, bbox_inches='tight', pad_inches=0)
    fig.savefig(str(savepath) + ".pdf", bbox_inches='tight', pad_inches=0)
    print(f"Image saved at {savepath}.png")

# Create save button
button = Button(ax_button, 'Save image')
button.on_clicked(save_image)
plt.show()