import sys
from pathlib import Path

repository_dir = Path.cwd()
sim_dir = repository_dir / "simulations"
sys.path.append(str(repository_dir))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button
from pysrc import DikeData
from py_scripts.utils import set_matplotlib_settings

set_matplotlib_settings(DEFAULT_SIZE=14, LEGEND_SIZE=14)


simIDs = [160, 161]
simLegends = [r"Квазидвумерная модель", r"Одномерная модель"]
step_rate = 100
timestep = 35
Tlim = (400, 900)
xlim = (-30, -10)
temperature_halfwidth_xlim = (0, 2.1)
grid_major_color = "0.35"
grid_minor_color = "0.5"

dikes = [
    DikeData(str(sim_dir / f"simID{simID}"), step_rate=step_rate)
    for simID in simIDs
]
frames = [dike.data[timestep] for dike in dikes]


def edge_halfwidths(frame):
    hw = frame.halfwidth
    hwb = (0.5 * (hw[:-1] + hw[1:])).tolist()
    hwb = np.array([hw[0], *hwb, hw[-1]])
    hwb[hwb < 1e-4] = 1e-4
    return hwb


def temperature_grid(frame):
    xb = frame.xb / 1e3
    yb = frame.yb
    hwb = edge_halfwidths(frame)
    x2db, _ = np.meshgrid(xb, yb, indexing="ij")
    y2db = np.array([hw * yb for hw in hwb])
    temperature = np.where(frame.open_mask[:, None], frame.temperature, np.nan)
    return y2db, x2db, temperature


hwb_list = [edge_halfwidths(frame) for frame in frames]


fig = plt.figure(figsize=(7, 9))
gs = gridspec.GridSpec(
    1,
    3,
    figure=fig,
    width_ratios=[1.0, 1.0, 0.08],
    left=0.08,
    right=0.90,
    bottom=0.13,
    top=0.88,
    wspace=0.30,
)

axT2d = fig.add_subplot(gs[0, 0])
axT1d = fig.add_subplot(gs[0, 1], sharey=axT2d)
caxT = fig.add_subplot(gs[0, 2])

temperature_axes = [axT2d, axT1d]
pcmT = None
for ax, frame, legend in zip(temperature_axes, frames, simLegends):
    y2db, x2db, temperature = temperature_grid(frame)
    pcmT = ax.pcolormesh(
        y2db,
        x2db,
        temperature,
        shading="flat",
        cmap="jet",
        vmin=Tlim[0],
        vmax=Tlim[1],
        zorder=2,
    )
    ax.set_title(rf"\bf {legend}")
    ax.set_xlabel(r"Halfwidth (m)")
    ax.set_xlim(temperature_halfwidth_xlim)
    ax.set_ylim(xlim)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.grid(
        which="major",
        linestyle="-",
        linewidth=0.8,
        color=grid_major_color,
        alpha=1.0,
        zorder=0,
    )
    ax.grid(
        which="minor",
        linestyle="-",
        linewidth=0.3,
        color=grid_minor_color,
        alpha=1.0,
        zorder=0,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axT2d.set_ylabel(r"Depth (km)")
axT1d.tick_params(axis="y", left=False, labelleft=False)
fig.colorbar(pcmT, cax=caxT, label=r"Temperature (C$^\circ$)")

# fig.suptitle(
#     rf"Temperature, simID{simIDs[0]} vs simID{simIDs[1]}, "
#     rf"$t = {frames[0].time / 3600:.1f}$ h"
# )

ax_button = fig.add_axes([0.72, 0.03, 0.2, 0.065])
ax_button.set_in_layout(False)


def save_image(event):
    savepath = (
        repository_dir
        / "images/phd"
        / f"2d_2dvs1d_{'_'.join(map(str, simIDs))}_timestep{timestep * step_rate}"
    )
    savepath.parent.mkdir(parents=True, exist_ok=True)
    ax_button.set_visible(False)
    fig.savefig(str(savepath) + ".png", dpi=600, bbox_inches="tight", pad_inches=0)
    fig.savefig(str(savepath) + ".pdf", bbox_inches="tight", pad_inches=0)
    print(f"Image saved at {savepath}.png")


button = Button(ax_button, "Save image")
button.on_clicked(save_image)

print(rf"$t = {frames[0].time / 3600:.1f}$ h")

plt.show()
