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
from py_scripts.utils import set_matplotlib_settings, create_layers_mask

set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)


simIDs = [160, 161]
simLegends = [r"quasi-2d", r"1d"]
step_rate = 100
timestep = 35
profile_depths_km = [-17.0, -20.0, -25.0]
depth_linestyles = ["-", "--", "-."]
model_colors = ["k", "r"]
Tlim = (400, 900)
xlim = (-30, -5)

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


def profile_indices(frame, depths_km):
    xb = frame.xb / 1e3
    return [int(np.argmin(np.abs(xb[1:] - depth))) for depth in depths_km]


def velocity_profile(frame, xi, Y, yid):
    A = frame.A[xi, :]
    C = frame.C[xi, :]
    Ay = np.array([A[iy] for iy in yid])
    Cy = np.array([C[iy] for iy in yid])
    return Ay * Y * Y + Cy


profile_xinds = profile_indices(frames[1], profile_depths_km)
Y = np.linspace(0, 1, 200)
yid = create_layers_mask(frames[0].yb, Y)
hwb_list = [edge_halfwidths(frame) for frame in frames]
temperature_halfwidth_xlim = (0, max(np.nanmax(hwb) for hwb in hwb_list))

velocity_data = []
for xi in profile_xinds:
    velocity_data.append([
        velocity_profile(frame, xi, Y, yid)
        for frame in frames
    ])

velocity_min = min(np.nanmin(v) for pair in velocity_data for v in pair)
velocity_max = max(np.nanmax(v) for pair in velocity_data for v in pair)


fig = plt.figure(figsize=(8.0, 5.6))
gs = gridspec.GridSpec(
    3,
    4,
    figure=fig,
    width_ratios=[1.0, 1.0, 1.15, 0.08],
    left=0.08,
    right=0.92,
    bottom=0.13,
    top=0.88,
    wspace=0.32,
    hspace=0.42,
)

axT2d = fig.add_subplot(gs[:, 0])
axT1d = fig.add_subplot(gs[:, 1], sharey=axT2d)
axU = [fig.add_subplot(gs[i, 2]) for i in range(len(profile_xinds))]
caxT = fig.add_subplot(gs[:, 3])

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
    )
    ax.set_title(rf"\bf {legend}")
    ax.set_xlabel(r"Halfwidth (m)")
    ax.set_xlim(temperature_halfwidth_xlim)
    ax.set_ylim(xlim)
    ax.grid(which="major", linestyle="-", linewidth=0.2)
    ax.grid(which="minor", linestyle="-", linewidth=0.1)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axT2d.set_ylabel(r"Depth (km)")
axT1d.tick_params(axis="y", left=False, labelleft=False)
fig.colorbar(pcmT, cax=caxT, label=r"Temperature (C$^\circ$)")

for ax in temperature_axes:
    for xi, linestyle in zip(profile_xinds, depth_linestyles):
        ax.axhline(frames[0].xb[xi + 1] / 1e3, lw=1.8, ls=linestyle, color="k")

for i, (ax, xi, velocities) in enumerate(zip(axU, profile_xinds, velocity_data)):
    for imodel, frame in enumerate(frames):
        x = Y * hwb_list[imodel][xi + 1]
        ax.plot(
            x,
            velocities[imodel],
            lw=2,
            color=model_colors[imodel],
            label=simLegends[imodel],
        )

    depth_km = frames[0].xb[xi + 1] / 1e3
    max_halfwidth = max(hwb[xi + 1] for hwb in hwb_list)
    ax.set_title(rf"$x = {depth_km:.1f}$ km")
    ax.set_xlim([0, max_halfwidth])
    ax.set_ylim([velocity_min, velocity_max])
    ax.grid(which="major", linestyle="-", linewidth=0.2)
    ax.grid(which="minor", linestyle="-", linewidth=0.1)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.set_ylabel(r"Velocity (m/s)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if i < len(axU) - 1:
        ax.tick_params(axis="x", labelbottom=False)
    else:
        ax.set_xlabel(r"$y$ (m)")
        ax.legend(loc="best").set_draggable(True)

fig.suptitle(
    rf"Temperature and velocity, simID{simIDs[0]} vs simID{simIDs[1]}, "
    rf"$t = {frames[0].time / 3600:.1f}$ h"
)

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

plt.show()
