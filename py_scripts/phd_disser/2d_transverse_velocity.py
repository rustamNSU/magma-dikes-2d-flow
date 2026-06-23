import sys
from pathlib import Path

repository_dir = Path.cwd()
sim_dir = repository_dir / "simulations"
sys.path.append(str(repository_dir))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import ScalarFormatter
from matplotlib.widgets import Button
from pysrc import DikeData
from py_scripts.utils import set_matplotlib_settings

set_matplotlib_settings(DEFAULT_SIZE=12, LEGEND_SIZE=12)


simID = 110
step_rate = 10
timestep = 20
depth_xlim = (-30, -5)
grid_major_color = "0.35"
grid_minor_color = "0.55"

dike = DikeData(str(sim_dir / f"simID{simID}"), step_rate=step_rate)
frame = dike.data[timestep]


def edge_halfwidths(data):
    hw = data.halfwidth
    hwb = (0.5 * (hw[:-1] + hw[1:])).tolist()
    hwb = np.array([hw[0], *hwb, hw[-1]])
    hwb[hwb < 1e-4] = 1e-4
    return hwb


def transverse_velocity_grid(data):
    xb = data.xb / 1e3
    yc = data.yc
    dx = data.xb[1] - data.xb[0]
    hwb = edge_halfwidths(data)

    x2db, _ = np.meshgrid(xb, yc, indexing="ij")
    y2db = np.array([hw * yc for hw in hwb])

    # qy is staggered: centered in x and located on y-boundaries.
    # Internal y-boundary values are drawn in cells bounded by yc.
    vy = data.qy[:, 1:-1] / dx
    vy = np.where(data.open_mask[:, None], vy, np.nan)
    return y2db, x2db, vy


y2db, x2db, vy = transverse_velocity_grid(frame)
vabs = np.nanmax(np.abs(vy))
if not np.isfinite(vabs) or vabs == 0:
    vabs = 1.0

fig = plt.figure(figsize=(4.8, 7.2))
gs = gridspec.GridSpec(
    1,
    2,
    figure=fig,
    width_ratios=[1.0, 0.07],
    left=0.13,
    right=0.88,
    bottom=0.10,
    top=0.90,
    wspace=0.22,
)

axV = fig.add_subplot(gs[0, 0])
caxV = fig.add_subplot(gs[0, 1])

pcmV = axV.pcolormesh(
    y2db,
    x2db,
    vy,
    shading="flat",
    cmap="seismic",
    norm=TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs),
    zorder=2,
)

axV.set_title(r"\bf Transverse velocity")
axV.set_xlabel(r"Halfwidth (m)")
axV.set_ylabel(r"Depth (km)")
axV.set_xlim(0, np.nanmax(y2db))
axV.set_ylim(depth_xlim)
axV.minorticks_on()
axV.set_axisbelow(True)
axV.grid(
    which="major",
    linestyle="-",
    linewidth=0.6,
    color=grid_major_color,
    alpha=1.0,
    zorder=0,
)
axV.grid(
    which="minor",
    linestyle="-",
    linewidth=0.3,
    color=grid_minor_color,
    alpha=1.0,
    zorder=0,
)
axV.spines["top"].set_visible(False)
axV.spines["right"].set_visible(False)

cbV = fig.colorbar(pcmV, cax=caxV, label=r"$v_y$ (m/s)")
cbV.formatter = ScalarFormatter(useMathText=True)
cbV.formatter.set_powerlimits((-2, 2))
cbV.update_ticks()

fig.suptitle(rf"simID{simID}, $t = {frame.time / 3600:.1f}$ h")

ax_button = fig.add_axes([0.62, 0.025, 0.25, 0.055])
ax_button.set_in_layout(False)


def save_image(event):
    savepath = (
        repository_dir
        / "images/phd"
        / f"2d_transverse_velocity_simID{simID}_timestep{timestep * step_rate}"
    )
    savepath.parent.mkdir(parents=True, exist_ok=True)
    ax_button.set_visible(False)
    fig.savefig(str(savepath) + ".png", dpi=600, bbox_inches="tight", pad_inches=0)
    fig.savefig(str(savepath) + ".pdf", bbox_inches="tight", pad_inches=0)
    print(f"Image saved at {savepath}.png")


button = Button(ax_button, "Save image")
button.on_clicked(save_image)

print(rf"$t = {frame.time / 3600:.1f}$ h")

plt.show()
