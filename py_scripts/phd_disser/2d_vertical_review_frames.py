import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable

repository_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sim_dir = os.path.join(repository_dir, "simulations")
sys.path.append(repository_dir)

from py_scripts.utils import set_matplotlib_settings, create_layers_mask
from pysrc import DikeData

set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)


DEFAULT_SIM_ID = 110
DEFAULT_STEP_RATE = 10
DEFAULT_TIMESTEP = 100
DEFAULT_DPI = 600

BASE_CONFIG = {
    "xind": [50, 150, 220],
    "colors": ["k", "k", "k"],
    "linestyles": ["-", "--", "-."],
    "markers": ["None", "None", "None"],
    "xlim": (-30, -5),
    "Tlim": (400, 900),
    "blim": (0, 1.0),
    "Mlim": (3, 10),
}

SIM_CONFIG_OVERRIDES = {
    113: {
        "Mlim": (3, 14),
    },
    154: {
        "blim": (0, 0.4),
    },
}


def get_config(sim_id):
    return {**BASE_CONFIG, **SIM_CONFIG_OVERRIDES.get(sim_id, {})}


def build_geometry(data):
    xb = data.xb / 1e3
    yb = data.yb
    hw = data.halfwidth
    hwb = (0.5 * (hw[:-1] + hw[1:])).tolist()
    hwb = np.array([hw[0], *hwb, hw[-1]])
    hwb[hwb < 1e-4] = 1e-4

    x2db, _ = np.meshgrid(xb, yb, indexing="ij")
    y2db = np.array([a * yb for a in hwb])
    return xb, yb, hwb, x2db, y2db


def build_velocity_reference(data, config):
    xb, yb, hwb, _, y2db = build_geometry(data)
    xind = config["xind"]

    ny = 200
    Y = np.linspace(0, 1, ny)
    yid = create_layers_mask(yb, Y)

    def velocity(A, C):
        Ay = np.array([A[iy] for iy in yid])
        Cy = np.array([C[iy] for iy in yid])
        return Ay * Y * Y + Cy

    profiles = []
    for ix, xi in enumerate(xind):
        profile = velocity(data.A[xi, :], data.C[xi, :])
        profiles.append({
            "index": ix,
            "xi": xi,
            "x": xb[xi + 1],
            "y": Y * hwb[xi + 1],
            "velocity": profile,
            "width": hwb[xi + 1],
        })

    common_ylim = [
        min(np.min(profile["velocity"]) for profile in profiles),
        max(np.max(profile["velocity"]) for profile in profiles),
    ]
    if np.isclose(common_ylim[0], common_ylim[1]):
        delta = 1e-12 if np.isclose(common_ylim[0], 0.0) else abs(common_ylim[0]) * 0.05
        common_ylim = [common_ylim[0] - delta, common_ylim[1] + delta]
    else:
        padding = 0.05 * (common_ylim[1] - common_ylim[0])
        common_ylim = [common_ylim[0] - padding, common_ylim[1] + padding]
    field_xlim = (0.0, float(np.nanmax(y2db)))

    return {
        "profiles": profiles,
        "common_ylim": common_ylim,
        "field_xlim": field_xlim,
    }


def add_inset_colorbar_axis(ax):
    cax = inset_axes(
        ax,
        width="10%",
        height="40%",
        loc="lower right",
        bbox_to_anchor=(-0.3, 0.02, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    cax.yaxis.set_ticks_position("right")
    cax.yaxis.set_label_position("right")
    rect = Rectangle(
        (0.55, 0.01),
        0.42,
        0.42,
        transform=ax.transAxes,
        color="white",
        zorder=10,
        clip_on=False,
    )
    ax.add_patch(rect)
    return cax


def create_frame_figure(frame_data, reference, config):
    fig = plt.figure(figsize=(6.5, 5.5), layout="constrained")
    col1 = 7
    col2 = 0
    col3 = 10
    nx = len(config["xind"])
    nc = col1 + col2

    gs = fig.add_gridspec(1, 3 * nc + col3)
    gs1 = gs[0, 0:3 * nc].subgridspec(1, 3 * nc)
    gs2 = gs[0, 3 * nc:].subgridspec(nx, 1)

    axT = fig.add_subplot(gs1[0:col1])
    caxT = add_inset_colorbar_axis(axT)

    axB = fig.add_subplot(gs1[nc:nc + col1])
    caxB = add_inset_colorbar_axis(axB)

    axM = fig.add_subplot(gs1[2 * nc:2 * nc + col1])
    caxM = add_inset_colorbar_axis(axM)

    axU1ds = [fig.add_subplot(gs2[-ix - 1, 0]) for ix in range(nx)]
    axU1ds[-1].set_title(r"Vertical velocity")

    axT.set_title(r"Temperature (C$^\circ$)")
    axT.set_xlabel(r"Halfwidth (m)")
    axT.set_ylabel(r"Depth (km)")
    axT.set_ylim(config["xlim"])
    axT.set_xlim(reference["field_xlim"])
    axT.grid(which="major", linestyle="-", linewidth=0.75)
    axT.grid(which="minor", linestyle="-", linewidth=0.5)
    axT.minorticks_on()
    axT.set_axisbelow(True)
    axT.spines["top"].set_visible(False)
    axT.spines["right"].set_visible(False)

    axB.set_title(r"Crystal concentration")
    axB.set_xlabel(r"Halfwidth (m)")
    axB.sharey(axT)
    axB.set_ylim(config["xlim"])
    axB.set_xlim(reference["field_xlim"])
    axB.grid(which="major", linestyle="-", linewidth=0.75)
    axB.grid(which="minor", linestyle="-", linewidth=0.5)
    axB.minorticks_on()
    axB.set_axisbelow(True)
    axB.tick_params(axis="y", left=False, labelleft=False)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)

    axM.set_title(r"Viscosity, $\log_{10}$ (Pa$\cdot$s)")
    axM.set_xlabel(r"Halfwidth (m)")
    axM.sharey(axT)
    axM.set_ylim(config["xlim"])
    axM.set_xlim(reference["field_xlim"])
    axM.grid(which="major", linestyle="-", linewidth=0.75)
    axM.grid(which="minor", linestyle="-", linewidth=0.5)
    axM.minorticks_on()
    axM.set_axisbelow(True)
    axM.tick_params(axis="y", left=False, labelleft=False)
    axM.spines["top"].set_visible(False)
    axM.spines["right"].set_visible(False)

    for ax in axU1ds:
        ax.set_xlabel(r"$y$ (m)")
        ax.set_ylabel(r"Velocity (m/s)")
        ax.grid()

    xb, _, _, x2db, y2db = build_geometry(frame_data)
    T = np.where(frame_data.open_mask[:, None], frame_data.temperature, np.nan)
    beta = frame_data.beta
    viscosity_mask = frame_data.open_mask[:, None] & (frame_data.viscosity > 0.0)
    viscosity = np.full(frame_data.viscosity.shape, np.nan, dtype=float)
    viscosity[viscosity_mask] = np.log10(frame_data.viscosity[viscosity_mask])

    pcmT = axT.pcolormesh(
        y2db, x2db, T, shading="flat", cmap="jet",
        vmin=config["Tlim"][0], vmax=config["Tlim"][1],
    )
    cbT = fig.colorbar(pcmT, cax=caxT)
    cbT.ax.set_zorder(11)

    pcmB = axB.pcolormesh(
        y2db, x2db, beta, shading="flat", cmap="jet",
        vmin=config["blim"][0], vmax=config["blim"][1],
    )
    fig.colorbar(pcmB, cax=caxB)

    pcmM = axM.pcolormesh(
        y2db, x2db, viscosity, shading="flat", cmap="jet",
        vmin=config["Mlim"][0], vmax=config["Mlim"][1],
    )
    fig.colorbar(pcmM, cax=caxM)

    # Keep the right-column profiles fixed so only the 2D fields change across frames.
    for profile in reference["profiles"]:
        ix = profile["index"]
        axT.axhline(
            profile["x"], lw=2,
            ls=config["linestyles"][ix],
            color=config["colors"][ix],
            marker=config["markers"][ix],
        )
        axB.axhline(
            profile["x"], lw=2,
            ls=config["linestyles"][ix],
            color=config["colors"][ix],
            marker=config["markers"][ix],
        )
        axM.axhline(
            profile["x"], lw=2,
            ls=config["linestyles"][ix],
            color=config["colors"][ix],
            marker=config["markers"][ix],
        )
        axU1ds[ix].plot(
            profile["y"], profile["velocity"], lw=2,
            color=config["colors"][ix],
            ls=config["linestyles"][ix],
            marker=config["markers"][ix],
        )
        axU1ds[ix].set_xlim([0, profile["width"]])
        axU1ds[ix].set_ylim(reference["common_ylim"])

    return fig


def save_frames(sim_id, step_rate, timestep, dpi):
    config = get_config(sim_id)
    sim_path = os.path.join(sim_dir, f"simID{sim_id}")
    dike = DikeData(sim_path, step_rate=step_rate)

    max_timestep = len(dike.data) - 1
    last_timestep = min(timestep, max_timestep)
    if last_timestep < 1:
        raise ValueError("No timesteps available for frame generation.")
    if timestep > max_timestep:
        print(f"Requested timestep {timestep}, using available maximum {last_timestep}.")

    frame_dir = os.path.join(repository_dir, "images", "phd", f"2d_frames_simID{sim_id}")
    os.makedirs(frame_dir, exist_ok=True)

    reference = build_velocity_reference(dike.data[last_timestep], config)

    for frame_index in tqdm(range(1, last_timestep + 1), desc=f"simID{sim_id}"):
        fig = create_frame_figure(dike.data[frame_index], reference, config)
        frame_path = os.path.join(
            frame_dir,
            f"frame-{frame_index:03d}.png",
        )
        fig.savefig(frame_path, dpi=dpi)
        plt.close(fig)

    print(f"Saved {last_timestep} frames to {frame_dir}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export 2D review frames with fixed velocity profiles from the last timestep."
    )
    parser.add_argument("--simID", type=int, default=DEFAULT_SIM_ID)
    parser.add_argument("--step-rate", type=int, default=DEFAULT_STEP_RATE)
    parser.add_argument("--timestep", type=int, default=DEFAULT_TIMESTEP)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    return parser.parse_args()


def main():
    args = parse_args()
    save_frames(
        sim_id=args.simID,
        step_rate=args.step_rate,
        timestep=args.timestep,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
