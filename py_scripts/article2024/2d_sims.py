import sys, os
from pathlib import Path
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_scripts.utils import set_matplotlib_settings, create_layers_mask
from matplotlib.widgets import Button
from pysrc import *
from multiprocessing import Process
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)

simIDs = [110, 113]  # Two different simulations
xind = [45, 160, 220]
colors = ["k", "k", "k"]
linestyles = ["-", "--", "-."]
markers = ["None", "None", "None"]
xlim = (-30, -2)

def plot_sim(simID):
    # Load data
    sim_path = os.path.join(sim_dir, f"simID{simID}")
    dike = DikeData(sim_path, step_rate=10)
    timestep = 90
    data = dike.data[timestep]

    # --- Prepare figure
    fig = plt.figure(figsize=(9, 5.5), layout="constrained")
    fig.canvas.manager.set_window_title(f"simID{simID}")
    col1, col2, col3 = 7, 1, 10
    gs = fig.add_gridspec(1, 3*(col1+col2)+col3)
    gs1 = gs[0, 0:3*(col1+col2)].subgridspec(1, 3*(col1+col2))
    gs2 = gs[0, 3*(col1+col2):].subgridspec(3, 1)

    axs1, caxs1 = [], []
    for i in range(3):
        nc = i*(col1 + col2)
        axs1.append(fig.add_subplot(gs1[nc:nc+col1]))
        caxs1.append(fig.add_subplot(gs1[nc+col1:nc+col1+col2]))

    axT, axB, axV = axs1
    caxT, caxB, caxV = caxs1
    axU1ds = [fig.add_subplot(gs2[-ix-1, 0]) for ix in range(len(xind))]

    # --- Common plotting helper
    def draw_step(frame):
        # Grid, scaling, etc.
        fig.suptitle(f"Time = {frame.time / 3600:.2f} h", fontsize=12)
        for ax in axs1:
            ax.clear()
            ax.set_ylim(xlim)
            ax.grid()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        for axU1d in axU1ds:
            axU1d.clear()
            axU1d.set_xlabel(r"$y$ (m)")
            axU1d.set_ylabel(r"Velocity (m/s)")
            axU1d.grid()

        # Precompute mesh coords
        xc, xb = frame.xc / 1e3, frame.xb / 1e3
        yb = frame.yb
        hw = frame.halfwidth
        hwb = np.array([hw[0], *0.5*(hw[:-1] + hw[1:]), hw[-1]])
        hwb[hwb < 1e-4] = 1e-4
        x2db, _ = np.meshgrid(xb, yb, indexing='ij')
        y2db = np.array([a * yb for a in hwb])

        # Plot fields
        T = np.where(frame.open_mask[:, None], frame.temperature, np.nan)
        B = np.where(frame.open_mask[:, None], frame.beta, np.nan)
        visc = np.where(frame.open_mask[:, None], frame.viscosity, np.nan)
        visc = np.where(visc <= 1, 1, visc)  # avoid log10(0) or negatives
        V = np.log10(visc)

        pcmT = axT.pcolormesh(y2db, x2db, T, shading='flat', cmap='jet')
        fig.colorbar(pcmT, cax=caxT)
        axT.set_title(r"\bf Temperature (C$^\circ$)")
        axT.set_xlabel("Halfwidth (m)")
        axT.set_ylabel("Depth (km)")

        pcmB = axB.pcolormesh(y2db, x2db, B, shading='flat', cmap='jet')
        fig.colorbar(pcmB, cax=caxB)
        axB.set_title(r"\bf Crystal concentration")
        axB.set_xlabel("Halfwidth (m)")

        pcmV = axV.pcolormesh(y2db, x2db, V, shading='flat', cmap='jet')
        fig.colorbar(pcmV, cax=caxV)
        axV.set_title(r"\bf Viscosity ($\log_{10}$ Pa$\cdot$s)")
        axV.set_xlabel("Halfwidth (m)")
        axV.set_ylabel("Depth (km)")

        # Velocity profiles
        As = [frame.A[xi, :] for xi in xind]
        Cs = [frame.C[xi, :] for xi in xind]
        Ny = 200
        Y = np.linspace(0, 1, Ny)
        yid = create_layers_mask(yb, Y)

        def velocity(A, C):
            Ay = np.array([A[iy] for iy in yid])
            Cy = np.array([C[iy] for iy in yid])
            return Ay * Y**2 + Cy

        v_arrays = [velocity(As[i], Cs[i]) for i in range(len(xind))]
        v_ylim = [min(np.min(v) for v in v_arrays), max(np.max(v) for v in v_arrays)]

        for ix, v in enumerate(v_arrays):
            x = xb[xind[ix]+1]
            axT.axhline(x, color=colors[ix], ls=linestyles[ix])
            axB.axhline(x, color=colors[ix], ls=linestyles[ix])
            axU1ds[ix].plot(Y * hwb[xind[ix]+1], v, color=colors[ix], ls=linestyles[ix])
            axU1ds[ix].set_xlim([0, hwb[xind[ix]+1]])
            axU1ds[ix].set_ylim(v_ylim)

    # Initial draw
    draw_step(data)

    # Slider
    ax_slider = fig.add_axes([0.15, 0.02, 0.7, 0.03])
    slider = Slider(ax_slider, "Timestep", 0, len(dike.data) - 1, valinit=timestep, valstep=1)

    def update(val):
        step = int(slider.val)
        draw_step(dike.data[step])
        fig.canvas.draw_idle()

    slider.on_changed(update)

    # Save button
    ax_button = plt.axes([0.7, 0.15, 0.2, 0.075])
    def call(event):
        print("Saving figure for simID", simID)
        ax_button.set_visible(False)
        plt.savefig(repository_dir + f"/images/article2024/vertical_sim{simID}.png", dpi=600)
    button = Button(ax_button, 'Save image')
    button.on_clicked(call)
    plt.show()  # <- Show this figure immediately


if __name__ == "__main__":
    from multiprocessing import set_start_method
    set_start_method("spawn")  # safer on some OSes like macOS or Windows

    simIDs = [110, 113]
    processes = []
    for simID in simIDs:
        p = Process(target=plot_sim, args=(simID,))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()