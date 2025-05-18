import sys
from pathlib import Path
repository_dir = Path.cwd()
sys.path.append(str(repository_dir))

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from py_scripts.utils import set_matplotlib_settings, create_save_button
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)

# Files to plot
filepaths = [repository_dir / "data/pinatubo" / file for file in [
    "dissolved_04.json",
    "dissolved_06.json",
    "dissolved_10.json",
]]

# Read isobars from CSV
isobar_data = pd.read_csv(repository_dir / "data/pinatubo/pinatubo-isobars.csv")
grouped_isobars = {int(p): df for p, df in isobar_data.groupby("Pressure")}

# Label, marker, color for each JSON
styles = [
    (r"$c_{H_2O} = 3.85\ \mathrm{wt\%}$", "o", "blue"),
    (r"$c_{H_2O} = 6.18\ \mathrm{wt\%}$", "^", "black"),
    (r"$c_{H_2O} = 9.57\ \mathrm{wt\%}$", "s", "red")
]

# Start figure
fig, ax = plt.subplots(figsize=(7, 4))

# Plot each dissolved volatile curve
for i, filepath in enumerate(filepaths):
    with open(filepath) as f:
        data = json.load(f)

    wth2o = np.array(data["wth2o"]) * 100
    wtco2 = np.array(data["wtco2"]) * 100

    label, marker, color = styles[i]
    ax.plot(wth2o, wtco2, linestyle='-', marker=marker, label=label,
            color=color, markersize=6, markerfacecolor='white',
            markeredgewidth=1.2, linewidth=3, markevery=5)

# Plot real isobars from CSV (in bar → MPa)
for p_bar, df in grouped_isobars.items():
    h2o = df["H2O_liq"]
    co2 = df["CO2_liq"]
    ax.plot(h2o, co2, "k--", linewidth=2)
    idx = co2.idxmax()
    ax.text(h2o.loc[idx], co2.loc[idx], f"{int(p_bar / 10)}",
            ha="center", va="bottom",
            bbox=dict(facecolor='wheat', edgecolor='grey', pad=2))

# Final touches
ax.set_xlabel(r"$H_2O$ content, wt\%", labelpad=10)
ax.set_ylabel(r"$CO_2$ content, wt\%", labelpad=10)
ax.set_xlim(0, 20)
ax.set_ylim(0, 0.6)
ax.legend().set_draggable(True)
ax.grid()
fig.tight_layout()

savepath = repository_dir / "images/article2024" / "mixed_saturation"
savepath.parent.mkdir(parents=True, exist_ok=True)
bsave = create_save_button(str(savepath))
plt.show()