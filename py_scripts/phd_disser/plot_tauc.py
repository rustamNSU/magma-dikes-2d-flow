import sys
from pathlib import Path
repository_dir = Path.cwd()
sim_dir = repository_dir / "simulations"
sys.path.append(str(repository_dir))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from py_scripts.utils import set_matplotlib_settings, create_save_button
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)

# Parameters
tau0 = 1e-6           # s
E = 210_000.0         # J/mol
R = 8.314             # J/mol/K

# Temperature range in Celsius
T_C = np.linspace(500, 1000, 500)   # [°C]
T_K = T_C + 273.15                  # Convert to Kelvin

# Relaxation time τ(T)
tau = tau0 * np.exp(E / (R * T_K))

# Plotting
plt.figure(figsize=(7, 4))
plt.plot(T_C, tau, lw=2)
plt.yscale('log')

plt.xlabel(r'Temperature $^\circ$C')
plt.ylabel(r'Characteristic time $\tau_c$ (s)')

plt.grid(True, which='major', linestyle='--', alpha=0.4)
plt.tight_layout()

savepath = repository_dir / "images/article2024" / "characteristic_crystal_growth_time"
savepath.parent.mkdir(parents=True, exist_ok=True)
bsave = create_save_button(str(savepath))
plt.show()
