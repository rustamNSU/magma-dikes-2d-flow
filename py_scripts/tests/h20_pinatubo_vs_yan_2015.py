import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

K = 273.15
pmin = 0 # MPa
pmax = 1000 #MPa
Tmin = 400 # Celcius
Tmax = 1100 # Celcius
P1d = np.linspace(pmin, pmax, 100)
T1d = np.linspace(Tmin, Tmax, 100)
H2O_yan = lambda p, t: (354.94*np.sqrt(p) + 9.623*p - 1.5233*np.sqrt(p)*p) / (t + 273.15) + 0.0012439*np.sqrt(p)*p



pinatuboPaths = [
    "./py_scripts/tests/pdpath0.4.csv",
    "./py_scripts/tests/pdpath0.6.csv",
    "./py_scripts/tests/pdpath0.8.csv",
]


Plist = []
h2olist = []
pinatuboLegend = [
    ".pdpath0.4",
    ".pdpath0.6",
    ".pdpath0.8",
]
for path in pinatuboPaths:
    df = pd.read_csv(path)
    P = df["Pressure_bars"].to_numpy() * 1e5 / 1e6
    h2o = df["H2O_liq"].to_numpy()
    Plist.append(P)
    h2olist.append(h2o)




fig = plt.figure(figsize=(8, 6), constrained_layout=True)
gs = fig.add_gridspec(1, 1)
ax = fig.add_subplot(gs[0, 0])

Tlist = [400, 600, 800, 1000]
for t in Tlist:
    ax.plot(P1d, H2O_yan(P1d, t), lw=3, ls="--", label=fr"T = {t}$^\circ$C, Yan V. (2015)")


for i in range(len(pinatuboPaths)):
    ax.plot(Plist[i], h2olist[i], lw=3, ls="-", label=pinatuboLegend[i])

ax.set_xlabel(r"Pressure (MPa)")
ax.set_ylabel(r"Dissolved $H_2O$ (wt$\%$)")
ax.legend()
ax.grid()
plt.show()