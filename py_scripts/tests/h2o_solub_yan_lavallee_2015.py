import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

K = 273.15
pmin = 0 # MPa
pmax = 1000 #MPa
Tmin = 400 # Celcius
Tmax = 1100 # Celcius
P1d = np.linspace(pmin, pmax, 100)
T1d = np.linspace(Tmin, Tmax, 100)
P, T = np.meshgrid(P1d, T1d)

h2o = lambda p, t: (354.94*np.sqrt(p) + 9.623*p - 1.5233*np.sqrt(p)*p) / (t + 273.15) + 0.0012439*np.sqrt(p)*p
Z = h2o(P, T)

fig = plt.figure(figsize=(12, 12), constrained_layout=True)
nax = 10
ncb = 1
gs = fig.add_gridspec(2, nax+ncb+nax)
ax11 = fig.add_subplot(gs[0, 0:nax])
ax12 = fig.add_subplot(gs[0, nax+ncb:], projection='3d')
ax21 = fig.add_subplot(gs[1, 0:nax])
ax22 = fig.add_subplot(gs[1, nax+ncb:])
ax11cb = fig.add_subplot(gs[0, nax])

Zcut = Z
Zcut[Zcut > 7] = np.nan
Zcut[Zcut < 0] = np.nan
pc11 = ax11.contourf(P, T, Z, cmap='viridis', levels=100)
cb11 = fig.colorbar(pc11, cax=ax11cb)
ax11.set_xlabel(r"Pressure (MPa)")
ax11.set_ylabel(r"Temperature ($^\circ$C)")

ax12.plot_surface(P, T, Z, cmap='viridis')
ax12.set_xlabel(r"Pressure (MPa)")
ax12.set_ylabel(r"Temperature ($^\circ$C)")
ax12.set_zlabel(r"Dissolved $H_2O$ (wt$\%$)")

colors = ["r", "b", "k"]
Tlist = [600, 800, 1000]
for color, t in zip(colors, Tlist):
    ax21.plot(P1d, h2o(P1d, t), lw=3, color=color, label=fr"T = {t}$^\circ$C")
    ax11.axhline(t, lw=3, color=color)
ax21.set_xlabel(r"Pressure (MPa)")
ax21.set_ylabel(r"Dissolved $H_2O$ (wt$\%$)")
ax21.legend()

colors = ["r", "b", "k"]
Plist = [400, 600, 800]
for color, p in zip(colors, Plist):
    ax22.plot(T1d, h2o(p, T1d), lw=3, ls="--", color=color, label=fr"p = {p} MPa")
    ax11.axvline(p, lw=3, ls="--", color=color)
ax22.set_xlabel(r"Temperature ($^\circ$C)")
ax22.set_ylabel(r"Dissolved $H_2O$ (wt$\%$)")
ax22.legend()
plt.show()