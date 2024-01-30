import pickle
import numpy as np
import itertools
from scipy.interpolate import RectBivariateSpline
from srcpython.core.fixed_channel import FixedChannel
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

sim_dir = "simulations"
simID = 140100
# simID = 21050
with open(sim_dir + "/simID{}.pkl".format(simID), 'rb') as inp:
    fc = pickle.load(inp)

fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
X = fc.xc2d
Y = fc.yc2d
T = fc.T2d
surf = ax.plot_surface(X, Y, T, cmap=cm.coolwarm, linewidth=0, antialiased=False)
fig.tight_layout()
plt.show()