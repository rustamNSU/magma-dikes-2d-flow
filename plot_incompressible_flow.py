import pickle
import numpy as np
import itertools
from scipy.interpolate import RectBivariateSpline
from src.core.fixed_channel import FixedChannel
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

sim_dir = "simulations"
simIDs = [10530, 11030, 12030, 12050, 120100, 140100]
simIDs = [12050]
data = []
for simID in simIDs:
    with open(sim_dir + "/simID{}.pkl".format(simID), 'rb') as inp:
        data.append(pickle.load(inp))

# fig, ax = plt.subplots(3,2)
# ax = list(itertools.chain.from_iterable(ax))
fig, ax = plt.subplots()
plh = ax.contourf(data[0].xc2d, data[0].yc2d, data[0].T2d)
clb = fig.colorbar(plh, ax=ax)
# for fc in data:
#     X = fc.xc2d
#     Y = fc.yc2d
#     T = fc.T2d
#     mask = X > 8000
#     xx, yy, zz = X[mask], Y[mask], T[mask]
#     xx = np.reshape(xx, (-1, X.shape[1]))
#     yy = np.reshape(yy, (-1, X.shape[1]))
#     zz = np.reshape(zz, (-1, X.shape[1]))
#     surf = ax.plot_surface(xx, yy, zz, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# for i in range(1, len(simIDs)):
#     fc2 = data[i]
#     fc1 = data[i-1]
    
#     interp_fc2 = RectBivariateSpline(fc2.xc, fc2.yc, fc2.T2d)
#     X = fc1.xc2d
#     Y = fc1.yc2d
#     T1 = fc1.T2d
#     T2 = interp_fc2(fc1.xc, fc1.yc)
#     surf = ax.plot_surface(X, Y, T2-T1, label="({})-({})".format(i, i-1))
#     surf._facecolors2d = surf._facecolor3d
#     surf._edgecolors2d = surf._edgecolor3d

# for i in range(1, len(simIDs)):
#     fc2 = data[i]
#     fc1 = data[i-1]
    
#     interp_fc2 = RectBivariateSpline(fc2.xc, fc2.yc, fc2.T2d)
#     X = fc1.xc2d
#     Y = fc1.yc2d
#     T1 = fc1.T2d
#     T2 = interp_fc2(fc1.xc, fc1.yc)
#     plh = ax[i-1].contourf(X, Y, T2-T1)
#     clb = fig.colorbar(plh, ax=ax[i-1])
#     ax[i-1].set_title("error {}--{}".format(simIDs[i-1], simIDs[i]))

fig.tight_layout()
plt.show()