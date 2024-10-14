import sys, os
repository_dir = os.path.abspath(os.getcwd())
sim_dir = repository_dir + "/simulations"
sys.path.append(repository_dir)


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
from py_scripts.utils import set_matplotlib_settings, create_layers_mask
from py_scripts.DikeData import DikeData
from scipy.interpolate import splrep, splev, UnivariateSpline
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10)

# simIDs = [3, 4, 5]
# simLegends = [r"$\Delta x = 100$", r"$\Delta x = 50$", r"$\Delta x = 25$"]
simIDs = [10, 11, 12]
simIDs = [12, 13]
simLegends = [str(simID) for simID in simIDs]
simPaths = [sim_dir + f"/simID{simID}" for simID in simIDs]
dikes = [DikeData(sim_path, step_rate=1) for sim_path in simPaths]
timeList = []
frontList = []
for dikedata in dikes:
    time, front = map(list, zip(*[(data["time"], data["tip_front"]) for data in dikedata.data]))
    timeList.append(np.array(time))
    frontList.append(np.array(front))


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.set_xlabel(r"Time (hour)")
ax1.set_ylabel(r"Front (km)")
ax1.grid()
ax2.set_xlabel(r"Time (hour)")
ax2.set_ylabel(r"Ascend velocity (m$/$s)")
ax2.grid()
for sid in range(len(simIDs)):
    time = timeList[sid]
    front = frontList[sid]
    v = np.zeros(front.shape)
    N = min(20, len(time)-1)
    for i in range(len(time)):
        if i < N:
            v[i] = (front[i+1] - front[i]) / (time[i+1] - time[i])
        else:
            df = front[i] - front[i-N]
            dt = time[i] - time[i-N]
            v[i] = df / dt
    
    ax1.plot(time / 3600.0, front / 1000.0, lw=3, label=simLegends[sid])
    ax2.plot(time / 3600.0, v, lw=3, label=simLegends[sid])
    
    # func = splrep(time, front, k=3, s=0)
    # t = np.linspace(time[0], time[-1], 1000)
    # f = splev(t, func)
    # # front = gaussian_filter1d(f, 11)
    # front2 = savgol_filter(f, window_length=21, polyorder=3)
    # func = splrep(t, front2, k=3, s=0)
    # ax1.plot(time / 3600.0, front / 1000.0, lw=3, label=simLegends[i])
    # ax2.plot(t / 3600.0, splev(t, func, der=1), lw=3, label=simLegends[i])
    # z = np.polyfit(time, front, deg=11)
    # p = np.poly1d(z)
    # dp = np.polyder(p)
    # ax1.plot(time / 3600.0, p(time) / 1000.0, lw=3, label=simLegends[i])
    # ax2.plot(time / 3600.0, dp(time), lw=3, label=simLegends[i])
    # f = savgol_filter(front, window_length=11, polyorder=2)
    # df = savgol_filter(front, window_length=11, polyorder=2, deriv=1)
    # ax1.plot(time / 3600.0, f / 1000.0, lw=3, label=fr"simID = {simIDs[i]}")
    # ax2.plot(time / 3600.0, df, lw=3, label=fr"simID = {simIDs[i]}")
    
ax1.legend().set_draggable(True)
ax2.legend().set_draggable(True)
plt.show()






# from matplotlib import pyplot as plt
# import numpy as np
# from scipy.interpolate import splrep, splev

# x = np.arange(0,2,0.008)
# data = np.polynomial.polynomial.polyval(x,[0,2,1,-2,-3,2.6,-0.4])
# noise = np.random.normal(0,0.1,250)
# noisy_data = data + noise

# f = splrep(x,noisy_data,k=5,s=3)
# #plt.plot(x, data, label="raw data")
# #plt.plot(x, noise, label="noise")
# plt.plot(x, noisy_data, label="noisy data")
# plt.plot(x, splev(x,f), label="fitted")
# plt.plot(x, splev(x,f,der=1)/10, label="1st derivative")
# #plt.plot(x, splev(x,f,der=2)/100, label="2nd derivative")
# plt.hlines(0,0,2)
# plt.legend(loc=0)
# plt.show()