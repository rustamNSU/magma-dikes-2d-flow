import matplotlib
from matplotlib import rc
import numpy as np


def set_matplotlib_settings(DEFAULT_SIZE=10, LEGEND_SIZE=10, MAC_OS=False):
    rc('text', usetex=True)
    rc('font', size=DEFAULT_SIZE)          # controls default text sizes
    rc('axes', titlesize=DEFAULT_SIZE)     # fontsize of the axes title
    rc('axes', labelsize=DEFAULT_SIZE)    # fontsize of the x and y labels
    rc('xtick', labelsize=DEFAULT_SIZE)    # fontsize of the tick labels
    rc('ytick', labelsize=DEFAULT_SIZE)    # fontsize of the tick labels
    rc('legend', fontsize=LEGEND_SIZE)    # legend fontsize
    rc('figure', titlesize=DEFAULT_SIZE)  # fontsize of the figure title
    matplotlib.rcParams['text.latex.preamble']= \
        r"\usepackage[utf8]{inputenc} \usepackage[russian]{babel} \usepackage{amsmath} \usepackage{amssymb} \usepackage{bm}"
    matplotlib.rcParams['font.family'] = 'serif'
    if (MAC_OS):
        matplotlib.use('MacOSX')


def createXData(XC, DX):
    x = []
    for xc, dx in zip(XC, DX):
        x.append(xc - 0.5*dx)
        x.append(xc + 0.5*dx)
    return np.array(x)


def createYData(Y):
    yy = []
    for y in Y:
        yy.append(y)
        yy.append(y)
    return np.array(yy)


def create_layers_mask(yb, Y):
    Ny = len(Y)
    yid = np.zeros(Ny, dtype=int)
    ind = 0
    for iy in range(Ny):
        if Y[iy] > yb[ind+1] and Y[iy] < 1:
            ind = ind + 1
        yid[iy] = ind
    return yid