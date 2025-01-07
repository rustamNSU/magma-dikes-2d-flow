import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.widgets import Button
import matplotlib
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
        r"\usepackage[utf8]{inputenc} \usepackage{amsmath} \usepackage{amssymb} \usepackage{bm}"
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


def save_pdf(pdf_final_filename):
    pdf_ext = '.pdf'
    pdf_final_filename_with_ext = pdf_final_filename + pdf_ext
    plt.savefig(pdf_final_filename_with_ext, bbox_inches='tight', pad_inches=0)


def save_png(png_final_filename):
    png_ext = '.png'
    png_final_filename_with_ext = png_final_filename + png_ext
    plt.savefig(png_final_filename_with_ext, dpi=300, bbox_inches='tight', pad_inches=0, format="pdf")
    

def save_eps(eps_final_filename):
    eps_ext = '.eps'
    eps_final_filename_with_ext = eps_final_filename + eps_ext
    plt.savefig(eps_final_filename_with_ext)   


def create_save_button(file_path, axs=None):
    ax_button = plt.axes([0.7, 0.15, 0.2, 0.075])

    def call(event):
        print("button clicked")
        ax_button.set_visible(False)
        if axs is not None:
            [ax.set_visible(False) for ax in axs]
        save_pdf(file_path)
        save_png(file_path)
        # save_eps(file_path)
    button = Button(ax_button, 'Save image')
    button.on_clicked(call)
    return button