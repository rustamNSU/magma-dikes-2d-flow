import json
import h5py
import numpy as np
import fnmatch, os
import re
from tqdm import tqdm


def hdf5_read_single(h5file, key):
    if key in h5file.keys():
        data_h5 = h5file[key]
        if (data_h5.size == 1):
            return float(np.array(data_h5))
    return None

def hdf5_read_int(h5file, key):
    if key in h5file.keys():
        data_h5 = h5file[key]
        if (data_h5.size == 1):
            return int(np.array(data_h5))
    return None


def hdf5_read_array(h5file, key):
    if key in h5file.keys():
        data_h5 = h5file[key]
        return np.array(data_h5)
    return None


def define_index(str):
    a = re.split(r'_|\.', str)
    return int(a[1])


class DikeData:
    def __init__(self, sim_path: str, step_rate=1, zero_closed_elements=True) -> None:
        self.sim_path = sim_path
        self.input_path = sim_path + "/input.json"
        self.input = json.load(open(self.input_path))
        self.data_dir = sim_path + "/data"
        self.filepaths = fnmatch.filter(os.listdir(self.data_dir), 'data_*.h5')
        self.filepaths.sort(key=define_index)
        self.ntimesteps = len(self.filepaths)
        self.step_rate = step_rate
        self.zero_closed_elements = zero_closed_elements
        self._parse_input()
        self._preload_data()
        
        
    def _parse_input(self) -> None:
        self.nx = int(self.input["meshProperties"]["n"])
        self.xmin = float(self.input["meshProperties"]["xmin"])
        self.xmax = float(self.input["meshProperties"]["xmax"])
    
    
    def _preload_data(self):
        self.timesteps = range(0, self.ntimesteps, self.step_rate)
        self.data = []
        self.time = []
        min_width = self.input["algorithmProperties"]["minWidth"]
        for itime in tqdm(self.timesteps):
            file = self.data_dir + "/" + self.filepaths[itime]
            data = h5py.File(file, 'r')
            self.keys = list(data.keys())
            xc = hdf5_read_array(data, "/mesh/xc")
            xl = hdf5_read_array(data, "/mesh/xl")
            xr = hdf5_read_array(data, "/mesh/xr")
            xb = np.append(xl, xr[-1])
            yc = hdf5_read_array(data, "/yc")
            yb = hdf5_read_array(data, "/yb")
            halfwidth = hdf5_read_array(data, "/halfwidth")
            width = hdf5_read_array(data, "/width")
            pressure = hdf5_read_array(data, "/pressure")
            overpressure = hdf5_read_array(data, "/overpressure")
            density = hdf5_read_array(data, "/density")
            viscosity = hdf5_read_array(data, "/viscosity")
            temperature = hdf5_read_array(data, "/temperature")
            Tmask = np.array(temperature)
            Tmask[halfwidth <= min_width, :] = np.nan
            Twall = hdf5_read_array(data, "/Twall")
            qx = hdf5_read_array(data, "/qx")
            qy = hdf5_read_array(data, "/qy")
            alpha = hdf5_read_array(data, "/alpha")
            beta = hdf5_read_array(data, "/beta")
            gamma = hdf5_read_array(data, "/gamma")
            betaeq = hdf5_read_array(data, "/betaeq")
            Tliquidus = hdf5_read_array(data, "/Tliquidus")
            Tsolidus = hdf5_read_array(data, "/Tsolidus")
            ux = qx / (yb[1] - yb[0])
            time = hdf5_read_single(data, "/time")
            A = hdf5_read_array(data, "/A")
            C = hdf5_read_array(data, "/C")
            G = hdf5_read_array(data, "/G")
            heatflux = hdf5_read_array(data, "/MagmaToRockHeatFlux")
            shearheat = hdf5_read_array(data, "/ShearHeat")
            shear_heat_rate = np.mean(shearheat, axis=1) / (xr - xl)
            self.time.append(time)
            result = dict(
                xc=xc,
                xb=xb,
                yc=yc,
                yb=yb,
                halfwidth=halfwidth,
                width=width,
                pressure=pressure,
                overpressure=overpressure,
                density=density,
                temperature=temperature,
                Tmask=Tmask,
                viscosity=viscosity,
                Twall=Twall,
                qx=qx,
                qy=qy,
                alpha=alpha,
                beta=beta,
                gamma=gamma,
                betaeq=betaeq,
                Tliquidus=Tliquidus,
                Tsolidus=Tsolidus,
                ux=ux,
                A=A,
                C=C,
                G=-G,
                time=time,
                Qx=hdf5_read_array(data, "/TotalFluxElements"),
                Mx=hdf5_read_array(data, "/TotalMassRateElements"),
                tip_element=hdf5_read_int(data, "/TipElement"),
                tip_front=hdf5_read_single(data, "/TipFront"),
                heatflux=heatflux,
                shear_heat_rate=shear_heat_rate,
            )
            self.data.append(result)