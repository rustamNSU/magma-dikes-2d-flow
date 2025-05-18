import json
import h5py
import numpy as np
import fnmatch, os
import re
from tqdm import tqdm
from dataclasses import dataclass, field
from typing import ClassVar

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


@dataclass
class DataFrame:
    xc: np.ndarray
    xb: np.ndarray
    yc: np.ndarray
    yb: np.ndarray
    halfwidth: np.ndarray
    width: np.ndarray
    pressure: np.ndarray
    overpressure: np.ndarray
    density: np.ndarray
    temperature: np.ndarray
    viscosity: np.ndarray
    wall_temperature: np.ndarray
    qx: np.ndarray
    qy: np.ndarray
    alpha: np.ndarray
    beta: np.ndarray
    gamma: np.ndarray
    betaeq: np.ndarray
    tau: np.ndarray
    liquidus_temperature: np.ndarray
    solidus_temperature: np.ndarray
    time: float
    A: np.ndarray
    C: np.ndarray
    G: np.ndarray
    Qx: np.ndarray
    Mx: np.ndarray
    tip_element: int
    tip_front: float
    heatflux: np.ndarray
    
    min_width: ClassVar[float] = 1e-5
    open_mask: np.ndarray = field(init=False)
    
    def __post_init__(self):
        # Compute a mask for open elements: True if halfwidth > min_width
        self.open_mask = self.halfwidth > DataFrame.min_width

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
        self.nx = int(self.input["mesh_properties"]["n"])
        self.xmin = float(self.input["mesh_properties"]["xmin"])
        self.xmax = float(self.input["mesh_properties"]["xmax"])
    
    
    def _preload_data(self):
        self.timesteps = range(0, self.ntimesteps, self.step_rate)
        self.data: list[DataFrame] = []
        self.time = []
        DataFrame.min_width = self.input["algorithm_properties"]["min_width"]
        for itime in tqdm(self.timesteps):
            file_path = f"{self.data_dir}/{self.filepaths[itime]}"
            with h5py.File(file_path, 'r') as d:
                df = DataFrame(
                    xc = hdf5_read_array(d, "/mesh/xc"),
                    xb = np.append(hdf5_read_array(d, "/mesh/xl"), hdf5_read_array(d, "/mesh/xr")[-1]),
                    yc = hdf5_read_array(d, "/yc"),
                    yb = hdf5_read_array(d, "/yb"),
                    halfwidth = hdf5_read_array(d, "/halfwidth"),
                    width = hdf5_read_array(d, "/width"),
                    pressure = hdf5_read_array(d, "/pressure"),
                    overpressure = hdf5_read_array(d, "/overpressure"),
                    density = hdf5_read_array(d, "/density"),
                    temperature = hdf5_read_array(d, "/temperature"),
                    viscosity = hdf5_read_array(d, "/viscosity"),
                    wall_temperature = hdf5_read_array(d, "/Twall"),
                    qx = hdf5_read_array(d, "/qx"),
                    qy = hdf5_read_array(d, "/qy"),
                    alpha = hdf5_read_array(d, "/alpha"),
                    beta = hdf5_read_array(d, "/beta"),
                    gamma = hdf5_read_array(d, "/gamma"),
                    betaeq = hdf5_read_array(d, "/betaeq"),
                    tau = hdf5_read_array(d, "/tau"),
                    liquidus_temperature = hdf5_read_array(d, "/Tliquidus"),
                    solidus_temperature = hdf5_read_array(d, "/Tsolidus"),
                    time = hdf5_read_single(d, "/time"),
                    A = hdf5_read_array(d, "/A"),
                    C = hdf5_read_array(d, "/C"),
                    G = -hdf5_read_array(d, "/G"),
                    Qx = hdf5_read_array(d, "/TotalFluxElements"),
                    Mx = hdf5_read_array(d, "/TotalMassRateElements"),
                    tip_element = hdf5_read_int(d, "/TipElement"),
                    tip_front = hdf5_read_single(d, "/TipFront"),
                    heatflux = hdf5_read_array(d, "/MagmaToRockHeatFlux")
                )
                self.data.append(df)
                self.time.append(df.time)