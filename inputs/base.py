import os, sys
repository_dir = os.path.abspath(os.getcwd())
sys.path.append(repository_dir)

import numpy as np
import json
import h5py
from py_scripts.input_utils import *

simID = 1
sim_dir = repository_dir + "/simulations/simID{}".format(simID)
inputDict = dict()
inputDict["simID"] = simID
inputDict["simDir"] = sim_dir
inputDict["algorithmProperties"] = {
    "isDebug": True,
    "numberOfLayers" : 30,
    "minWidth" : 1e-5,
    "minMobilityWidth" : 1e-4,
    "viscosityApproximation" : "harmonic", # "min", "harmonic", "mean"
    "shearHeating" : True,
    "latentHeatCrystallization" : True,
    "solverName" : "umfpack", # "denselu", "umfpack", "pardiso"
    "isSparseElasticity" : True,
    "isCohesiveStress" : False,
}
inputDict["reservoirProperties"] = {
    "E": 20e9,
    "nu": 0.25,
    "g": 10.0,
    "KIc": 1.0e6,
    "specificHeatCapacity" : 1200,
    "thermalConductivity" : 2.0,
    "densityModel" : ReservoirDensity.constant,
    "temperatureModel" : ReservoirTemperature.constant_gradient,
    "numberOfLayers" : 30,
    "reservoirWidth" : 3.0,
    "meshRefinementAlgorithm" : "cosine",
    "constantDensity" : {
        "rho" : 2700.0
    },
    "constantTemperatureGradient" : {
        "dT" : 30e-3,
        "maximum_temperature" : 900,
        "minimum_temperature" : 0
    },
}
inputDict["magmaProperties"] = {
    "thermalConductivity" : 2,
    "specificHeatCapacity" : 1200,
    "latentHeat" : 350000,
    "densityModel" : MagmaDensity.mixed_h2o_co2,
    "viscosityModel" : ViscosityModel.grdmodel08,
    "crystallizationModel" : CrystallizationModel.const_relaxation,
    "saturationModel": MagmaSaturationModel.mixed_h2o_co2,
    "constantDensity" : {
        "rho" : 2000,
    },
    "mixed_h2o_co2" : {
        "rhom0" : 2360.0,
        "rhoh2o0" : 852.0,
        "rhoco20" : 2124,
        "rhoc0" : 2700.0,
        "dissolved_data_path" : "./data/pinatubo/dissolved_06.json",
        "gas_density_data_path" : "./data/pinatubo/h2o_co2_gas_density.json"
    },
    "vftConstantViscosity" : {
        "A" : -4.55,
        "B" : 7455,
        "C" : 180,
        "muMaxLimit" : 1e14
    },
    "constantViscosity" : {
        "mu" : 1000
    },
    "grdmodel08" : {
        # SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 F2O-1
        "composition" : [64.6, 0.53, 16.5, 4.47, 0.01, 2.39, 5.23, 4.49, 1.54, 0.04, 0],
        "muMaxLimit" : 1e14
    },
    "linearViscosity" : {
        "much" : 50,
        "musurf" : 1000,
    },
    "constantRelaxationCrystallization" : {
        "tau" : 24 * 3600 * 0.5,
    }
}
inputDict["scheduleProperties"] = {
    "Q": [1.0, 0.0],
    "t": [0.0, 10000],
    "rho": 2000.0,
    "T" : 800,
    "beta" : 0.0,
}
inputDict["timestepProperties"] = {
    "startTime" : 0.0,
    "endTime" : 200000.0,
    "dtList" : [1.0, 2.0, 5.0],
    "dtTime" : [0.0, 20000.0, 40000.0],
    "saverateList" : [100, 100, 100]
}
inputDict["meshProperties"] = {
    "n" : 300,
    "xmin" : -30000.0,
    "xmax" : 0.0
}

input_file = repository_dir + "/input.json"
with open(input_file, "w") as f: json.dump(inputDict, f, indent=4)