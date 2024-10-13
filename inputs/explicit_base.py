import os, sys
repository_dir = os.path.abspath(os.getcwd())
sys.path.append(repository_dir)

import numpy as np
import json
import h5py
from py_scripts.input_utils import *

simID = 322
sim_dir = repository_dir + "/simulations/simID{}".format(simID)
inputDict = dict()
inputDict["simID"] = simID
inputDict["simDir"] = sim_dir
inputDict["algorithmProperties"] = {
    "isDebug": True,
    "flowModel" : "dike", # "channel", "dike"
    "timestepScheme" : "explicit",
    "numberOfLayers" : 30,
    "cutoffVelocity" : 1e-5,
    "lubricationCflScheme" : "spectral",
    "lubricationCflFactor" : 0.001,
    "massBalanceMinMobilityWidth" : 1e-10,
    "viscosityApproximation" : "mean", # "min", "harmonic", "mean"
    "shearHeating" : False,
    "latentHeatCrystallization" : True,
    "highOrderApproximation" : False,
    "solverName" : "umfpack", # "denselu", "umfpack", "pardiso"
}
inputDict["reservoirProperties"] = {
    "E": 20e9,
    "nu": 0.25,
    "g": 10.0,
    "KIc": 0.0,
    "specificHeatCapacity" : 1200,
    "thermalConductivity" : 2,
    "densityModel" : ReservoirDensity.constant,
    "temperatureModel" : ReservoirTemperature.constant_gradient,
    "numberOfLayers" : 10,
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
    "densityModel" : MagmaDensity.water_saturated,
    "viscosityModel" : ViscosityModel.grdmodel08,
    "crystallizationModel" : CrystallizationModel.const_relaxation,
    "saturationModel": MagmaSaturationModel.lavallee2015,
    "constantDensity" : {
        "rho" : 2000
    },
    "waterSaturatedDensity" : {
        "rhom0" : 2300.0,
        "rhow0" : 852.0,
        "rhoc0" : 2700.0,
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
        "tau" : 24 * 3600 * 1,
    }
}
inputDict["scheduleProperties"] = {
    "Q": [1.0, 0.0],
    "t": [0.0, 10000],
    "rho": 2000.0,
    "T" : 900,
    "beta" : 0.0,
}
inputDict["timestepProperties"] = {
    "startTime" : 0.0,
    "endTime" : 200000.0,
    "dtList" : [1.0/2.0, 2.0/2.0, 5.0/2.0],
    "dtTime" : [0.0, 10000.0, 30000.0],
    "outputSaveRate" : 100
}
inputDict["meshProperties"] = {
    "n" : 400,
    "xmin" : -30000.0,
    "xmax" : -10000.0
}

input_file = repository_dir + "/input.json"
with open(input_file, "w") as f: json.dump(inputDict, f, indent=4)