import os, sys
repository_dir = os.path.abspath(os.getcwd())
sys.path.append(repository_dir)

from pysrc import *
import numpy as np
import json
import h5py
from py_scripts.input_utils import *

hour = 3600.0
day = 24 * hour
week = 7 * day

simulation_input = SimulationInput(
    sim_id=157,
    repository_dir=repository_dir,
    algorithm_properties=AlgorithmProperties(
        is_debug=True,
        number_of_layers=30,
        min_width=1e-5,
        min_mobility_width=1e-4,
        viscosity_approximation=ViscosityApproximation.harmonic,
        shear_heating=True,
        latent_heat_crystallization=True,
        solver_name=SolverName.umfpack,
        is_sparse_elasticity=True,
        is_cohesive_stress=True
    ),
    reservoir_properties=ReservoirProperties(
        youngs_modulus=20e9,
        poisson_ratio=0.25,
        gravitational_acceleration=9.81,
        KIc=1e6,
        specific_heat_capacity=1200,
        thermal_conductivity=2.0,
        number_of_layers=30,
        reservoir_width=5.0,
        density_model=ReservoirProperties.ConstantDensity(
            rho=2700.0
        ),
        temperature_model=ReservoirProperties.ConstantTemperatureGradient(
            temperature_gradient=30e-3,
            maximum_temperature=900,
            minimum_temperature=0
        ),
    ),
    magma_properties=MagmaProperties(
        thermal_conductivity=2.0,
        specific_heat_capacity=1200,
        latent_heat=350000,
        density_model=MagmaProperties.MixedH2OCO2(
            melt_density=2360.0,
            dissolved_h2o_density=900.0,
            dissolved_co2_density=1400.0,
            crystal_density=2700.0,
            dissolved_data_path="./data/pinatubo/dissolved_06.json",
            gas_density_data_path="./data/pinatubo/h2o_co2_gas_density.json"
        ),
        viscosity_model=MagmaProperties.GRDModel08(
            composition=[64.6, 0.53, 16.5, 4.47, 0.01, 2.39, 5.23, 4.49, 1.54, 0.04, 0]
        ),
        crystallization_model=MagmaProperties.ConstantRelaxationCrystallization(
            tau=week
        )
    ),
    schedule_properties=ScheduleProperties(
        volume_rate=[1.0, 0.0, 1.0, 0.0],
        time=[0.0, 7500, week, week + 7500],
        density=2000.0,
        temperature=900,
        beta=0.0
    ),
    timestep_properties=TimestepProperties(
        start_time=0.0,
        end_time=week+5*day,
        dt_list=[1.0, 2.0, 4.0, 1.0, 2.0, 4.0],
        dt_time=[0.0, 4*hour, 2*day, week, week + 4*hour, week + 2*day],
        saverate_list=[100, 100, 100, 100, 100, 100]
    ),
    mesh_properties=MeshProperties(
        n=300,
        xmin=-30000.0,
        xmax=0.0
    )
)
simulation_input.write_data(True)