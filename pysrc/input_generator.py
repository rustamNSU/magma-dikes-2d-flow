from pathlib import Path
import json
from dataclasses import dataclass, asdict
from typing import ClassVar
from enum import StrEnum
from abc import ABC, abstractmethod
import shutil

class ViscosityApproximation(StrEnum):
    minimum = "min"
    harmonic = "harmonic"
    mean = "mean"

class SolverName(StrEnum):
    denselu = "denselu"
    umfpack = "umfpack"
    pardiso = "pardiso"


@dataclass
class AlgorithmProperties:
    is_debug: bool
    number_of_layers: int
    min_width: float
    min_mobility_width: float
    viscosity_approximation: str
    shear_heating: bool
    latent_heat_crystallization: bool
    solver_name: str
    is_sparse_elasticity: bool
    is_cohesive_stress: bool
    
    def create_dict(self) -> dict:
        return {
            "is_debug": self.is_debug,
            "number_of_layers": self.number_of_layers,
            "min_width": self.min_width,
            "min_mobility_width": self.min_mobility_width,
            "viscosity_approximation": self.viscosity_approximation,
            "shear_heating": self.shear_heating,
            "latent_heat_crystallization": self.latent_heat_crystallization,
            "solver_name": self.solver_name,
            "is_sparse_elasticity": self.is_sparse_elasticity,
            "is_cohesive_stress": self.is_cohesive_stress
        }


@dataclass
class ConstantDensity:
    rho: float


@dataclass
class ConstantTemperatureGradient:
    dT: float
    maximum_temperature: int
    minimum_temperature: int


@dataclass
class ReservoirProperties:
    youngs_modulus: float              # E in Pa
    poisson_ratio: float               # ν (nu)
    gravitational_acceleration: float  # Gravitational acceleration (m/s²)
    KIc: float                         # Fracture toughness (Pa·m^0.5)
    specific_heat_capacity: float
    thermal_conductivity: float
    number_of_layers: int
    reservoir_width: float
    density_model: 'ReservoirProperties.DensityModel'
    temperature_model: 'ReservoirProperties.TemperatureModel'
    mesh_refinement_algorithm: str = "cosine"
    
    def create_dict(self) -> dict:
        return {
            "youngs_modulus": self.youngs_modulus,
            "poisson_ratio": self.poisson_ratio,
            "gravitational_acceleration": self.gravitational_acceleration,
            "KIc": self.KIc,
            "specific_heat_capacity": self.specific_heat_capacity,
            "thermal_conductivity": self.thermal_conductivity,
            "number_of_layers": self.number_of_layers,
            "reservoir_width": self.reservoir_width,
            "density_model": self.density_model.create_dict(),
            "temperature_model": self.temperature_model.create_dict(),
            "mesh_refinement_algorithm": self.mesh_refinement_algorithm
        }
    
    
    class DensityModel(ABC):
        name: str

        @abstractmethod
        def create_dict(self) -> dict:
            pass
    
    @dataclass
    class ConstantDensity(DensityModel):
        rho: float
        name: ClassVar[str] = "constant_density"
        
        def create_dict(self) -> dict:
            return {"rho": self.rho, "name": self.name}
    
    
    class TemperatureModel(ABC):
        name: str

        @abstractmethod
        def create_dict(self) -> dict:
            pass
    
    @dataclass
    class ConstantTemperatureGradient(TemperatureModel):
        temperature_gradient: float
        maximum_temperature: float
        minimum_temperature: float
        name: ClassVar[str] = "constant_temperature_gradient"
        
        def create_dict(self) -> dict:
            return {
                "temperature_gradient": self.temperature_gradient,
                "maximum_temperature": self.maximum_temperature,
                "minimum_temperature": self.minimum_temperature,
                "name": self.name
            }
            

@dataclass
class MagmaProperties:
    thermal_conductivity: float
    specific_heat_capacity: float
    latent_heat: float
    density_model: 'MagmaProperties.DensityModel'
    viscosity_model: 'MagmaProperties.ViscosityModel'
    crystallization_model: 'MagmaProperties.CrystallizationModel'

    def create_dict(self) -> dict:
        return {
            "thermal_conductivity": self.thermal_conductivity,
            "specific_heat_capacity": self.specific_heat_capacity,
            "latent_heat": self.latent_heat,
            "density_model": self.density_model.create_dict(),
            "viscosity_model": self.viscosity_model.create_dict(),
            "crystallization_model": self.crystallization_model.create_dict(),
        }
    
    class DensityModel(ABC):
        name: str

        @abstractmethod
        def create_dict(self) -> dict:
            pass
    
    @dataclass
    class ConstantDensity(DensityModel):
        rho: float
        name: ClassVar[str] = "constant_density"
        
        def create_dict(self) -> dict:
            return {"rho": self.rho, "name": self.name}
    
    @dataclass
    class MixedH2OCO2(DensityModel):
        melt_density: float
        dissolved_h2o_density: float
        dissolved_co2_density: float
        crystal_density: float
        dissolved_data_path: str
        gas_density_data_path: str
        name: ClassVar[str] = "mixed_h2o_co2"
        
        def create_dict(self) -> dict:
            return {
                "melt_density": self.melt_density,
                "dissolved_h2o_density": self.dissolved_h2o_density,
                "dissolved_co2_density": self.dissolved_co2_density,
                "crystal_density": self.crystal_density,
                "dissolved_data_path": self.dissolved_data_path,
                "gas_density_data_path": self.gas_density_data_path,
                "name": self.name
            }
    
    
    class ViscosityModel(ABC):
        name: str

        @abstractmethod
        def create_dict(self) -> dict:
            pass
    
    @dataclass
    class VftConstantViscosity(ViscosityModel):
        A: float
        B: float
        C: float
        maximum_viscosity: float
        name: ClassVar[str] = "vft_const_coeff"
        
        def create_dict(self) -> dict:
            return {
                "A": self.A,
                "B": self.B,
                "C": self.C,
                "maximum_viscosity": self.maximum_viscosity,
                "name": self.name
            }
    
    @dataclass
    class ConstantViscosity(ViscosityModel):
        mu: float
        name: ClassVar[str] = "constant_viscosity"
        
        def create_dict(self) -> dict:
            return {"mu": self.mu, "name": self.name}
    
    @dataclass
    class GRDModel08(ViscosityModel):
        composition: list[float]
        maximum_viscosity: float = 1e14
        name: ClassVar[str] = "grdmodel08"
        
        def create_dict(self) -> dict:
            return {
                "composition": self.composition,
                "maximum_viscosity": self.maximum_viscosity,
                "name": self.name
            }
    
    
    class CrystallizationModel(ABC):
        name: str

        @abstractmethod
        def create_dict(self) -> dict:
            pass
    
    @dataclass
    class EquilibriumCrystallization(CrystallizationModel):
        """
        β = β_eq
        """
        name: ClassVar[str] = "equilibrium_crystallization"
        
        def create_dict(self) -> dict:
            return {"name": self.name}
        
        
    @dataclass
    class ConstantRelaxationCrystallization(CrystallizationModel):
        """
        dβ/dt = (β_eq - β) / τ
        """
        tau: float
        name: ClassVar[str] = "constant_relaxation_crystallization"
        
        def create_dict(self) -> dict:
            return {"tau": self.tau, "name": self.name}
        
    @dataclass
    class ArrheniusRelaxationCrystallization(CrystallizationModel):
        """
        Arrhenius-based crystallization model:
            τ = τ₀ exp(E/(R*T))
        where:
            tau0: Pre-exponential factor,
            E   : Activation energy [J/mol],
            R   : Gas constant [J/(mol*K)],
            T   : Temperature in Kelvin.
        """
        tau0: float
        E: float  # Activation energy [J/mol]
        R: float = 8.314  # Gas constant [J/(mol*K)]
        name: ClassVar[str] = "arrhenius_relaxation_crystallization"
        
        def create_dict(self) -> dict:
            return {
                "tau0": self.tau0, 
                "E" : self.E,
                "R" : self.R,
                "name": self.name
            }


@dataclass
class ScheduleProperties:
    volume_rate: list[float]
    time: list[float]
    density: float
    temperature: int
    beta: float

    def create_dict(self) -> dict:
        return {
            "volume_rate": self.volume_rate,
            "time": self.time,
            "density": self.density,
            "temperature": self.temperature,
            "beta": self.beta,
        }


@dataclass
class TimestepProperties:
    start_time: float
    end_time: float
    dt_list: list[float]
    dt_time: list[float]
    saverate_list: list[int]

    def create_dict(self) -> dict:
        return {
            "start_time": self.start_time,
            "end_time": self.end_time,
            "dt_list": self.dt_list,
            "dt_time": self.dt_time,
            "saverate_list": self.saverate_list,
        }


@dataclass
class MeshProperties:
    n: int
    xmin: float
    xmax: float

    def create_dict(self) -> dict:
        return {
            "n": self.n,
            "xmin": self.xmin,
            "xmax": self.xmax,
        }


class SimulationInput:
    def __init__(
        self, 
        sim_id: int, 
        repository_dir: str,
        algorithm_properties: AlgorithmProperties,
        reservoir_properties: ReservoirProperties,
        magma_properties: MagmaProperties,
        schedule_properties: ScheduleProperties,
        timestep_properties: TimestepProperties,
        mesh_properties: MeshProperties,
    ):
        self.sim_id = sim_id
        self.repository_dir = repository_dir
        self.simulation_dir = repository_dir + f'/simulations/simID{sim_id}'
        self.input_path = self.simulation_dir + '/input.json'
        self.algorithm_properties = algorithm_properties
        self.reservoir_properties = reservoir_properties
        self.magma_properties = magma_properties
        self.schedule_properties = schedule_properties
        self.timestep_properties = timestep_properties
        self.mesh_properties = mesh_properties
    
    def write_data(self, copy_input_for_example: bool=False):
        # Ensure the simulation directory exists.
        sim_dir = Path(self.simulation_dir)
        sim_dir.mkdir(parents=True, exist_ok=True)

        input_dict = {
            "sim_id": self.sim_id,
            "simulation_dir": self.simulation_dir,
            "repository_dir": self.repository_dir,
            "algorithm_properties": self.algorithm_properties.create_dict(),
            "reservoir_properties": self.reservoir_properties.create_dict(),
            "magma_properties": self.magma_properties.create_dict(),
            "schedule_properties": self.schedule_properties.create_dict(),
            "timestep_properties": self.timestep_properties.create_dict(),
            "mesh_properties": self.mesh_properties.create_dict()
        }
        
        # Write the JSON file with an indent for readability.
        with open(self.input_path, "w") as f:
            json.dump(input_dict, f, indent=4)
        
        if copy_input_for_example:
            example_dir = Path(self.repository_dir) / "example"
            example_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy(self.input_path, example_dir / f"input_simID{self.sim_id}.json")
            
        runner_file = self.repository_dir + "/runner.json"
        with open(runner_file, "w") as f:
            runnerDict = {
                "input_paths" : [
                    str(self.input_path)
                ]
            }
            json.dump(runnerDict, f, indent=4)