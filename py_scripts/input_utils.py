from enum import StrEnum


class ReservoirDensity(StrEnum):
    constant = "constant"
    
    
class ReservoirTemperature(StrEnum):
    constant_gradient = "constant_gradient"
    
    
class MagmaDensity(StrEnum):
    constant = "constant"
    water_saturated = "mixed_h2o_co2"
    

class ViscosityModel(StrEnum):
    constant = "constant"
    vft_const_coeff = "vft_const_coeff"
    vft_const_coeff_avg = "vft_const_coeff_avg"
    vft_const_coeff_cryst = "vft_const_coeff_cryst"
    grdmodel08 = "grdmodel08" # Giordano D, Russell JK, & Dingwell DB (2008), Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science, Letters, v. 271, 123-134.


class CrystallizationModel(StrEnum):
    equilibrium = "equilibrium"
    const_relaxation = "const_relaxation"


class MagmaSaturationModel(StrEnum):
    lavallee2015 = "lavallee2015"