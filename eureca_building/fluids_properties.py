"""
This module includes dictionaries with fluid properties (air, water, vapour, fuels)
To be used all over the code
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

#########################################
# RECOMMENDATION: Use unit in units.py  #
#########################################
air_properties = {
    "density": 1.2,  # kg/m3
    "specific_heat": 1000.,  # J/kgK
    "atmospheric_pressure": 101325.,  # Pa
}

water_properties = {
    "density": 1000.,  # kg/m3
    "specific_heat": 4182.,  # J/kgK
}

vapour_properties = {
    "latent_heat": 2501000.,  # J/kg
    "specific_heat": 1875.,  # J/kgK
}

fuels_pci = {
    "Natural Gas": 9970,   # Wh/Nm3
    "Oil": 9908, # Wh/L
    "Wood": 5277, # Wh/kg
    "Pellets": 4722, # Wh/kg
    "Coal": 7152, # Wh/kg
    "Gasoline": 14206, # Wh/L
    "LPG": 12791, # Wh/kg
}

gravitational_acceleration = 9.81  # m/s2