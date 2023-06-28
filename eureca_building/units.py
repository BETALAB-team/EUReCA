"""
List of physic properties with their unit
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

units = {
    "length": "[m]",
    "conductivity": "[W/(m K)]",
    "density": "[kg/m3]",
    "specific_heat": "[J/(kg K)]",
    "latent_heat": "[J/kg]",
    "specific_mass_flow_rate": "[kg/(s m2)]",
    "temperature": "[°C]",
    "thermal_resistance": "[(m2 K)/W]",
    "absorptance": "[-]",
    "U_value": "[W/(m2 K)]",
    "solar_heat_gain_coefficient": "[-]",
    "non_dimensional_coefficient": "[-]",
}

material_limits = {
    "thickness": [0.0, 1.0],
    "conductivity": [0.0, 100.0],
    "density": [0.0, 10000.0],
    "specific_heat": [0.0, 3000.0],
    "starting_temperature": [0.0, 99.0],
    "thermal_resistance": [0, 20],
    "absorptance": [0, 1],
}

window_material_limits = {
    "window_u_value": [1.0, 7.0],
    "solar_heat_gain_coefficient": [0.0, 1.0],
    "non_dimensional_coefficient": [0.0, 1.0],
}
