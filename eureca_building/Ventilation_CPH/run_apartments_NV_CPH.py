"""
List of custom exceptions
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import os
import time

import datetime as dt
import matplotlib
matplotlib.use('TkAgg')
matplotlib.interactive(True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

#########################################################
# Config loading
# Loads a global config object
from eureca_building.config import load_config

config_path = os.path.join('.', 'config_NV_CPH.json')
load_config(config_path)
from eureca_building.config import CONFIG

#########################################################

from eureca_building.weather import WeatherFile
from eureca_building.material import Material
from eureca_building.surface import Surface, SurfaceInternalMass
from eureca_building.internal_load import People, Lights, ElectricLoad
from eureca_building.ventilation import Infiltration, MechanicalVentilation
from eureca_building.thermal_zone import ThermalZone
from eureca_building.air_handling_unit import AirHandlingUnit
from eureca_building.schedule import Schedule
from eureca_building.construction_dataset import ConstructionDataset
from eureca_building.construction import Construction
from eureca_building.setpoints import SetpointDualBand
from eureca_building.building import Building
from eureca_building.domestic_hot_water import DomesticHotWater
from eureca_building.window import SimpleWindow
from eureca_building.ventilation import NaturalVentilation

import subprocess

#########################################################
# Running consecutively the simulations for the apartments in Copenhagen
print("Starting apartments' simulations")
start = time.time()

# print("Simulation for apartment L07")
# subprocess.call("./eureca_example_NV_CPH_L07.py", shell=True)
# exec(open("eureca_example_NV_CPH_L07.py").read())
# os.system('python .\\eureca_example_NV_CPH_L07.py')

print("Simulation for apartment L09")
# subprocess.call("./eureca_example_NV_CPH_L09.py", shell=True)
# exec(open("eureca_example_NV_CPH_L09.py").read())
os.system('python .\\eureca_example_NV_CPH_L09.py')

# print("Simulation for apartment L10")
# subprocess.call("./eureca_example_NV_CPH_L10.py", shell=True)
# exec(open("eureca_example_NV_CPH_L10.py").read())
# os.system('python .\\eureca_example_NV_CPH_L10.py')

print("Simulation for apartment L12")
# subprocess.call("./eureca_example_NV_CPH_L12.py", shell=True)
# exec(open("eureca_example_NV_CPH_L12.py").read())
os.system('python .\\eureca_example_NV_CPH_L12.py')

print("Simulation for apartment L13")
# subprocess.call("./eureca_example_NV_CPH_L13.py", shell=True)
# exec(open("eureca_example_NV_CPH_L13.py").read())
os.system('python .\\eureca_example_NV_CPH_L13.py')

print("Simulation for apartment L14")
# subprocess.call("./eureca_example_NV_CPH_L14.py", shell=True)
# exec(open("eureca_example_NV_CPH_L14.py").read())
os.system('python .\\eureca_example_NV_CPH_L14.py')

print("Simulation for apartment L15")
# subprocess.call("./eureca_example_NV_CPH_L15.py", shell=True)
# exec(open("eureca_example_NV_CPH_L15.py").read())
os.system('python .\\eureca_example_NV_CPH_L15.py')

# print("Simulation for apartment L17")
# subprocess.call("./eureca_example_NV_CPH_L17.py", shell=True)
# exec(open("eureca_example_NV_CPH_L17.py").read())
# os.system('python .\\eureca_example_NV_CPH_L17.py')

# print("Simulation for apartment L18")
# subprocess.call("./eureca_example_NV_CPH_L18.py", shell=True)
# exec(open("eureca_example_NV_CPH_L18.py").read())
# os.system('python .\\eureca_example_NV_CPH_L18.py')

# print("Simulation for apartment L19")
# subprocess.call("./eureca_example_NV_CPH_L19.py", shell=True)
# exec(open("eureca_example_NV_CPH_L19.py").read())
# os.system('python .\\eureca_example_NV_CPH_L19.py')

print("Simulation for apartment L20 (Ali)")
# # subprocess.run("./eureca_example_NV_CPH_L20_Ali.py", shell=True)
# # exec(open("eureca_example_NV_CPH_L20_Ali.py").read())
os.system('python .\\eureca_example_NV_CPH_L20_Ali.py')

# print("Simulation for apartment L20 (Ali), 2C model")
# subprocess.run(".\\eureca_example_NV_CPH_L20_Ali_2C.py", shell=True)
# exec(open(".\\eureca_example_NV_CPH_L20_Ali_2C.py").read())
# os.system('python .\\eureca_example_NV_CPH_L20_Ali_2C.py')

stop = time.time()
print(f'Total calibration time: {(stop-start):.1f} s')