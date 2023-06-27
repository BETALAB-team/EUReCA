
import os
import time
import random


import matplotlib
matplotlib.use('TkAgg')
matplotlib.interactive(True)
import matplotlib.pyplot as plt
import numpy as np
plt.style.use("ggplot")
import pandas as pd

#########################################################
# Config loading
# Loads a global config object
from eureca_building.config import load_config

config_path = os.path.join('.', 'config.json')
load_config(config_path)
from eureca_building.config import CONFIG

#########################################################

from eureca_building.weather import WeatherFile
from eureca_building.surface import Surface, SurfaceInternalMass
from eureca_building.internal_load import People, Lights, ElectricLoad
from eureca_building.ventilation import Infiltration, MechanicalVentilation, NaturalVentilation
from eureca_building.thermal_zone import ThermalZone
from eureca_building.air_handling_unit import AirHandlingUnit
from eureca_building.schedule import Schedule
from eureca_building.construction_dataset import ConstructionDataset
from eureca_building.construction import Construction
from eureca_building.setpoints import SetpointDualBand
from eureca_building.building import Building

#########################################################
# Epw loading
epw_path = os.path.join('..', 'example_scripts', 'ITA_Venezia-Tessera.161050_IGDG.epw')
weather_file = WeatherFile(epw_path,
                           time_steps=CONFIG.ts_per_hour,
                           azimuth_subdivisions=CONFIG.azimuth_subdivisions,
                           height_subdivisions=CONFIG.height_subdivisions, )
#########################################################
path = os.path.join(
    "..",
    "example_scripts",
    "materials_and_construction_test.xlsx",
)
# Define some constructions
dataset = ConstructionDataset.read_excel(path)
roof_cs = dataset.constructions_dict[13]
ceiling_cs = dataset.constructions_dict[17]
floor_cs = dataset.constructions_dict[16]
ext_wall_cs = dataset.constructions_dict[14]
int_wall_cs = dataset.constructions_dict[18]
window_cs = dataset.windows_dict[2]
mat_cs = dataset.materials_dict[1]

ext_wall_from_U = Construction.from_U_value("ExtWall from U", 0.7,weight_class="Medium", construction_type="ExtWall")
#########################################################

# Definition of surfaces
wall_south = Surface(
    "Wall 1",
    vertices=((0, 0, 0), (21.36, 0, 0), (21.36, 0, 6), (0, 0, 6)),
    wwr=0.125,
    surface_type="ExtWall",
    construction=ext_wall_from_U,
    window=window_cs,
    n_window_layers=2
)
wall_east = Surface(
    "Wall 2",
    vertices=((21.36, 0, 0), (21.36, 20.42, 0), (21.36, 20.42, 6), (21.36, 0, 6)),
    wwr=0.125,
    surface_type="ExtWall",
    construction=ext_wall_cs,
    window=window_cs,
    n_window_layers=2
)
wall_north = Surface(
    "Wall 2",
    vertices=((21.36, 20.42, 0), (0, 20.42, 0), (0, 20.42, 6), (21.36, 20.42, 6)),
    wwr=0.125,
    surface_type="ExtWall",
    construction=ext_wall_cs,
    window=window_cs,
)
wall_west = Surface(
    "Wall 2",
    vertices=((0, 20.42, 0), (0, 0, 0), (0, 0, 6), (0, 20.42, 6)),
    wwr=0.125,
    surface_type="ExtWall",
    construction=ext_wall_cs,
    window=window_cs,
)
floor = Surface(
    "Floor",
    vertices=((0, 0, 0), (0, 20.42, 0), (21.36, 20.42, 0), (21.36, 0, 0)),
    wwr=0.0,
    surface_type="GroundFloor",
    construction=floor_cs,
)
roof = Surface(
    "Roof",
    vertices=((0, 0, 6), (21.36, 0, 6), (21.36, 20.42, 6), (0, 20.42, 6)),
    wwr=0.0,
    surface_type="Roof",
    construction=roof_cs,
)
intwall = SurfaceInternalMass(
    "IntWall",
    area=floor._area*2.5*2,
    surface_type="IntWall",
    construction=int_wall_cs,
)
intceiling = SurfaceInternalMass(
    "IntCeiling",
    area=floor._area,
    surface_type="IntCeiling",
    construction=ceiling_cs,
)

#########################################################
# Ventilation

natural_vent_sched = Schedule(
    "nat_vent_sche",
    "dimensionless",
    np.array(([.3] * 8 * 2 + [.5] * 2 * 2 + [.3] * 4 * 2 + [.5] * 10 * 2) * 365)[:-1],
)

# Option 1
inf_obj = NaturalVentilation(
    name='nat_vent',
    unit='%',
    nominal_value=90., # 90% moltiplicato per il vettore della schedule
    schedule=natural_vent_sched,
    weather = weather_file,
    surfaces_with_opening = [wall_south, wall_east]
)

# Option 2
inf_obj = NaturalVentilation(
    name='nat_vent',
    unit='%',
    nominal_value=90., # 90% moltiplicato per il vettore della schedule
    schedule=natural_vent_sched,
)

weather_file.hourly_data["wind_direction"] = np.linspace(0, 360, 17519)

inf_obj.define_pressure_coef(
    weather = weather_file,
    surfaces_with_opening = [wall_south, wall_east]
)

###### Solution
# Summer week

vect_t_zona = []
vect_v = []
t_zona = 23.
z_n_tot = []
for t in range(5000 * 2, 5000 * 2 + 100):
    t_zona += 1.*random.randint(-100,100)/100
    z_n = inf_obj.get_timestep_ventilation_mass_flow(t, t_zona, weather_file)

    vect_t_zona.append(t_zona)
    z_n_tot.append(z_n)

z_n_tot = np.array(z_n_tot)

fig, [ax11, ax12, ax13, ax14] = plt.subplots(ncols=1, nrows=4)
ax11_ = ax11.twinx()
ax11.plot(vect_t_zona)
ax11.plot(weather_file.hourly_data["out_air_db_temperature"][5000 * 2: 5000 * 2 + 100])
ax11_.plot(weather_file.hourly_data["wind_speed"][5000 * 2: 5000 * 2 + 100])
ax11.set_ylabel("Zone temperature [°C]")
ax12.plot(z_n_tot, color = "r")
ax12.set_ylabel("Neutral plane height\n[-----]")


ax13.plot(wall_south.wind_pressure_coeff, color = "r")
ax13.plot(wall_east.wind_pressure_coeff, color = "b")
ax13.set_ylabel("Natural ventilation\nmass flow rate [-----]")
