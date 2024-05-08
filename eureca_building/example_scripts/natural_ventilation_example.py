
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
    n_window_layers=2
)
wall_west = Surface(
    "Wall 2",
    vertices=((0, 20.42, 0), (0, 0, 0), (0, 0, 6), (0, 20.42, 6)),
    wwr=0.125,
    surface_type="ExtWall",
    construction=ext_wall_cs,
    window=window_cs,
    n_window_layers=2
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

#weather_file.hourly_data["wind_speed"] = weather_file.hourly_data["wind_speed"]*0

natural_vent_sched = Schedule(
    "nat_vent_sched",
    "dimensionless",
    np.array(([.3] * 8 * 2 + [.5] * 2 * 2 + [.3] * 4 * 2 + [.5] * 10 * 2) * 365)[:-1],
)

# Option 1
inf_obj = NaturalVentilation(
    name='nat_vent',
    unit='%',
    nominal_value=20., # 90% moltiplicato per il vettore della schedule
    schedule=natural_vent_sched,
    weather = weather_file,
    # surfaces_with_opening = [wall_west, wall_east],
    surfaces_with_opening = [wall_north, wall_south],
)

# # Option 2
# inf_obj = NaturalVentilation(
#     name='nat_vent',
#     unit='%',
#     nominal_value=.9, # 90% moltiplicato per il vettore della schedule
#     schedule=natural_vent_sched,
# )

# weather_file.hourly_data["wind_direction"] = np.linspace(0, 360, 17519)

# inf_obj.define_pressure_coef(
#     weather = weather_file,
#     # surfaces_with_opening = [wall_west, wall_east],
#     surfaces_with_opening = [wall_north, wall_south],
# )

###### Solution
# Summer week

vect_t_zona = []
vect_v = []
t_zona = 23.
z_n_tot = []
mass_flow_rate_tot = []
vol_inflow_tot = []
vol_outflow_tot = []
vol_flow_aggr = []
ts_to_sim = 25
for t in range(5000 * 2, 5000 * 2 + ts_to_sim):
    t_zona += 1.*random.randint(-100,100)/100
    # t_zona = weather_file.hourly_data["out_air_db_temperature"][t]
    x = inf_obj.get_timestep_ventilation_mass_flow(t, t_zona, weather_file, 5)
    mass_flow_rate = x[0]
    z_n = x[1]
    vol_inflow = x[2]
    vol_outflow = x[3]
    vol_flow = x[4]
    print(mass_flow_rate)
    print(z_n)
    print(vol_inflow)
    print(vol_outflow)
    print(vol_flow)
    print('')

    vect_t_zona.append(t_zona)
    z_n_tot.append(z_n)
    mass_flow_rate_tot.append(mass_flow_rate)
    vol_inflow_tot.append(vol_inflow)
    vol_outflow_tot.append(vol_outflow)
    vol_flow_aggr.append(vol_flow)
   
wall_south._h_bottom_windows
z_n_tot = np.array(z_n_tot)
mass_flow_rate_tot = np.array(mass_flow_rate_tot)

x_lim = [0,24]
x_lim_flow = [0,24]

fig, [ax11, ax12, ax13, ax14] = plt.subplots(ncols=1, nrows=4)
ax11_ = ax11.twinx()
ax11.plot(vect_t_zona)
ax11.plot(weather_file.hourly_data["out_air_db_temperature"][5000 * 2: 5000 * 2 + ts_to_sim])
ax11_.plot(weather_file.hourly_data["wind_speed"][5000 * 2: 5000 * 2 + ts_to_sim])
ax11.set_ylabel("Zone temperature [°C]")
ax11.set_xlim(x_lim)
ax12.hlines(wall_south._h_bottom_windows, xmin=0, xmax=len(z_n_tot)-1,color = "b")
ax12.hlines(wall_south._h_top_windows, xmin=0, xmax=len(z_n_tot)-1,color = "k")
ax12.plot(z_n_tot, color = "r")
ax12.set_ylabel("Neutral plane height\n[-----]")
ax12.set_xlim(x_lim)


#ax13.plot(vol_flow_rate_tot[:], color = "b")
#ax13.plot(vol_flow_rate_tot[:,1], color = "r")
#ax13.set_ylim([-0.04,0.04])
ax13.set_xlim(x_lim)
ax13.plot(mass_flow_rate_tot, color = "y")
ax13.set_ylabel("Natural ventilation\nmass flow rate [kg/s]")

ax14.set_xlim(x_lim)
ax14.plot(wall_north.wind_pressure_coeff[5000 * 2: 5000 * 2 + ts_to_sim], color = "b")
ax14.plot(wall_south.wind_pressure_coeff[5000 * 2: 5000 * 2 + ts_to_sim], color = "k")
ax14.set_ylabel("Natural ventilation\nmass flow rate [-----]")

# print(np.max(np.abs(vol_flow_rate_sopra_tot)))

fig, [ax11, ax12] = plt.subplots(ncols=1, nrows=2)
ax11.plot(vol_inflow_tot)
ax11.set_xlim(x_lim_flow)
ax12.plot(vol_outflow_tot)
ax12.set_xlim(x_lim_flow)

fig, [ax11, ax12] = plt.subplots(ncols=1, nrows=2)
ax11_ = ax11.twinx()
ax11.plot(vect_t_zona)
ax11.plot(weather_file.hourly_data["out_air_db_temperature"][5000 * 2: 5000 * 2 + ts_to_sim])
ax11_.plot(weather_file.hourly_data["wind_speed"][5000 * 2: 5000 * 2 + ts_to_sim])
ax11.set_ylabel("Zone temperature [°C]")
ax11.set_xlim(x_lim)
# ax12.hlines(wall_south._h_bottom_windows, xmin=0, xmax=len(z_n_tot)-1,color = "b")
# ax12.hlines(wall_south._h_top_windows, xmin=0, xmax=len(z_n_tot)-1,color = "k")
# ax12.plot(z_n_tot, color = "r")
# ax12.set_ylabel("Neutral plane height\n[-----]")
# ax12.set_xlim(x_lim)


#ax13.plot(vol_flow_rate_tot[:], color = "b")
#ax13.plot(vol_flow_rate_tot[:,1], color = "r")
#ax13.set_ylim([-0.04,0.04])
ax12.set_xlim(x_lim)
ax12.plot(mass_flow_rate_tot, color = "y")
ax12.set_ylabel("Natural ventilation\nmass flow rate [kg/s]")


# ax14.plot(wall_west.wind_pressure_coeff, color = "b")
# ax14.plot(wall_east.wind_pressure_coeff, color = "k")
# ax14.set_ylabel("Natural ventilation\nmass flow rate [-----]")
