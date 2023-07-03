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


import matplotlib
matplotlib.use('TkAgg')
matplotlib.interactive(True)
import matplotlib.pyplot as plt
import numpy as np
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

#########################################################
# Epw loading
epw_path = os.path.join('..', 'example_scripts', 'ITA_Venezia-Tessera.161050_IGDG.epw')
weather_file = WeatherFile(epw_path,
                           time_steps=CONFIG.ts_per_hour,
                           azimuth_subdivisions=CONFIG.azimuth_subdivisions,
                           height_subdivisions=CONFIG.height_subdivisions, )
#########################################################



fig, ax = plt.subplots()
ax.set_ylabel("Specific Humidity [" + "$g_{v}/kg_{da}$"+"]")
ax.set_xlabel("Temperature ["+"$°C$"+"]")
ax.set_xlim(-10,45)
ax.set_ylim(0,30)
ax.set_title("Psychrometric chart")
t = np.arange(-10,46,1)
p_sat = 6.1094*np.exp(17.625*t/(t+243.04))*100
for ur in np.arange(0.1,1.1,0.1):
    p = p_sat * ur
    sh = 0.622*(p/(101325 - p))*1000
    ax.plot(t,sh, 'k-', linewidth=0.3)
    x_text = min([t[-1] - 5, 35])
    y_text = min([sh[-1] - 2, 28])

ax.text(27.7, 26.9, f"100%", backgroundcolor = "white", fontsize = 6, ma="center")
ax.text(30.5, 23.7, f"80%", backgroundcolor = "white", fontsize = 6, ma="center")
ax.text(32.5, 20., f"60%", backgroundcolor = "white", fontsize = 6, ma="center")
ax.text(34.5, 14.5, f"40%", backgroundcolor = "white", fontsize = 6, ma="center")
ax.text(35.8, 8, f"20%", backgroundcolor = "white", fontsize = 6, ma="center")

x_text, y_text = -9.,3.5

for h in np.arange(0.,200.,10.):
    x = (h - 1.006*t)/(1.86*t+2501)*1000
    ax.plot(t,x, 'k:', linewidth=0.3)
    if y_text < 30.:
        ax.text(x_text, y_text, f"{h:.0f}"+ " ["+"$kJ/kg_{da}$"+"]", backgroundcolor = "white", fontsize = 6)
    x_text += 3
    y_text += 2.7












"""
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
)
wall_east = Surface(
    "Wall 2",
    vertices=((21.36, 0, 0), (21.36, 20.42, 0), (21.36, 20.42, 6), (21.36, 0, 6)),
    wwr=0.125,
    surface_type="ExtWall",
    construction=ext_wall_cs,
    window=window_cs,
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
# Loads


# A schedule
people_sched = Schedule(
    "PeopleOccupancy1",
    "percent",
    np.array(([0.1] * 7 * 2 + [0.6] * 2 * 2 + [0.4] * 5 * 2 + [0.6] * 10 * 2) * 365)[:-1],
)

# Loads
people = People(
    name='occupancy_tz',
    unit='px',
    nominal_value=1.2,
    schedule=people_sched,
    fraction_latent=0.45,
    fraction_radiant=0.3,
    fraction_convective=0.7,
    metabolic_rate=150,
)

lights = Lights(
    name='lights_tz',
    unit='W/m2',
    nominal_value=5.,
    schedule=people_sched,
    fraction_radiant=0.7,
    fraction_convective=0.3,
)

pc = ElectricLoad(
    name='pc_tz',
    unit='W',
    nominal_value=300.,
    schedule=people_sched,
    fraction_radiant=0.2,
    fraction_convective=0.8,
)
#########################################################
# Setpoints
heat_t = Schedule.from_daily_schedule(
    name = "t_heat",
    schedule_type="temperature",
    schedule_week_day=np.array([18] * 7 * 2 + [21] * 2 * 2 + [18] * 5 * 2 + [21] * 10 * 2),
    schedule_saturday=np.array([18] * 7 * 2 + [21] * 2 * 2 + [18] * 5 * 2 + [21] * 10 * 2) -5,
    schedule_sunday=np.array([18] * 7 * 2 + [21] * 2 * 2 + [18] * 5 * 2 + [21] * 10 * 2) - 10,
    schedule_holiday=np.array([18] * 7 * 2 + [21] * 2 * 2 + [18] * 5 * 2 + [21] * 10 * 2) * 0,
    holidays = (10,11,12,13,14),
    starting_day = 3,
)

heat_t = Schedule.from_constant_value(
    name = "t_heat",
    schedule_type="temperature",
    value = 15.
)

heat_t = Schedule(
    "t_heat",
    "temperature",
    np.array(([18] * 7 * 2 + [21] * 2 * 2 + [18] * 5 * 2 + [21] * 10 * 2) * 365)[:-1],
)

cool_t = Schedule(
    "t_cool",
    "temperature",
    np.array(([28] * 8 * 2 + [26] * 2 * 2 + [28] * 4 * 2 + [26] * 10 * 2) * 365)[:-1],
)

heat_h = Schedule(
    "h_heat",
    "dimensionless",
    np.array(([0.1] * 7 * 2 + [0.3] * 2 * 2 + [.1] * 5 * 2 + [.3] * 10 * 2) * 365)[:-1],
)

cool_h = Schedule(
    "h_cool",
    "dimensionless",
    np.array(([.9] * 8 * 2 + [.5] * 2 * 2 + [.9] * 4 * 2 + [.5] * 10 * 2) * 365)[:-1],
)

t_sp = SetpointDualBand(
    "t_sp",
    "temperature",
    schedule_lower=heat_t,
    schedule_upper=cool_t,
)
h_sp = SetpointDualBand(
    "h_sp",
    "relative_humidity",
    schedule_lower=heat_h,
    schedule_upper=cool_h,
)
#########################################################
# Ventilation

infiltration_sched = Schedule(
    "inf_sched",
    "dimensionless",
    np.array(([.3] * 8 * 2 + [.5] * 2 * 2 + [.3] * 4 * 2 + [.5] * 10 * 2) * 365)[:-1],
)

inf_obj = Infiltration(
    name='inf_obj',
    unit='Vol/h',
    nominal_value=1.,
    schedule=infiltration_sched,
)


#########################################################
# Mechanical ventilation
vent_sched = Schedule(
    "vent_sched",
    "dimensionless",
    np.array(([.0] * 8 * 2 + [1] * 2 * 2 + [0] * 4 * 2 + [1] * 10 * 2) * 365)[:-1],
)

vent_obj = MechanicalVentilation(
    name='vent_obj',
    unit='Vol/h',
    nominal_value=1,
    schedule=vent_sched,
)

T_supply_sched = Schedule(
    "T_supply_sched",
    "temperature",
    np.array(([23.] * 8 * 2 + [23.] * 2 * 2 + [23.] * 4 * 2 + [23.] * 10 * 2) * 365)[:-1],
)

x_supply_sched = Schedule(
    "x_supply_sched",
    "dimensionless",
    np.array(([0.0101] * 8 * 2 + [0.0101] * 2 * 2 + [0.0101] * 4 * 2 + [0.0101] * 10 * 2) * 365)[:-1]*0.7,
)

availability_sched = np.array(([0] * 8 * 2 + [1] * 2 * 2 + [0] * 4 * 2 + [1] * 10 * 2) * 365)[:-1]
availability_sched[120*24*2:273*24*2] = -1*availability_sched[120*24*2:273*24*2]
ahu_availability_sched = Schedule(
    "ahu_availability_sched",
    "availability",
    availability_sched,
)



#########################################################

zones = []

for i in range(1):
    print(i)

    # Create zone
    tz1 = ThermalZone(
        name="Zone 1",
        surface_list=[wall_south, wall_east, wall_north, wall_west, roof, floor, intwall, intceiling],
        net_floor_area=floor._area*2,
        volume=floor._area*3.3*2)

    zones.append(tz1)
    tz1._ISO13790_params()
    tz1._VDI6007_params()

    tz1.add_internal_load(people)
    tz1.add_internal_load(lights, pc)

    # IHG preprocessing
    tz_loads = tz1.extract_convective_radiative_latent_load()
    tz1.calculate_zone_loads_ISO13790(weather_file)
    # tz1._plot_ISO13790_IHG()

    # 2C model
    tz1.calculate_zone_loads_VDI6007(weather_file)
    # tz1._plot_VDI6007_IHG(weather_file)

    tz1.add_temperature_setpoint(t_sp)
    tz1.add_humidity_setpoint(h_sp)

    # Natural Ventilation preprocessing
    tz1.add_infiltration(inf_obj)
    tz_inf = tz1.calc_infiltration(weather_file)

    ahu = AirHandlingUnit(
    "ahu",
    vent_obj,
    T_supply_sched,
    x_supply_sched,
    ahu_availability_sched,
    True,
    0.5,
    0.5,
    0.9,
    weather_file,
    tz1,
)

    cooling_1C_peak_load = tz1.design_sensible_cooling_load(weather_file, model = "1C")
    heating_peak_load = tz1.design_heating_load(-5.)

    bd = Building("Bd 1", thermal_zones_list=[tz1], model = "1C")
    start = time.time()
    bd.set_hvac_system("CondensingBoiler", "SplitAirConditioner")
    bd.set_hvac_system_capacity(weather_file)
    df_res = bd.simulate(weather_file, output_folder="Results")
    print(f"2C model: \n\t{8760 * 2 - 1} time steps\n\t{(time.time() - start):.2f} s")
"""