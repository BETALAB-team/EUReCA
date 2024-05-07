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

config_path = os.path.join('.', 'config_VDI.json')
load_config(config_path)
from eureca_building.config import CONFIG

#########################################################

from eureca_building.weather import WeatherFile
from eureca_building.surface import Surface, SurfaceInternalMass
from eureca_building.thermal_zone import ThermalZone
from eureca_building.internal_load import People, Lights, ElectricLoad
from eureca_building.ventilation import Infiltration, MechanicalVentilation
from eureca_building.schedule import Schedule
from eureca_building.construction_dataset import ConstructionDataset
from eureca_building.setpoints import SetpointDualBand
from eureca_building.air_handling_unit import AirHandlingUnit

TESTS = [1,3,
         2,4,
         5,
         6,
         7,
         12
    ]

#########################################################
# Epw loading
epw_path = os.path.join('..', 'example_scripts', 'ITA_Venezia-Tessera.161050_IGDG.epw')
weather_file = WeatherFile(epw_path,
                           time_steps=1, )
#########################################################
path = os.path.join(
    "..",
    "example_scripts",
    "materials_and_construction_test_VDI.xlsx",
)
# Define some constructions
dataset = ConstructionDataset.read_excel(path)
# TEST_L
extWall_cs_L = dataset.constructions_dict[6]
door_cs_L = dataset.constructions_dict[7]
int_wall_cs_L = dataset.constructions_dict[8]
ceiling_cs_L = dataset.constructions_dict[10]
floor_cs_L = dataset.constructions_dict[9]
window_cs_L = dataset.windows_dict[3]
# TEST_S
extWall_cs_S = dataset.constructions_dict[1]
door_cs_S = dataset.constructions_dict[2]
int_wall_cs_S = dataset.constructions_dict[3]
ceiling_cs_S = dataset.constructions_dict[5]
floor_cs_S = dataset.constructions_dict[4]
window_cs_S = dataset.windows_dict[2]
#########################################################

# Definition of surfaces
wall_south_L = Surface(
    "Wall 1",
    vertices=((0, 0, 0), (10.5, 0, 0), (10.5, 0, 1), (0, 0, 1)),
    wwr=7/10.5,
    surface_type="ExtWall",
    construction=extWall_cs_L,
    window=window_cs_L,
)
intwall_L = SurfaceInternalMass(
    "IntWall",
    area=38.5,
    surface_type="IntWall",
    construction=int_wall_cs_L,
)
intdoor_L = SurfaceInternalMass(
    "IntWall",
    area=2.,
    surface_type="IntWall",
    construction=door_cs_L,
)
intceiling_L = SurfaceInternalMass(
    "IntCeiling",
    area=17.5,
    surface_type="IntCeiling",
    construction=ceiling_cs_L,
)
intfloor_L = SurfaceInternalMass(
    "IntFloor",
    area=17.5,
    surface_type="IntFloor",
    construction=floor_cs_L,
)

# TEST_S
wall_south_S = Surface(
    "Wall 1",
    vertices=((0, 0, 0), (10.5, 0, 0), (10.5, 0, 1), (0, 0, 1)),
    wwr=7/10.5,
    surface_type="ExtWall",
    construction=extWall_cs_S,
    window=window_cs_S,
)
intwall_S = SurfaceInternalMass(
    "IntWall",
    area=38.5,
    surface_type="IntWall",
    construction=int_wall_cs_S,
)
intdoor_S = SurfaceInternalMass(
    "IntWall",
    area=2.,
    surface_type="IntWall",
    construction=door_cs_S,
)
intceiling_S = SurfaceInternalMass(
    "IntCeiling",
    area=17.5,
    surface_type="IntCeiling",
    construction=ceiling_cs_S,
)
intfloor_S = SurfaceInternalMass(
    "IntFloor",
    area=17.5,
    surface_type="IntFloor",
    construction=floor_cs_S,
)
#########################################################

# Create zone
tz_S = ThermalZone(
    name="Zone S",
    surface_list=[wall_south_S, intwall_S, intdoor_S, intceiling_S, intfloor_S],
    net_floor_area=17.5,
    volume=52.5)

tz_L = ThermalZone(
    name="Zone L",
    surface_list=[wall_south_L, intwall_L, intdoor_L, intceiling_L, intfloor_L],
    net_floor_area=17.5,
    volume=52.5)

tz_S._ISO13790_params()
tz_S._VDI6007_params()
tz_L._ISO13790_params()
tz_L._VDI6007_params()



########################################################################################################################
#######################             ####################################################################################
#######################  TEST 1-3   ####################################################################################
#######################             ####################################################################################
########################################################################################################################

weather_file.hourly_data['out_air_db_temperature'] = np.array(
    [22.]*365*24
)
weather_file.general_data['average_dt_air_sky'] = 0.

weather_file.hourly_data_irradiances[0][90]['global'] = np.array([0.]*365*24)
weather_file.hourly_data_irradiances[0][90]['direct'] = np.array([0.]*365*24)
weather_file.hourly_data_irradiances[0][90]['AOI'] = np.array([0.]*365*24)

# Test 1 schedule
app_sched = Schedule(
    "Appliances",
    "percent",
    np.array(([0.] * 6 + [1.0] * 12 + [0.0] * 6) * 365),
)

# Loads
app_sens = ElectricLoad(
    name='AppliancesSens',
    unit='W',
    nominal_value=1000,
    schedule=app_sched,
    fraction_radiant=0.0,
    fraction_convective=1.0,
)

tz_S.add_internal_load(app_sens)
tz_L.add_internal_load(app_sens)

# IHG preprocessing
tz_S.extract_convective_radiative_latent_electric_load()
tz_S.calculate_zone_loads_ISO13790(weather_file)
tz_S.calculate_zone_loads_VDI6007(weather_file)
tz_L.extract_convective_radiative_latent_electric_load()
tz_L.calculate_zone_loads_ISO13790(weather_file)
tz_L.calculate_zone_loads_VDI6007(weather_file)


tz_L._air_thermal_capacity = 0.
tz_S._air_thermal_capacity = 0.

tz_S.theta_eq_tot = np.array(
    [22.]*365*24
)
tz_L.theta_eq_tot = np.array(
    [22.]*365*24
)

# Setpoints
heat_t = Schedule(
    "t_heat",
    "temperature",
    np.array([10.] * 365*24),
)

cool_t = Schedule(
    "t_cool",
    "temperature",
    np.array([80.] * 365*24),
)

heat_h = Schedule(
    "h_heat",
    "dimensionless",
    np.array([0.01] * 24 * 365),
)

cool_h = Schedule(
    "h_cool",
    "dimensionless",
    np.array([0.99] * 24 * 365),
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

tz_S.add_temperature_setpoint(t_sp)
tz_L.add_temperature_setpoint(t_sp)
tz_S.add_humidity_setpoint(h_sp)
tz_L.add_humidity_setpoint(h_sp)

# Ventilation

infiltration_sched = Schedule(
    "inf_sched",
    "dimensionless",
    np.array([0.0] * 24 * 365),
)

inf_obj = Infiltration(
    name='inf_obj',
    unit='Vol/h',
    nominal_value=1.,
    schedule=infiltration_sched,
)

# Natural Ventilation preprocessing
tz_S.add_infiltration(inf_obj)
tz_L.add_infiltration(inf_obj)
tz_L.calc_infiltration(weather_file)
tz_S.calc_infiltration(weather_file)




vent_sched = Schedule(
    "vent_sched",
    "dimensionless",
    np.array(([.0] * 8 + [0.] * 2 + [0] * 4 + [0.] * 10) * 365),
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
    np.array(([23.] * 8 + [23.] * 2 + [23.] * 4  + [23.] * 10 ) * 365),
)

x_supply_sched = Schedule(
    "x_supply_sched",
    "dimensionless",
    np.array(([0.0101] * 8 + [0.0101] * 2 + [0.0101] * 4 + [0.0101] * 10) * 365)*0.7,
)

availability_sched = np.array(([0] * 8 + [1] * 2  + [0] * 4 + [1] * 10) * 365)
availability_sched[120*24:273*24] = -1*availability_sched[120*24:273*24]
ahu_availability_sched = Schedule(
    "ahu_availability_sched",
    "availability",
    availability_sched,
)

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
    tz_S,
)

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
    tz_L,
)


# %% Simulation
preprocessing_timsteps = 0

if 1 in TESTS:
    df_S = pd.DataFrame(index=range(-preprocessing_timsteps, 8760),
                      columns=pd.MultiIndex.from_product(
                          [['1C', '2C'],
                           ['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'Ta_set_lower [°C]', 'Ta_set_upper [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                            'RH_set_lower [-]', 'RH_set_upper [-]', 'Sens Load [W]',
                            'Lat Load [W]']]))
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='1C')
        df_S.loc[t]['1C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"1C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='2C')
        df_S.loc[t]['2C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"2C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")

    df_S.to_csv(os.path.join('..', 'example_scripts', 'VDI_tests_results','Tests_1.csv'))

if 3 in TESTS:
    df_L = pd.DataFrame(index=range(-preprocessing_timsteps, 8760),
                      columns=pd.MultiIndex.from_product(
                          [['1C', '2C'],
                           ['Ta [°C]', 'Top [°C]', 'Tmr [°C]','Ta_set_lower [°C]', 'Ta_set_upper [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                            'RH_set_lower [-]', 'RH_set_upper [-]', 'Sens Load [W]',
                            'Lat Load [W]']]))
    start = time.time()
    tz_L.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_L.solve_timestep(t, weather_file, model='1C')
        df_L.loc[t]['1C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"1C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")
    start = time.time()
    tz_L.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_L.solve_timestep(t, weather_file, model='2C')
        df_L.loc[t]['2C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"2C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")

    df_L.to_csv(os.path.join('..', 'example_scripts', 'VDI_tests_results','Tests_3.csv'))

# fig, [[ax11, ax12], [ax21, ax22], [ax31, ax32]] = plt.subplots(ncols=2, nrows=3)
# for model, axes in zip(['1C', '2C'], [[ax11, ax21, ax31], [ax12, ax22, ax32]]):
#     ts_in = 1
#     ts_fine = 24
#     df_L.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Ta [°C]']].plot(ax=axes[0])
#     df_L.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['RH [-]']].plot(ax=axes[1],linestyle='-')
#     df_L.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Sens Load [W]', 'Lat Load [W]']].plot(ax=axes[2])
#
# fig, [[ax11, ax12], [ax21, ax22], [ax31, ax32]] = plt.subplots(ncols=2, nrows=3)
# for model, axes in zip(['1C', '2C'], [[ax11, ax21, ax31], [ax12, ax22, ax32]]):
#     ts_in = 1
#     ts_fine = 24
#     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Ta [°C]']].plot(ax=axes[0])
#     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['RH [-]']].plot(ax=axes[1],linestyle='-')
#     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Sens Load [W]', 'Lat Load [W]']].plot(ax=axes[2])



########################################################################################################################
#######################             ####################################################################################
#######################  TEST 2-4   ####################################################################################
#######################             ####################################################################################
########################################################################################################################

# Create zone
tz_S = ThermalZone(
    name="Zone S",
    surface_list=[wall_south_S, intwall_S, intdoor_S, intceiling_S, intfloor_S],
    net_floor_area=17.5,
    volume=52.5)

tz_L = ThermalZone(
    name="Zone L",
    surface_list=[wall_south_L, intwall_L, intdoor_L, intceiling_L, intfloor_L],
    net_floor_area=17.5,
    volume=52.5)

tz_S._ISO13790_params()
tz_S._VDI6007_params()
tz_L._ISO13790_params()
tz_L._VDI6007_params()

weather_file.hourly_data['out_air_db_temperature'] = np.array(
    [22.]*365*24
)
weather_file.general_data['average_dt_air_sky'] = 0.

weather_file.hourly_data_irradiances[0][90]['global'] = np.array([0]*365*24)
weather_file.hourly_data_irradiances[0][90]['direct'] = np.array([0]*365*24)
weather_file.hourly_data_irradiances[0][90]['AOI'] = np.array([0]*365*24)

# Test 1 schedule
app_sched = Schedule(
    "Appliances",
    "percent",
    np.array(([0.] * 6 + [1.0] * 12 + [0.0] * 6) * 365),
)

# Loads
app_sens = ElectricLoad(
    name='AppliancesSens',
    unit='W',
    nominal_value=1000,
    schedule=app_sched,
    fraction_radiant=1.,
    fraction_convective=0.0,
)

tz_S.add_internal_load(app_sens)
tz_L.add_internal_load(app_sens)

# IHG preprocessing
tz_S.extract_convective_radiative_latent_electric_load()
tz_S.calculate_zone_loads_ISO13790(weather_file)
tz_S.calculate_zone_loads_VDI6007(weather_file)
tz_L.extract_convective_radiative_latent_electric_load()
tz_L.calculate_zone_loads_ISO13790(weather_file)
tz_L.calculate_zone_loads_VDI6007(weather_file)


tz_L._air_thermal_capacity = 0.
tz_S._air_thermal_capacity = 0.

tz_S.theta_eq_tot = np.array(
    [22.]*365*24
)
tz_L.theta_eq_tot = np.array(
    [22.]*365*24
)

# Setpoints
heat_t = Schedule(
    "t_heat",
    "temperature",
    np.array([10] * 365*24),
)

cool_t = Schedule(
    "t_cool",
    "temperature",
    np.array([80] * 365*24),
)

heat_h = Schedule(
    "h_heat",
    "dimensionless",
    np.array([0.01] * 24 * 365),
)

cool_h = Schedule(
    "h_cool",
    "dimensionless",
    np.array([0.99] * 24 * 365),
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

tz_S.add_temperature_setpoint(t_sp)
tz_L.add_temperature_setpoint(t_sp)
tz_S.add_humidity_setpoint(h_sp)
tz_L.add_humidity_setpoint(h_sp)

# Ventilation

infiltration_sched = Schedule(
    "inf_sched",
    "dimensionless",
    np.array([0.0] * 24 * 365),
)

inf_obj = Infiltration(
    name='inf_obj',
    unit='Vol/h',
    nominal_value=1.,
    schedule=infiltration_sched,
)

# Natural Ventilation preprocessing
tz_S.add_infiltration(inf_obj)
tz_L.add_infiltration(inf_obj)
tz_L.calc_infiltration(weather_file)
tz_S.calc_infiltration(weather_file)




vent_sched = Schedule(
    "vent_sched",
    "dimensionless",
    np.array(([.0] * 8 + [0.] * 2 + [0] * 4 + [0.] * 10) * 365),
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
    np.array(([23.] * 8 + [23.] * 2 + [23.] * 4  + [23.] * 10 ) * 365),
)

x_supply_sched = Schedule(
    "x_supply_sched",
    "dimensionless",
    np.array(([0.0101] * 8 + [0.0101] * 2 + [0.0101] * 4 + [0.0101] * 10) * 365)*0.7,
)

availability_sched = np.array(([0] * 8 + [1] * 2  + [0] * 4 + [1] * 10) * 365)
availability_sched[120*24:273*24] = -1*availability_sched[120*24:273*24]
ahu_availability_sched = Schedule(
    "ahu_availability_sched",
    "availability",
    availability_sched,
)

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
    tz_S,
)

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
    tz_L,
)

# %% Simulation

preprocessing_timsteps = 0

if 2 in TESTS:
    df_S = pd.DataFrame(index=range(-preprocessing_timsteps, 8760),
                      columns=pd.MultiIndex.from_product(
                          [['1C', '2C'],
                           ['Ta [°C]', 'Top [°C]', 'Tmr [°C]','Ta_set_lower [°C]', 'Ta_set_upper [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                            'RH_set_lower [-]', 'RH_set_upper [-]', 'Sens Load [W]',
                            'Lat Load [W]']]))
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='1C')
        df_S.loc[t]['1C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"1C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='2C')
        df_S.loc[t]['2C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"2C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")

    df_S.to_csv(os.path.join('..', 'example_scripts', 'VDI_tests_results','Tests_2.csv'))

if 4 in TESTS:
    df_L = pd.DataFrame(index=range(-preprocessing_timsteps, 8760),
                      columns=pd.MultiIndex.from_product(
                          [['1C', '2C'],
                           ['Ta [°C]', 'Top [°C]', 'Tmr [°C]','Ta_set_lower [°C]', 'Ta_set_upper [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                            'RH_set_lower [-]', 'RH_set_upper [-]', 'Sens Load [W]',
                            'Lat Load [W]']]))
    start = time.time()
    tz_L.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_L.solve_timestep(t, weather_file, model='1C')
        df_L.loc[t]['1C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"1C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")
    start = time.time()
    tz_L.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_L.solve_timestep(t, weather_file, model='2C')
        df_L.loc[t]['2C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"2C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")

    df_L.to_csv(os.path.join('..', 'example_scripts', 'VDI_tests_results','Tests_4.csv'))

    # fig, [[ax11, ax12], [ax21, ax22], [ax31, ax32]] = plt.subplots(ncols=2, nrows=3)
    # for model, axes in zip(['1C', '2C'], [[ax11, ax21, ax31], [ax12, ax22, ax32]]):
    #     ts_in = 1
    #     ts_fine = 24
    #     df_L.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Ta [°C]']].plot(ax=axes[0])
    #     df_L.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['RH [-]']].plot(ax=axes[1],linestyle='-')
    #     df_L.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Sens Load [W]', 'Lat Load [W]']].plot(ax=axes[2])
    #
    # fig, [[ax11, ax12], [ax21, ax22], [ax31, ax32]] = plt.subplots(ncols=2, nrows=3)
    # for model, axes in zip(['1C', '2C'], [[ax11, ax21, ax31], [ax12, ax22, ax32]]):
    #     ts_in = 1
    #     ts_fine = 24
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Ta [°C]']].plot(ax=axes[0])
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['RH [-]']].plot(ax=axes[1],linestyle='-')
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Sens Load [W]', 'Lat Load [W]']].plot(ax=axes[2])

########################################################################################################################
#######################             ####################################################################################
#######################  TEST 5     ####################################################################################
#######################             ####################################################################################
########################################################################################################################

# Create zone
tz_S = ThermalZone(
    name="Zone S",
    surface_list=[wall_south_S, intwall_S, intdoor_S, intceiling_S, intfloor_S],
    net_floor_area=17.5,
    volume=52.5)

tz_S._ISO13790_params()
tz_S._VDI6007_params()

Out_T = np.array([18.8,17.1,16.5,16.1,16.5,17.8,
                  20.3,22.8,24.8,26.7,28.1,29,
                  29.7,30.4,30.9,31,30.8,30.1,
                  28.9,27,24.7,22.9,21.9,20.9,]*365)
Dir_irr = np.array([0,0,0,0,0,0,
                  0,18,13.05,25.8,35.1,38.4,
                  35.1,25.8,13.05,18,0,0,
                  0,0,0,0,0,0,]*365)
Diff_irr = np.array([0,0,0,0,17,38,
                     59,80,14.85,17.25,18.75,19.35,
                     18.75,17.25,14.85,80,59,38,
                     17,0,0,0,0,0,]*365)

weather_file.hourly_data['out_air_db_temperature'] = Out_T
weather_file.general_data['average_dt_air_sky'] = 0.

weather_file.hourly_data_irradiances[0][90]['global'] = Diff_irr + Dir_irr
weather_file.hourly_data_irradiances[0][90]['direct'] = Dir_irr
weather_file.hourly_data_irradiances[0][90]['AOI'] = np.array([0]*365*24)


tz_S._air_thermal_capacity = 0.

# Test 1 schedule
app_sched = Schedule(
    "Appliances",
    "percent",
    np.array(([0.] * 7 + [1.0] * 10 + [0.0] * 7) * 365),
)

# Loads
app_conv = ElectricLoad(
    name='AppliancesConv',
    unit='W',
    nominal_value=280,
    schedule=app_sched,
    fraction_radiant=0.0,
    fraction_convective=1.0,
)
app_rad = ElectricLoad(
    name='AppliancesRad',
    unit='W',
    nominal_value=80,
    schedule=app_sched,
    fraction_radiant=1.,
    fraction_convective=0.0,
)

tz_S.add_internal_load(app_conv)
tz_S.add_internal_load(app_rad)

# IHG preprocessing
tz_S.extract_convective_radiative_latent_electric_load()
tz_S.calculate_zone_loads_ISO13790(weather_file)
tz_S.calculate_zone_loads_VDI6007(weather_file)

# Setpoints
heat_t = Schedule(
    "t_heat",
    "temperature",
    np.array([10] * 365*24),
)

cool_t = Schedule(
    "t_cool",
    "temperature",
    np.array([80] * 365*24),
)

heat_h = Schedule(
    "h_heat",
    "dimensionless",
    np.array([0.01] * 24 * 365),
)

cool_h = Schedule(
    "h_cool",
    "dimensionless",
    np.array([0.99] * 24 * 365),
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

tz_S.add_temperature_setpoint(t_sp)
tz_S.add_humidity_setpoint(h_sp)

# Ventilation

infiltration_sched = Schedule(
    "inf_sched",
    "dimensionless",
    np.array([0.0] * 24 * 365),
)

inf_obj = Infiltration(
    name='inf_obj',
    unit='Vol/h',
    nominal_value=1.,
    schedule=infiltration_sched,
)

# Natural Ventilation preprocessing
tz_S.add_infiltration(inf_obj)
tz_S.calc_infiltration(weather_file)




vent_sched = Schedule(
    "vent_sched",
    "dimensionless",
    np.array(([.0] * 8 + [0.] * 2 + [0] * 4 + [0.] * 10) * 365),
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
    np.array(([23.] * 8 + [23.] * 2 + [23.] * 4  + [23.] * 10 ) * 365),
)

x_supply_sched = Schedule(
    "x_supply_sched",
    "dimensionless",
    np.array(([0.0101] * 8 + [0.0101] * 2 + [0.0101] * 4 + [0.0101] * 10) * 365)*0.7,
)

availability_sched = np.array(([0] * 8 + [1] * 2  + [0] * 4 + [1] * 10) * 365)
availability_sched[120*24:273*24] = -1*availability_sched[120*24:273*24]
ahu_availability_sched = Schedule(
    "ahu_availability_sched",
    "availability",
    availability_sched,
)

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
    tz_S,
)

# %% Simulation

preprocessing_timsteps = 0

if 5 in TESTS:
    df_S = pd.DataFrame(index=range(-preprocessing_timsteps, 8760),
                      columns=pd.MultiIndex.from_product(
                          [['1C', '2C'],
                           ['Ta [°C]', 'Top [°C]', 'Tmr [°C]','Ta_set_lower [°C]', 'Ta_set_upper [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                            'RH_set_lower [-]', 'RH_set_upper [-]', 'Sens Load [W]',
                            'Lat Load [W]']]))
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='1C')
        df_S.loc[t]['1C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"1C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='2C')
        df_S.loc[t]['2C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"2C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")

    df_S.to_csv(os.path.join('..', 'example_scripts', 'VDI_tests_results','Tests_5.csv'))

    # fig, [[ax11, ax12], [ax21, ax22], [ax31, ax32]] = plt.subplots(ncols=2, nrows=3)
    # for model, axes in zip(['1C', '2C'], [[ax11, ax21, ax31], [ax12, ax22, ax32]]):
    #     ts_in = 1
    #     ts_fine = 24
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Ta [°C]']].plot(ax=axes[0])
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['RH [-]']].plot(ax=axes[1],linestyle='-')
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Sens Load [W]', 'Lat Load [W]']].plot(ax=axes[2])

########################################################################################################################
#######################             ####################################################################################
#######################  TEST 6     ####################################################################################
#######################             ####################################################################################
########################################################################################################################

# Create zone
tz_S = ThermalZone(
    name="Zone S",
    surface_list=[wall_south_S, intwall_S, intdoor_S, intceiling_S, intfloor_S],
    net_floor_area=17.5,
    volume=52.5)

tz_S._ISO13790_params()
tz_S._VDI6007_params()

weather_file.hourly_data['out_air_db_temperature'] = np.array(
    [22.]*365*24
)
weather_file.general_data['average_dt_air_sky'] = 0.

weather_file.hourly_data_irradiances[0][90]['global'] = np.array([0]*365*24)
weather_file.hourly_data_irradiances[0][90]['direct'] = np.array([0]*365*24)
weather_file.hourly_data_irradiances[0][90]['AOI'] = np.array([0]*365*24)

# Test 1 schedule
app_sched = Schedule(
    "Appliances",
    "percent",
    np.array(([0.] * 6 + [1.0] * 12 + [0.0] * 6) * 365),
)

# Loads
app_sens = ElectricLoad(
    name='AppliancesSens',
    unit='W',
    nominal_value=1000,
    schedule=app_sched,
    fraction_radiant=1.,
    fraction_convective=0.0,
)

tz_S.add_internal_load(app_sens)

# IHG preprocessing
tz_S.extract_convective_radiative_latent_electric_load()
tz_S.calculate_zone_loads_ISO13790(weather_file)
tz_S.calculate_zone_loads_VDI6007(weather_file)

tz_S._air_thermal_capacity = 0.


tz_S.theta_eq_tot = np.array(
    [22.]*365*24
)

# Setpoints
sp_sched = np.array(([22.] * 6 + [27.] * 12 + [22.] * 6) * 365)
heat_t = Schedule(
    "t_heat",
    "temperature",
    sp_sched,
)

cool_t = Schedule(
    "t_cool",
    "temperature",
    sp_sched + 0.001,
)

heat_h = Schedule(
    "h_heat",
    "dimensionless",
    np.array([0.01] * 24 * 365),
)

cool_h = Schedule(
    "h_cool",
    "dimensionless",
    np.array([0.99] * 24 * 365),
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

tz_S.add_temperature_setpoint(t_sp)
tz_S.add_humidity_setpoint(h_sp)

# Ventilation

infiltration_sched = Schedule(
    "inf_sched",
    "dimensionless",
    np.array([0.0] * 24 * 365),
)

inf_obj = Infiltration(
    name='inf_obj',
    unit='Vol/h',
    nominal_value=1.,
    schedule=infiltration_sched,
)

# Natural Ventilation preprocessing
tz_S.add_infiltration(inf_obj)
tz_S.calc_infiltration(weather_file)




vent_sched = Schedule(
    "vent_sched",
    "dimensionless",
    np.array(([.0] * 8 + [0.] * 2 + [0] * 4 + [0.] * 10) * 365),
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
    np.array(([23.] * 8 + [23.] * 2 + [23.] * 4  + [23.] * 10 ) * 365),
)

x_supply_sched = Schedule(
    "x_supply_sched",
    "dimensionless",
    np.array(([0.0101] * 8 + [0.0101] * 2 + [0.0101] * 4 + [0.0101] * 10) * 365)*0.7,
)

availability_sched = np.array(([0] * 8 + [1] * 2  + [0] * 4 + [1] * 10) * 365)
availability_sched[120*24:273*24] = -1*availability_sched[120*24:273*24]
ahu_availability_sched = Schedule(
    "ahu_availability_sched",
    "availability",
    availability_sched,
)

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
    tz_S,
)

# %% Simulation

preprocessing_timsteps = 0

if 6 in TESTS:
    df_S = pd.DataFrame(index=range(-preprocessing_timsteps, 8760),
                      columns=pd.MultiIndex.from_product(
                          [['1C', '2C'],
                           ['Ta [°C]', 'Top [°C]', 'Tmr [°C]','Ta_set_lower [°C]', 'Ta_set_upper [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                            'RH_set_lower [-]', 'RH_set_upper [-]', 'Sens Load [W]',
                            'Lat Load [W]']]))
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='1C')
        df_S.loc[t]['1C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"1C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='2C')
        df_S.loc[t]['2C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"2C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")

    df_S.to_csv(os.path.join('..', 'example_scripts', 'VDI_tests_results','Tests_6.csv'))

    # fig, [[ax11, ax12], [ax21, ax22], [ax31, ax32]] = plt.subplots(ncols=2, nrows=3)
    # for model, axes in zip(['1C', '2C'], [[ax11, ax21, ax31], [ax12, ax22, ax32]]):
    #     ts_in = 1
    #     ts_fine = 24
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Ta [°C]']].plot(ax=axes[0])
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['RH [-]']].plot(ax=axes[1],linestyle='-')
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Sens Load [W]', 'Lat Load [W]']].plot(ax=axes[2])

########################################################################################################################
#######################             ####################################################################################
#######################  TEST 7     ####################################################################################
#######################             ####################################################################################
########################################################################################################################

# Create zone
tz_S = ThermalZone(
    name="Zone S",
    surface_list=[wall_south_S, intwall_S, intdoor_S, intceiling_S, intfloor_S],
    net_floor_area=17.5,
    volume=52.5)

tz_S.design_heating_system_power = 500
tz_S.design_cooling_system_power = -500

tz_S._ISO13790_params()
tz_S._VDI6007_params()

weather_file.hourly_data['out_air_db_temperature'] = np.array(
    [22.]*365*24
)
weather_file.general_data['average_dt_air_sky'] = 0.

weather_file.hourly_data_irradiances[0][90]['global'] = np.array([0]*365*24)
weather_file.hourly_data_irradiances[0][90]['direct'] = np.array([0]*365*24)
weather_file.hourly_data_irradiances[0][90]['AOI'] = np.array([0]*365*24)

# Test 1 schedule
app_sched = Schedule(
    "Appliances",
    "percent",
    np.array(([0.] * 6 + [1.0] * 12 + [0.0] * 6) * 365),
)

# Loads
app_sens = ElectricLoad(
    name='AppliancesSens',
    unit='W',
    nominal_value=1000,
    schedule=app_sched,
    fraction_radiant=1.,
    fraction_convective=0.0,
)

tz_S.add_internal_load(app_sens)

# IHG preprocessing
tz_S.extract_convective_radiative_latent_electric_load()
tz_S.calculate_zone_loads_ISO13790(weather_file)
tz_S.calculate_zone_loads_VDI6007(weather_file)

tz_S._air_thermal_capacity = 0.


tz_S.theta_eq_tot = np.array(
    [22.]*365*24
)

# Setpoints
sp_sched = np.array(([22.] * 6 + [27.] * 12 + [22.] * 6) * 365)
heat_t = Schedule(
    "t_heat",
    "temperature",
    sp_sched,
)

cool_t = Schedule(
    "t_cool",
    "temperature",
    sp_sched + 0.001,
)

heat_h = Schedule(
    "h_heat",
    "dimensionless",
    np.array([0.01] * 24 * 365),
)

cool_h = Schedule(
    "h_cool",
    "dimensionless",
    np.array([0.99] * 24 * 365),
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

tz_S.add_temperature_setpoint(t_sp)
tz_S.add_humidity_setpoint(h_sp)

# Ventilation

infiltration_sched = Schedule(
    "inf_sched",
    "dimensionless",
    np.array([0.0] * 24 * 365),
)

inf_obj = Infiltration(
    name='inf_obj',
    unit='Vol/h',
    nominal_value=1.,
    schedule=infiltration_sched,
)

# Natural Ventilation preprocessing
tz_S.add_infiltration(inf_obj)
tz_S.calc_infiltration(weather_file)




vent_sched = Schedule(
    "vent_sched",
    "dimensionless",
    np.array(([.0] * 8 + [0.] * 2 + [0] * 4 + [0.] * 10) * 365),
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
    np.array(([23.] * 8 + [23.] * 2 + [23.] * 4  + [23.] * 10 ) * 365),
)

x_supply_sched = Schedule(
    "x_supply_sched",
    "dimensionless",
    np.array(([0.0101] * 8 + [0.0101] * 2 + [0.0101] * 4 + [0.0101] * 10) * 365)*0.7,
)

availability_sched = np.array(([0] * 8 + [1] * 2  + [0] * 4 + [1] * 10) * 365)
availability_sched[120*24:273*24] = -1*availability_sched[120*24:273*24]
ahu_availability_sched = Schedule(
    "ahu_availability_sched",
    "availability",
    availability_sched,
)

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
    tz_S,
)

# %% Simulation

preprocessing_timsteps = 0

if 7 in TESTS:
    df_S = pd.DataFrame(index=range(-preprocessing_timsteps, 8760),
                      columns=pd.MultiIndex.from_product(
                          [['1C', '2C'],
                           ['Ta [°C]', 'Top [°C]', 'Tmr [°C]','Ta_set_lower [°C]', 'Ta_set_upper [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                            'RH_set_lower [-]', 'RH_set_upper [-]', 'Sens Load [W]',
                            'Lat Load [W]']]))
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='1C')
        df_S.loc[t]['1C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"1C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='2C')
        df_S.loc[t]['2C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"2C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")

    df_S.to_csv(os.path.join('..', 'example_scripts', 'VDI_tests_results','Tests_7.csv'))

    # fig, [[ax11, ax12], [ax21, ax22], [ax31, ax32]] = plt.subplots(ncols=2, nrows=3)
    # for model, axes in zip(['1C', '2C'], [[ax11, ax21, ax31], [ax12, ax22, ax32]]):
    #     ts_in = 1
    #     ts_fine = 24
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Ta [°C]']].plot(ax=axes[0])
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['RH [-]']].plot(ax=axes[1],linestyle='-')
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Sens Load [W]', 'Lat Load [W]']].plot(ax=axes[2])


########################################################################################################################
#######################             ####################################################################################
#######################  TEST 12    ####################################################################################
#######################             ####################################################################################
########################################################################################################################

# Create zone
tz_S = ThermalZone(
    name="Zone S",
    surface_list=[wall_south_S, intwall_S, intdoor_S, intceiling_S, intfloor_S],
    net_floor_area=17.5,
    volume=52.5)

tz_S.design_heating_system_power = 500
tz_S.design_cooling_system_power = -500

tz_S._ISO13790_params()
tz_S._VDI6007_params()

Out_T = np.array([18.8,17.1,16.5,16.1,16.5,17.8,
                  20.3,22.8,24.8,26.7,28.1,29,
                  29.7,30.4,30.9,31,30.8,30.1,
                  28.9,27,24.7,22.9,21.9,20.9,]*365)
Dir_irr = np.array([0,0,0,0,0,0,
                  0,18,13.05,25.8,35.1,38.4,
                  35.1,25.8,13.05,18,0,0,
                  0,0,0,0,0,0,]*365)
Diff_irr = np.array([0,0,0,0,17,38,
                     59,80,14.85,17.25,18.75,19.35,
                     18.75,17.25,14.85,80,59,38,
                     17,0,0,0,0,0,]*365)


weather_file.hourly_data['out_air_db_temperature'] = Out_T
weather_file.general_data['average_dt_air_sky'] = 0.

weather_file.hourly_data_irradiances[0][90]['global'] = Diff_irr + Dir_irr
weather_file.hourly_data_irradiances[0][90]['direct'] = Dir_irr
weather_file.hourly_data_irradiances[0][90]['AOI'] = np.array([0]*365*24)


# Test 1 schedule
app_sched = Schedule(
    "Appliances",
    "percent",
    np.array(([0.] * 7 + [1.0] * 10 + [0.0] * 7) * 365),
)

# Loads
app_conv = ElectricLoad(
    name='AppliancesConv',
    unit='W',
    nominal_value=280,
    schedule=app_sched,
    fraction_radiant=0.0,
    fraction_convective=1.0,
)
app_rad = ElectricLoad(
    name='AppliancesRad',
    unit='W',
    nominal_value=80,
    schedule=app_sched,
    fraction_radiant=1.,
    fraction_convective=0.0,
)

tz_S.add_internal_load(app_conv)
tz_S.add_internal_load(app_rad)

# IHG preprocessing
tz_S.extract_convective_radiative_latent_electric_load()
tz_S.calculate_zone_loads_ISO13790(weather_file)
tz_S.calculate_zone_loads_VDI6007(weather_file)

tz_S._air_thermal_capacity = 0.

tz_S.theta_eq_tot = np.array(
    [22.]*365*24
)

# Setpoints
heat_t = Schedule(
    "t_heat",
    "temperature",
    np.array(([10] * 24) * 365),
)

cool_t = Schedule(
    "t_cool",
    "temperature",
    np.array(([80] * 24) * 365),
)

heat_h = Schedule(
    "h_heat",
    "dimensionless",
    np.array([0.01] * 24 * 365),
)

cool_h = Schedule(
    "h_cool",
    "dimensionless",
    np.array([0.99] * 24 * 365),
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

tz_S.add_temperature_setpoint(t_sp)
tz_S.add_humidity_setpoint(h_sp)

# Ventilation

infiltration_sched = Schedule(
    "inf_sched",
    "dimensionless",
    np.array(([1.] * 7 + [0.5] * 10 + [1.] * 7) * 365),
)

inf_obj = Infiltration(
    name='inf_obj',
    unit="m3/s",
    nominal_value=100./3600,
    schedule=infiltration_sched,
)

# Natural Ventilation preprocessing
tz_S.add_infiltration(inf_obj)
tz_S.calc_infiltration(weather_file)



vent_sched = Schedule(
    "vent_sched",
    "dimensionless",
    np.array(([.0] * 8 + [0.] * 2 + [0] * 4 + [0.] * 10) * 365),
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
    np.array(([23.] * 8 + [23.] * 2 + [23.] * 4  + [23.] * 10 ) * 365),
)

x_supply_sched = Schedule(
    "x_supply_sched",
    "dimensionless",
    np.array(([0.0101] * 8 + [0.0101] * 2 + [0.0101] * 4 + [0.0101] * 10) * 365)*0.7,
)

availability_sched = np.array(([0] * 8 + [1] * 2  + [0] * 4 + [1] * 10) * 365)
availability_sched[120*24:273*24] = -1*availability_sched[120*24:273*24]
ahu_availability_sched = Schedule(
    "ahu_availability_sched",
    "availability",
    availability_sched,
)

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
    tz_S,
)
# %% Simulation

preprocessing_timsteps = 0

if 12 in TESTS:
    df_S = pd.DataFrame(index=range(-preprocessing_timsteps, 8760),
                      columns=pd.MultiIndex.from_product(
                          [['1C', '2C'],
                           ['Ta [°C]', 'Top [°C]', 'Tmr [°C]','Ta_set_lower [°C]', 'Ta_set_upper [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                            'RH_set_lower [-]', 'RH_set_upper [-]', 'Sens Load [W]',
                            'Lat Load [W]']]))
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='1C')
        df_S.loc[t]['1C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"1C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")
    start = time.time()
    tz_S.reset_init_values_VDI()
    for t in range(-preprocessing_timsteps, 8760):
        [ta,top,tmr], rha, [tmaw, tmiw], pot, lat_pot = tz_S.solve_timestep(t, weather_file, model='2C')
        df_S.loc[t]['2C'][['Ta [°C]', 'Top [°C]', 'Tmr [°C]', 'TmAW [°C]', 'TmIW [°C]', 'RH [-]',
                        'Sens Load [W]',
                        'Lat Load [W]']] = [ta,top,tmr, tmaw, tmiw, rha, pot, lat_pot]
    print(f"2C model: \n\t{8760} time steps\n\t{(time.time() - start):.2f} s")

    df_S.to_csv(os.path.join('..', 'example_scripts', 'VDI_tests_results','Tests_12.csv'))

    # fig, [[ax11, ax12], [ax21, ax22], [ax31, ax32]] = plt.subplots(ncols=2, nrows=3)
    # for model, axes in zip(['1C', '2C'], [[ax11, ax21, ax31], [ax12, ax22, ax32]]):
    #     ts_in = 1
    #     ts_fine = 24
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Ta [°C]']].plot(ax=axes[0])
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['RH [-]']].plot(ax=axes[1],linestyle='-')
    #     df_S.iloc[preprocessing_timsteps + ts_in: ts_fine][model][['Sens Load [W]', 'Lat Load [W]']].plot(ax=axes[2])