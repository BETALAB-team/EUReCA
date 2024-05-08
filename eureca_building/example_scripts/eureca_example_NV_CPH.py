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

ext_wall = Construction.from_U_value("ExtWall", 0.7, weight_class="Medium", construction_type="ExtWall")    # Modificare con i dati di trasmittanza per le strutture
ext_wall_from_U = Construction.from_U_value("ExtWall from U", 0.7,weight_class="Medium", construction_type="ExtWall")
ext_wall_from_U = Construction.from_U_value("ExtWall from U", 0.7,weight_class="Medium", construction_type="ExtWall")

window = SimpleWindow(
                name="Window",
                u_value=1.6,
                solar_heat_gain_coef=0.4,
                visible_transmittance=0.9,
                frame_factor=0.10, # TODO: To be verified
                shading_coef_int=0.46,
                shading_coef_ext=0.99,
)
#########################################################

# Definition of surfaces
wall_south = Surface(
    "Wall 1",
    vertices=((0, 0, 0), (21.36, 0, 0), (21.36, 0, 6), (0, 0, 6)),
    wwr=0.125,
    surface_type="ExtWall",
    construction=ext_wall_from_U,
    window=window,
    n_window_layers=1
)
intwall = SurfaceInternalMass(
    "IntWall",
    area=floor._area*2.5*2,  # Area delle pareti interne
    surface_type="IntWall",
    construction=int_wall_cs,
)


#########################################################
# Loads

ts_h = CONFIG.ts_per_hour
delay_ts = 8760*ts_h+1-ts_h

# A schedule
people_sched = Schedule(
    "PeopleOccupancy1",
    "percent",
    np.array(([0.1] * 7 * ts_h + [0.6] * 2 * ts_h + [0.4] * 5 * ts_h + [0.6] * 10 * ts_h) * 365)[:delay_ts],
)

app_sched = Schedule(
    "AppliancesSchedule",
    "dimensionless",
    np.array(([0.1] * 7 * ts_h + [0.6] * 2 * ts_h + [0.4] * 5 * ts_h + [0.6] * 10 * ts_h) * 365)[:delay_ts],
)

# Loads
people = People(
    name='occupancy_tz',
    unit='px',
    nominal_value=1,
    schedule=people_sched,
    fraction_latent=0.45,
    fraction_radiant=0.3,
    fraction_convective=0.7,
    metabolic_rate=150,
)

electric_devices = ElectricLoad(
    name='ElectricLoads',
    unit='W',
    nominal_value=300.,  # Un valore indicativo di 1.5 kW
    schedule=people_sched,
    fraction_radiant=0.3,
    fraction_convective=0.7,
    fraction_to_zone=1
)

#########################################################
# Setpoints

heat_t = Schedule.from_constant_value(
    name = "t_heat",
    schedule_type="temperature",
    value = 0.
)

cool_t = Schedule.from_constant_value(
    name = "t_cool",
    schedule_type="temperature",
    value = 50.
)

heat_h = Schedule.from_constant_value(
    name = "h_heat",
    schedule_type = "dimensionless",
    value = 0
)

cool_h = Schedule.from_constant_value(
    name = "h_cool",
    schedule_type = "dimensionless",
    value = 0.99
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

natural_vent_sched = Schedule(
    "nat_vent_sche",
    "dimensionless",
    np.array(([.3] * 8 * ts_h + [.5] * 2 * ts_h + [.3] * 4 * ts_h + [.5] * 10 * ts_h) * 365)[:delay_ts],
)

nv_obj = NaturalVentilation(
    name='nat_vent',
    unit='%',
    nominal_value=99, # 90% moltiplicato per il vettore della schedule
    schedule=natural_vent_sched,
)

#########################################################
# Mechanical ventilation
vent_sched = Schedule.from_constant_value(
    name = "vent_sched",
    schedule_type = "dimensionless",
    value = 1,
)

mech_vent_obj = MechanicalVentilation(
    name='vent_obj',
    unit='m3/s',
    nominal_value=1,   # Portata di ventilazione in m3/s
    schedule=vent_sched,
)

T_supply_sched = Schedule(
    "T_supply_sched",
    "temperature",
    np.array(([23.] * 8 * ts_h + [23.] * 2 * ts_h + [23.] * 4 * ts_h + [23.] * 10 * ts_h) * 365)[:delay_ts],
)

x_supply_sched = Schedule(
    "x_supply_sched",
    "dimensionless",
    np.array(([0.0101] * 8 * ts_h + [0.0101] * 2 * ts_h + [0.0101] * 4 * ts_h + [0.0101] * 10 * ts_h) * 365)[:delay_ts]*0.7,
)

ahu_availability_sched = Schedule.from_constant_value(
    name = "ahu_availability_sched",
    schedule_type = "availability",
    value = 0, # 0 equivale a free-cooling
)


# DHW

dhw_flow_rate = Schedule.from_constant_value(
    name = "dhw_flow_rate",
    schedule_type = "mass_flow_rate",
    value = 0,
)

dhw_1 = DomesticHotWater(
    "dhw_1",
    calculation_method="Schedule",
    unit = "L/(m2 h)",
    schedule=dhw_flow_rate,

)

#########################################################

# Create zone
tz1 = ThermalZone(
    name="Zone 1",
    surface_list=[wall_south, intwall],
    net_floor_area=78,
    volume=78*2.7)

tz1._ISO13790_params()
# tz1._VDI6007_params()

tz1.add_internal_load(people,electric_devices)

# IHG preprocessing
tz_loads = tz1.extract_convective_radiative_latent_electric_load()
tz1.calculate_zone_loads_ISO13790(weather_file)


# 2C model
# tz1.calculate_zone_loads_VDI6007(weather_file)


tz1.add_temperature_setpoint(t_sp)
tz1.add_humidity_setpoint(h_sp)

# Natural Ventilation preprocessing
tz1.add_natural_ventilation(nv_obj, weather_file)

ahu = AirHandlingUnit(
"ahu",
mech_vent_obj,
T_supply_sched,
x_supply_sched,
ahu_availability_sched,
True,
0,
0,
1,
weather_file,
tz1,
)

cooling_1C_peak_load = tz1.design_sensible_cooling_load(weather_file, model = "1C")
heating_peak_load = tz1.design_heating_load(-5.)

tz1.add_domestic_hot_water(weather_file, dhw_1)

bd = Building("Bd 1", thermal_zones_list=[tz1], model = "1C")
bd.set_hvac_system("Traditional Gas Boiler, Centralized, Low Temp Radiator", "A-W chiller, Centralized, Radiant surface")
bd.set_hvac_system_capacity(weather_file)
start = time.time()
df_res = bd.simulate(weather_file, output_folder="Results")
print(f"1C model: \n\t{8760 * 2 - 1} time steps\n\t{(time.time() - start):.2f} s")


hourly = pd.read_csv(os.path.join('Results','Results Bd 1.csv'), header = [0,1], delimiter = ";")
hourly.set_index(pd.date_range(start=CONFIG.start_date,end=CONFIG.final_date,periods=CONFIG.number_of_time_steps), inplace = True)
hourly.columns = hourly.columns.droplevel(level = 1)
hourly_res = hourly[['TZ sensible load [W]','TZ latent load [W]','TZ AHU pre heater load [W]','TZ AHU post heater load [W]']]
hourly_res['AHU'] = hourly_res[['TZ AHU pre heater load [W]','TZ AHU post heater load [W]']].sum(axis=1)
hourly_res.drop(['TZ AHU pre heater load [W]','TZ AHU post heater load [W]'],axis = 1,inplace = True)
hourly_res = hourly_res.resample('1M').sum()/1000

hourly_res.plot(kind = 'bar')

fig, [ax1, ax2,ax3] = plt.subplots(nrows = 3)

for ax, i in zip ([ax1, ax2,ax3],range(3)):
    pd.DataFrame({
        hourly_res.columns[i] : hourly_res[hourly_res.columns[i]].values,
    }).plot(ax = ax, kind = 'bar')