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

config_path = os.path.join('config.json')
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

#########################################################
# Epw loading
epw_path = os.path.join('../eureca_building', 'example_scripts', 'ITA_Venezia-Tessera.161050_IGDG.epw')
weather_file = WeatherFile(epw_path,
                           time_steps=CONFIG.ts_per_hour,
                           azimuth_subdivisions=CONFIG.azimuth_subdivisions,
                           height_subdivisions=CONFIG.height_subdivisions, )
#########################################################
path = os.path.join(
    "../eureca_building",
    "example_scripts",
    "materials_and_construction_test.xlsx",
)
# Define some constructions

dataset = ConstructionDataset.read_excel(path)
roof_cs = dataset.constructions_dict[13]
ceiling_cs = dataset.constructions_dict[17]
floor_cs = dataset.constructions_dict[16]
ext_wall_cs = dataset.constructions_dict[70]
int_wall_cs = dataset.constructions_dict[70]
window_cs = dataset.windows_dict[2]
mat_cs = dataset.materials_dict[1]

ext_wall_from_U = Construction.from_U_value("ExtWall from U", 0.7, weight_class="Medium", construction_type="ExtWall")
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
    area=floor._area * 2.5 * 2,
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

ts_h = CONFIG.ts_per_hour
delay_ts = 8760 * ts_h + 1 - ts_h

# A schedule
people_sched = Schedule(
    "PeopleOccupancy1",
    "percent",
    np.array(([0.1] * 7 * ts_h + [0.6] * 2 * ts_h + [0.4] * 5 * ts_h + [0.6] * 10 * ts_h) * 365)[:delay_ts],
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
    name="t_heat",
    schedule_type="temperature",
    schedule_week_day=np.array([18] * 7 * ts_h + [21] * 2 * ts_h + [18] * 5 * ts_h + [21] * 10 * ts_h),
    schedule_saturday=np.array([18] * 7 * ts_h + [21] * 2 * ts_h + [18] * 5 * ts_h + [21] * 10 * ts_h) - 5,
    schedule_sunday=np.array([18] * 7 * ts_h + [21] * 2 * ts_h + [18] * 5 * ts_h + [21] * 10 * ts_h) - 10,
    schedule_holiday=np.array([18] * 7 * ts_h + [21] * 2 * ts_h + [18] * 5 * ts_h + [21] * 10 * ts_h) * 0,
    holidays=(10, 11, 12, 13, 14),
    starting_day=3,
)

heat_t = Schedule.from_constant_value(
    name="t_heat",
    schedule_type="temperature",
    value=15.
)

heat_t = Schedule(
    "t_heat",
    "temperature",
    np.array(([18] * 7 * ts_h + [21] * 2 * ts_h + [18] * 5 * ts_h + [21] * 10 * ts_h) * 365)[:delay_ts],
)

cool_t = Schedule(
    "t_cool",
    "temperature",
    np.array(([28] * 8 * ts_h + [26] * 2 * ts_h + [28] * 4 * ts_h + [26] * 10 * ts_h) * 365)[:delay_ts],
)

heat_h = Schedule(
    "h_heat",
    "dimensionless",
    np.array(([0.1] * 7 * ts_h + [0.3] * 2 * ts_h + [.1] * 5 * ts_h + [.3] * 10 * ts_h) * 365)[:delay_ts],
)

cool_h = Schedule(
    "h_cool",
    "dimensionless",
    np.array(([.9] * 8 * ts_h + [.5] * 2 * ts_h + [.9] * 4 * ts_h + [.5] * 10 * ts_h) * 365)[:delay_ts],
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
    np.array(([.3] * 8 * ts_h + [.5] * 2 * ts_h + [.3] * 4 * ts_h + [.5] * 10 * ts_h) * 365)[:delay_ts],
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
    np.array(([.0] * 8 * ts_h + [1] * 2 * ts_h + [0] * 4 * ts_h + [1] * 10 * ts_h) * 365)[:delay_ts],
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
    np.array(([23.] * 8 * ts_h + [23.] * 2 * ts_h + [23.] * 4 * ts_h + [23.] * 10 * ts_h) * 365)[:delay_ts],
)

x_supply_sched = Schedule(
    "x_supply_sched",
    "dimensionless",
    np.array(([0.0101] * 8 * ts_h + [0.0101] * 2 * ts_h + [0.0101] * 4 * ts_h + [0.0101] * 10 * ts_h) * 365)[
    :delay_ts] * 0.7,
)

availability_sched = np.array(([0] * 8 * ts_h + [1] * 2 * ts_h + [0] * 4 * ts_h + [1] * 10 * ts_h) * 365)[:delay_ts]
availability_sched[120 * 24 * ts_h:273 * 24 * ts_h] = -1 * availability_sched[120 * 24 * ts_h:273 * 24 * ts_h]
ahu_availability_sched = Schedule(
    "ahu_availability_sched",
    "availability",
    availability_sched,
)

# DHW

dhw_flow_rate = Schedule(
    "dhw_flow_rate",
    "mass_flow_rate",
    np.array(([.5] * 8 * ts_h + [.0] * 2 * ts_h + [.5] * 4 * ts_h + [.0] * 10 * ts_h) * 365)[:delay_ts] * 10,
)

dhw_1 = DomesticHotWater(
    "dhw_1",
    calculation_method="Schedule",
    unit="L/(m2 h)",
    schedule=dhw_flow_rate,

)

dhw_2 = DomesticHotWater(
    "dhw_2",
    calculation_method="UNI-TS 11300-2",
)

#########################################################

zones = []

for i in range(1):
    # Create zone
    tz1 = ThermalZone(
        name="Zone 1",
        surface_list=[wall_south, wall_east, wall_north, wall_west, roof, floor, intwall, intceiling],
        net_floor_area=floor._area * 2,
        volume=floor._area * 3.3 * 2)

    zones.append(tz1)
    tz1._ISO13790_params()
    tz1._VDI6007_params()

    tz1.add_internal_load(people)
    tz1.add_internal_load(lights, pc)

    # IHG preprocessing
    tz_loads = tz1.extract_convective_radiative_latent_electric_load()
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

    cooling_1C_peak_load = tz1.design_sensible_cooling_load(weather_file, model="1C")
    heating_peak_load = tz1.design_heating_load(-5.)

    tz1.add_domestic_hot_water(weather_file, dhw_1, dhw_2)

    bd = Building("Bd 1", thermal_zones_list=[tz1], model="1C")
    bd.set_hvac_system("Traditional Gas Boiler, Centralized, Low Temp Radiator",
                       "A-W chiller, Centralized, Radiant surface")
    bd.set_hvac_system_capacity(weather_file)

    tz2 = ThermalZone(
        name="Zone 2",
        surface_list=[wall_north, wall_west, roof, floor, intwall, intceiling],
        net_floor_area=floor._area * 2,
        volume=floor._area * 3.3 * 2)

    zones.append(tz1)
    tz2._ISO13790_params()
    tz2._VDI6007_params()

    tz2.add_internal_load(people)
    tz2.add_internal_load(lights, pc)

    # IHG preprocessing
    tz_loads = tz2.extract_convective_radiative_latent_electric_load()
    tz2.calculate_zone_loads_ISO13790(weather_file)
    # tz1._plot_ISO13790_IHG()

    # 2C model
    tz2.calculate_zone_loads_VDI6007(weather_file)
    # tz1._plot_VDI6007_IHG(weather_file)

    tz2.add_temperature_setpoint(t_sp)
    tz2.add_humidity_setpoint(h_sp)

    # Natural Ventilation preprocessing
    tz2.add_infiltration(inf_obj)
    tz_inf = tz2.calc_infiltration(weather_file)

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
        tz2,
    )

    cooling_1C_peak_load = tz2.design_sensible_cooling_load(weather_file, model="1C")
    heating_peak_load = tz2.design_heating_load(-5.)

    tz2.add_domestic_hot_water(weather_file, dhw_1, dhw_2)

    bd2 = Building("Bd 2", thermal_zones_list=[tz2], model="1C")
    bd2.set_hvac_system("Traditional Gas Boiler, Centralized, Low Temp Radiator",
                       "A-W chiller, Centralized, Radiant surface")
    bd2.set_hvac_system_capacity(weather_file)
    # tz1.print_modelica_class_with_output_csv("C:\\Users\\pratenr82256\\Desktop\\ModelicaEUReCA\\ClassiDaEUReCA")


    weather_mos = os.path.join(os.getcwd(), "ExampleFiles", "ITA_Venezia-Tessera.161050_IGDG.mos")
    path_example = os.path.join(os.getcwd(), "ExampleFiles")

    mod_str1 = bd.get_modelica_bd_model(
        weather_mos,path_example,x = 50, y = 50)

    mod_str2 = bd2.get_modelica_bd_model(
        weather_mos,path_example,x = 50, y = -70)


    modelica_model = f"""within ;
model EurecaModelicaModel
{mod_str1}
{mod_str2}
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=31536000,
      Interval={CONFIG.time_step},
      __Dymola_Algorithm="Dassl"))
end EurecaModelicaModel;"""


    with open(os.path.join(".","ExampleFiles","EurecamodelicaModel.mo"), "w") as file:
        file.writelines(modelica_model)


    start = time.time()

    # df_res = bd.simulate(weather_file, output_folder="Results")
    print(f"2C model: \n\t{8760 * 2 - 1} time steps\n\t{(time.time() - start):.2f} s")