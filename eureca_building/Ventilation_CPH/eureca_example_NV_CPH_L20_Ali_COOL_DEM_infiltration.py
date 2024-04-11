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



def simulation(
        *args,
        ):
       
    C_m, TotIntGain, weather, start_time_step, end_time_step = args
        
    #########################################################
    # Definition of opaque constructions and windows
    
    # path = os.path.join(
    #     "..",
    #     "example_scripts",
    #     "materials_and_construction_test.xlsx",
    # )
    
    ext_wall_North = Construction.from_U_value("NorthExtWall", 0.19, weight_class = "Medium", construction_type = "ExtWall")  # At the moment, I left "medium" weight class
    ext_wall_East = Construction.from_U_value("EastExtWall", 0.19, weight_class = "Medium", construction_type = "ExtWall")
    roof_constr = Construction.from_U_value("RoofConstr", 0.09, weight_class = "Medium", construction_type = "Roof")
    internal_wall = Construction.from_U_value("InternalWall", 0.3, weight_class = "Medium", construction_type = "IntWall")  # I took a reasonable U-value for the internal walls between the unit and the other flats or common areas
    internal_partition = Construction.from_U_value("InternalPartition", 1.8, weight_class = "Medium", construction_type = "IntWall")  # I took a reasonable U-value for the internal partitions
    internal_floor = Construction.from_U_value("InternalFloor", 0.4, weight_class = "Medium", construction_type = "IntFloor")  # I took a reasonbale U-value for the internal slab
    
    
    window = SimpleWindow(
                    name="Window",
                    u_value = 0.918,
                    solar_heat_gain_coef=0.53,
                    visible_transmittance=0.9,
                    frame_factor=0.10,     # I left a reasonble value
                    shading_coef_int=0.99, # Actually, I don't know this detail, I consider no shading
                    shading_coef_ext=0.99,
    )
    
    #########################################################
    # Definition of surfaces
    wall_North = Surface(
        "Wall North",
        vertices=((0, 12.5, 12.), (0, 12.5, 14.7), (4.2, 12.5, 14.7), (4.2, 12.5, 12.)),  # The third coordinate takes into account the flat floor number (each floor of 3.0 m - external dimensions)
        wwr=0.48,
        surface_type="ExtWall",
        construction=ext_wall_North,
        window=window,
        n_window_layers=1,
        h_window=2.1,
        h_bottom=0
    )
    
    wall_East = Surface(
        "Wall East",
        vertices=((4.2, 0, 12.), (4.2, 12.5, 12.), (4.2, 12.5, 14.7), (4.2, 0, 14.7)),  # The third coordinate takes into account the flat floor number (each floor of 3.0 m - external dimensions)
        wwr=0.32,
        surface_type="ExtWall",
        construction=ext_wall_East,
        window=window,
        n_window_layers=1,
        h_window=2.1,
        h_bottom=0
    )
    
    roof = Surface(
        "Roof",
        vertices=((4.2, 0, 14.7), (4.2, 12.5, 14.7), (0, 12.5, 14.7), (0, 6.5, 14.7), (-2.9, 6.5, 14.7), (-2.9, 2.7, 14.7), (0.9, 2.7, 14.7), (0.9, 0, 14.7)),
        wwr=0,
        surface_type="Roof",
        construction=roof_constr,
        window=window,
        n_window_layers=1
    )
    
    int_wall = SurfaceInternalMass(
        "IntWall",
        area=60,
        surface_type="IntWall",
        construction=internal_wall
    )
    
    int_part = SurfaceInternalMass(
        "IntPart",
        area=80,     # Area is doubled since partitions inside the occupied zone present two surfaces
        surface_type="IntWall",
        construction=internal_partition
    )
    
    int_floor = SurfaceInternalMass(
        "IntFloor",
        area=61,  # Floor area based the average value between internal and external dimensions
        surface_type="IntFloor",
        construction=internal_floor
    )
    
    
    #########################################################
    # Loads
    
    ts_h = CONFIG.ts_per_hour
    delay_ts = 8760*ts_h+1-ts_h
    
    # A schedule for occupants' presence
    people_sched = Schedule(
        "PeopleOccupancy",
        "dimensionless",
        np.array(([1] * 7 * ts_h + [0] * 10 * ts_h + [1] * 7 * ts_h) * 365)[:delay_ts],
    )
    
    # app_sched = Schedule(          # This schedule actually considers all the sensible loads (including people)
    #     "AppliancesSchedule",
    #     "dimensionless",
    #     np.tile(Loads_schedule, 365)[:delay_ts],
    # )
    
    # app_sched = Schedule.from_constant_value(
    #     name = "app_sched",
    #     schedule_type="dimensionless",
    #     value = Loads_schedule
    # )
    
    # A schedule
    app_sched = Schedule(
        "AppSched",
        "dimensionless",
        np.array(([0.43] * 17 * ts_h + [0.84] * 5 * ts_h + [0.43] * 2 * ts_h) * 365)[:delay_ts],
    )
    
    # Setting Internal Heat Loads
    people = People(
        name='occupancy_tz',
        unit='px',
        nominal_value=1,
        schedule=people_sched,
        fraction_latent=0.45,
        fraction_radiant=0.3,
        fraction_convective=0.7,
        metabolic_rate=120,
    )
    
    electric_devices = ElectricLoad(
        name='ElectricLoads',
        unit='W',
        nominal_value=1000*TotIntGain,  # An indicative value of 1 kW could be reasonable
        schedule=app_sched,
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
    
    # Active cooling if Ta in the thermal zone is above 26°C (constant setpoint)
    cool_t = Schedule.from_constant_value(
        name = "t_cool",
        schedule_type="temperature",
        value = 26.
    )
    
    heat_h = Schedule.from_constant_value(
        name = "h_heat",
        schedule_type = "dimensionless",
        value = 0
    )
    
    # Activation of dehumidification above 50% RH (constant setpoint)
    cool_h = Schedule.from_constant_value(
        name = "h_cool",
        schedule_type = "dimensionless",
        value = 0.50
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
    # Infiltration
    
    # inf_sched = Schedule.from_constant_value(
    #     name = "nat_vent_sched",
    #     schedule_type = "dimensionless",
    #     value = 4,
    # )
    
    # Infiltration schedule to have 4 vol/h of constant infiltration when Ali is at home and 0.1 the rest of the day
    inf_sched = Schedule(
        name = "infiltration_sched",
        schedule_type = "dimensionless",
        schedule = np.array(([1] * 7 * ts_h + [0.04] * 10 * ts_h + [1] * 7 * ts_h) * 365)[:delay_ts],
        )
    
    inf_obj = Infiltration(
        name='inf',
        unit='Vol/h',
        nominal_value=5,
        schedule=inf_sched,
    )
    
    #########################################################
    # Natural ventilation
    
    natural_vent_sched = Schedule(
        "nat_vent_sched",
        "dimensionless",
        np.array(([1] * 24 * ts_h)*365)[:delay_ts],
    )
    
    nv_obj = NaturalVentilation(
        name='nat_vent',
        unit='%',
        nominal_value=0, # Windows are considered as closed during the whole period
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
        nominal_value=0.0365,   # Ventilation flow rate in m3/s
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
        value = 0, # 0 equivalent to free-cooling
    )
    
    #########################################################
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
    
    # Creating thermal zones
    tz1 = ThermalZone(
        name="Zone 1",
        surface_list=[wall_North, wall_East, roof, int_wall, int_part, int_floor],
        net_floor_area=56, # Net floor area based on internal dimensions
        volume=56*2.7,
        )
    
    # Calculation of thermal zone RC network parameters
    tz1._ISO13790_params()
    # tz1._VDI6007_params()

    # tz1.Cm = C_m*tz1.Cm
    # Fixed value for thermal capacity found from calibration process
    tz1.Cm = C_m*tz1.Cm
    
    
    # Adding internal heat loads to thermal zone
    tz1.add_internal_load(people, electric_devices)
    # IHG preprocessing
    tz_loads = tz1.extract_convective_radiative_latent_electric_load()
    tz1.calculate_zone_loads_ISO13790(weather)
    
    
    # 2C model
    # tz1.calculate_zone_loads_VDI6007(weather_file)
    
    
    # Adding temperature and humidity setpoints to thermal zone
    tz1.add_temperature_setpoint(t_sp)
    tz1.add_humidity_setpoint(h_sp)
    
    
    # Adding natural ventilation and infiltration to thermal zone
    # Natural Ventilation preprocessing
    tz1.add_natural_ventilation(nv_obj, weather)
    tz1.add_infiltration(inf_obj)
    tz1.calc_infiltration(weather)
    
    
    # Adding air handling unit (mechanical ventilation) to thermal zone
    ahu = AirHandlingUnit(
    "ahu",
    mech_vent_obj,
    T_supply_sched,
    x_supply_sched,
    ahu_availability_sched,
    True,
    0.,
    0.,
    1.,
    weather,
    tz1,
    )
    
    # tz1.add_air_handling_unit(ahu, weather)
    
    
    # Calculation of design heating and sensible cooling systems' power
    cooling_1C_peak_load = tz1.design_sensible_cooling_load(weather, model = "1C")
    heating_peak_load = tz1.design_heating_load(-5.)
    
    
    # Adding DHW to thermal zone
    tz1.add_domestic_hot_water(weather, dhw_1)
    
    
    # Creating building with thermal zones and HVAC systems and arranging the simulation
    bd = Building("Bd L20", thermal_zones_list=[tz1], model = "1C")
    # bd.set_hvac_system("Traditional Gas Boiler, Centralized, Low Temp Radiator", "A-W chiller, Centralized, Radiant surface")
    bd.set_hvac_system("IdealLoad", "IdealLoad")
    bd.set_hvac_system_capacity(weather)
    path_output_folder = "C:\\Users\\gecky\\Desktop\\EUReCA\\eureca_building\\Ventilation_CPH\\Results"
    df_res = bd.simulate(weather, t_start=start_time_step, t_stop=end_time_step, preprocessing_ts=100, output_folder=path_output_folder, output_type="csv")
    nat_vent_dict = tz1.nat_vent_info
    phi_int_tot = tz_loads["convective [W]"] + tz_loads["radiative [W]"]
    outdoor_temp_profile = weather.hourly_data["out_air_db_temperature"][start_time_step:end_time_step]
    T_calc = df_res["TZ Ta [°C]"]["Zone 1"].values[:(end_time_step - start_time_step)]
    temp_profiles = np.array([T_calc, outdoor_temp_profile]).T
    
      
    return temp_profiles, nat_vent_dict, phi_int_tot

    
#########################################################
#########################################################
# Cooling load calculation (CASE 1 - Infiltration fixed daily schedule with maximum of 5 Vol/h)
#########################################################
# Input simulation
# Simulation days
start_date = dt.datetime(year=2023, month=6, day=1)
start_day = start_date.timetuple().tm_yday - 1
start_time_step = start_day*24*CONFIG.ts_per_hour
# start_date = "06/01/2023" #m/d/y American format
stop_date = dt.datetime(year=2023, month=9, day=30)
stop_day = stop_date.timetuple().tm_yday - 1
end_time_step = stop_day*24*CONFIG.ts_per_hour + 24*CONFIG.ts_per_hour

# Input parameters (thermal capacity multiplier and maximum total internal heat gain)
C_m = 2.0553625
TotIntGain= 0.703775

#########################################################
# Epw loading
epw_path = os.path.join('..', 'Ventilation_CPH', 'Weather_CPH_2023.epw')
weather_file = WeatherFile(epw_path,
                           time_steps=CONFIG.ts_per_hour,
                           azimuth_subdivisions=CONFIG.azimuth_subdivisions,
                           height_subdivisions=CONFIG.height_subdivisions,
                           )

#########################################################
# Initializing vectors of results
temperatures = np.array([[],[]]).T
NV_fr = np.array([])

start = time.time()
                                        
temp_profiles, NV_dict, IHG = simulation(C_m,
                                         TotIntGain,
                                         weather_file,
                                         start_time_step,
                                         end_time_step,
                                         )
    
    
temperatures = np.append(temperatures, temp_profiles, axis=0)
NV_fr = np.append(NV_fr, NV_dict['airflow_rate']['vol/h'][start_time_step:end_time_step], axis=0)
    
stop = time.time()
print(f'Total calculation time from {start_date.strftime("%m/%d/%Y") + " 00:00:00"} to {stop_date.strftime("%m/%d/%Y") + " 23:00:00"}: {(stop-start):.1f} s')

#########################################################
# Final results and graph printing

# Temperature profiles (outdoor and simulated)
starting_timestep = (int(start_date.strftime("%j"))-1)*24
n_timesteps = end_time_step - starting_timestep
time_interval = np.array([range(n_timesteps)]).T
fig1, ax1 = plt.subplots()
ax1.plot(time_interval, temperatures[:,0], time_interval, temperatures[:,1])

# Window opening schedule and NV air change rates
fig2, ax2 = plt.subplots(nrows=1)
ax2.plot(NV_fr, 'b')

fig3, ax3 = plt.subplots(nrows=1)
ax3.plot(IHG, 'r')