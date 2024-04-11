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
        x0,
        *args,
        flag: str = "cal"
        ):
       
    weather, T_meas, start_time_step, end_time_step = args
    
    C_m = x0[0]
    Loads_schedule = x0[1]
    # Opening_NV_schedule = np.array([0]*15+[1]*5+[0]*4) #x0[25:]
    Opening_NV_schedule = x0[2:]
    
    #########################################################
    # Definition of opaque constructions and windows
    
    # path = os.path.join(
    #     "..",
    #     "example_scripts",
    #     "materials_and_construction_test.xlsx",
    # )
    
    ext_wall_North = Construction.from_U_value("NorthExtWall", 0.19, weight_class = "Medium", construction_type = "ExtWall")  # At the moment, I left "medium" weight class; for this wall the solar radiation is almost completely shaded, since it is part of a tunnel
    ext_wall_East = Construction.from_U_value("EastExtWall", 0.19, weight_class = "Medium", construction_type = "ExtWall")
    ext_wall_West = Construction.from_U_value("WestExtWall", 0.19, weight_class = "Medium", construction_type = "ExtWall")
    internal_wall = Construction.from_U_value("InternalWall", 0.3, weight_class = "Medium", construction_type = "IntWall")  # I took a reasonable U-value for the internal walls between the unit and the other flats or common areas
    internal_partition = Construction.from_U_value("InternalPartition", 1.8, weight_class = "Medium", construction_type = "IntWall")  # I took a reasonable U-value for the internal partitions
    internal_floor = Construction.from_U_value("InternalFloor", 0.4, weight_class = "Medium", construction_type = "IntFloor")  # I took a reasonbale U-value for the internal slab
    internal_ceiling = Construction.from_U_value("InternalCeiling", 0.4, weight_class = "Medium", construction_type = "IntCeiling")  # I took a reasonbale U-value for the internal slab
    
    
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
        vertices=((0, 9.5, 3.), (0, 9.5, 5.7), (12.3, 9.5, 5.7), (12.3, 9.5, 3.)),  # The third coordinate takes into account the flat floor number (each floor of 3.0 m - external dimensions)
        wwr=0,
        surface_type="ExtWall",
        construction=ext_wall_North,
        window=window,
        n_window_layers=1,
        h_window=2.1,
        h_bottom=0
    )
    
    wall_East = Surface(
        "Wall East",
        vertices=((12.3, 3.7, 3.), (12.3, 9.5, 3.), (12.3, 9.5, 5.7), (12.3, 3.7, 5.7)),  # The third coordinate takes into account the flat floor number (each floor of 3.0 m - external dimensions)
        wwr=0.35,
        surface_type="ExtWall",
        construction=ext_wall_East,
        window=window,
        n_window_layers=1,
        h_window=2.1,
        h_bottom=0
    )
    
    wall_West = Surface(
        "Wall West",
        vertices=((0., 0., 3.0), (0., 0., 5.7), (0., 9.5, 5.7), (0., 9.5, 3.0)),  # The third coordinate takes into account the flat floor number (each floor of 3.0 m - external dimensions)
        wwr=0.32,
        surface_type="ExtWall",
        construction=ext_wall_West,
        window=window,
        n_window_layers=1,
        h_window=2.1,
        h_bottom=0
    )
        
    int_wall = SurfaceInternalMass(
        "IntWall",
        area=43,
        surface_type="IntWall",
        construction=internal_wall
    )
    
    int_part = SurfaceInternalMass(
        "IntPart",
        area=172,     # Area is doubled since partitions inside the occupied zone present two surfaces
        surface_type="IntWall",
        construction=internal_partition
    )
    
    int_floor = SurfaceInternalMass(
        "IntFloor",
        area=97,  # Floor area based the average value between internal and external dimensions
        surface_type="IntFloor",
        construction=internal_floor
    )
    
    int_ceil = SurfaceInternalMass(
        "IntCeil",
        area=97,  # Ceiling area based the average value between internal and external dimensions
        surface_type="IntCeiling",
        construction=internal_ceiling
    )
    
    #########################################################
    # Loads
    
    ts_h = CONFIG.ts_per_hour
    delay_ts = 8760*ts_h+1-ts_h
    
    # A schedule
    # people_sched = Schedule(
    #     "PeopleOccupancy1",
    #     "percent",
    #     np.array(([0.1] * 7 * ts_h + [0.6] * 2 * ts_h + [0.4] * 5 * ts_h + [0.6] * 10 * ts_h) * 365)[:delay_ts],
    # )
    
    # app_sched = Schedule(          # This schedule actually considers all the sensible loads (including people)
    #     "AppliancesSchedule",
    #     "dimensionless",
    #     np.tile(Loads_schedule, 365)[:delay_ts],
    # )
    
    app_sched = Schedule.from_constant_value(
        name = "app_sched",
        schedule_type="dimensionless",
        value = Loads_schedule
    )
    
    # Loads
    # people = People(
    #     name='occupancy_tz',
    #     unit='px',
    #     nominal_value=1,
    #     schedule=people_sched,
    #     fraction_latent=0.45,
    #     fraction_radiant=0.3,
    #     fraction_convective=0.7,
    #     metabolic_rate=150,
    # )
    
    electric_devices = ElectricLoad(
        name='ElectricLoads',
        unit='W',
        nominal_value=1000,  # An indicative value of 1 kW could be reasonable
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
    # Infiltration
    
    inf_sched = Schedule.from_constant_value(
        name = "nat_vent_sched",
        schedule_type = "dimensionless",
        value = 0,
    )
    
    inf_obj = Infiltration(
        name='inf',
        unit='Vol/h',
        nominal_value=0,
        schedule=inf_sched,
    )
    
    #########################################################
    # Natural ventilation
    
    natural_vent_sched = Schedule(
        "nat_vent_sched",
        "dimensionless",
        np.tile(Opening_NV_schedule, 365)[:delay_ts],
    )
    
    nv_obj = NaturalVentilation(
        name='nat_vent',
        unit='%',
        nominal_value=50, # Windows are considered as closed during the whole period
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
        nominal_value=0.0365,   # Ventilation flow rate in m3/s (volumetric flow rate assumed equal to Ali's apartment at the moment)
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
    
    # Create thermal zone
    tz1 = ThermalZone(
        name="Zone 1",
        surface_list=[wall_North, wall_East, wall_West, int_wall, int_part, int_floor, int_ceil],
        net_floor_area=87, # Net floor area based on internal dimensions
        volume=87*2.7)
    
    tz1._ISO13790_params()
    # tz1._VDI6007_params()

    tz1.Cm = C_m*tz1.Cm
    # tz1.Cm = 0.65*tz1.Cm
    
    tz1.add_internal_load(electric_devices)
    
    # IHG preprocessing
    tz_loads = tz1.extract_convective_radiative_latent_electric_load()
    tz1.calculate_zone_loads_ISO13790(weather)
    
    
    # 2C model
    # tz1.calculate_zone_loads_VDI6007(weather_file)
    
    
    tz1.add_temperature_setpoint(t_sp)
    tz1.add_humidity_setpoint(h_sp)
    
    # Natural Ventilation preprocessing
    tz1.add_natural_ventilation(nv_obj, weather)
    tz1.add_infiltration(inf_obj)
    tz1.calc_infiltration(weather)
    
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
    
    cooling_1C_peak_load = tz1.design_sensible_cooling_load(weather, model = "1C")
    heating_peak_load = tz1.design_heating_load(-5.)
    
    tz1.add_domestic_hot_water(weather, dhw_1)
    
    bd = Building("BdL09", thermal_zones_list=[tz1], model = "1C")
    bd.set_hvac_system("Traditional Gas Boiler, Centralized, Low Temp Radiator", "A-W chiller, Centralized, Radiant surface")
    bd.set_hvac_system_capacity(weather)
    #start = time.time()
    df_res = bd.simulate(weather, t_start = start_time_step, t_stop = end_time_step, preprocessing_ts=100)
    #nat_vent_mass_flow_rate = tz1.nat_vent_air_flow_rate[start_time_step:end_time_step]
    nat_vent_dict = tz1.nat_vent_info
    phi_int_tot = tz_loads["convective [W]"] + tz_loads["radiative [W]"]
    outdoor_temp_profile = weather.hourly_data["out_air_db_temperature"][start_time_step:end_time_step]
    #print(f"1C model: \n\t{8760 * 2 - 1} time steps\n\t{(time.time() - start):.2f} s")
    T_model = df_res["TZ Ta [°C]"]["Zone 1"].values[:len(T_meas)]
    RMSE = np.sqrt(np.sum((T_meas-T_model)**2)/len(T_meas))
    #RMSE_np = np.sqrt(np.square(np.subtract(T_meas,T_model)).mean())
    #global temp_profiles
    temp_profiles = np.array([T_meas, T_model, outdoor_temp_profile]).T
    # print(f'RMSE: {RMSE} °C')
    
    # plt.plot(tz1.phi_ia[start_time_step:end_time_step+100])
    # plt.plot(tz1.phi_m[start_time_step:end_time_step+100])
    # plt.plot(tz1.phi_st[start_time_step:end_time_step+100])
    
    if flag == "cal":
        return RMSE #, RMSE_np, temp_profiles
    else:
        return RMSE, temp_profiles, nat_vent_dict, phi_int_tot # nat_vent_mass_flow_rate

    
    
    
    # hourly = pd.read_csv(os.path.join('Results','Results Bd 1.csv'), header = [0,1], delimiter = ";")
    # hourly.set_index(pd.date_range(start=CONFIG.start_date,end=CONFIG.final_date,periods=CONFIG.number_of_time_steps), inplace = True)
    # hourly.columns = hourly.columns.droplevel(level = 1)
    # hourly_res = hourly[['TZ Ta [°C]','TZ To [°C]','TZ Tmr [°C]','TZ RH [-]']]
    # hourly_res["T_ext"] = weather_file.hourly_data["out_air_db_temperature"][CONFIG.start_time_step:CONFIG.final_time_step]
    # fig, ax1 = plt.subplots(nrows = 1)
    # hourly_res.plot(ax = ax1)
    
    # T_e = pd.DataFrame({"T_ext": weather_file.hourly_data["out_air_db_temperature"][CONFIG.start_time_step:CONFIG.final_time_step]})
    # T_e.plot(ax = ax1)
    # fig, [ax1, ax2,ax3] = plt.subplots(nrows = 3)
    
    # for ax, i in zip ([ax1,ax2,ax3],range(3)):
    #     pd.DataFrame({
    #         hourly_res.columns[i] : hourly_res[hourly_res.columns[i]].values,
    #     }).plot(ax = ax, kind = 'line')


# #########################################################
# #########################################################
# # Calibration process
# #########################################################
# # Input simulation
# # Simulation day
# date_to_calibrate = "06/23/2023" #m/d/y American format

# # Initial vector and boundaries
# C_0 = np.array([1])
# InternalLoad_0 = np.array([0.3]*24)
# x0 = np.hstack([C_0, InternalLoad_0])

# C_0_lb = np.array([0.5])  #0.5
# InternalLoad_0_lb = np.array([0]*24)
# C_0_ub = np.array([3])  #1.5
# InternalLoad_0_ub = np.array([1]*24)
# x0_lb = np.hstack([C_0_lb, InternalLoad_0_lb])
# x0_ub = np.hstack([C_0_ub, InternalLoad_0_ub])

# #########################################################
# # Epw loading
# epw_path = os.path.join('..', 'Ventilation_CPH', 'Weather_CPH_2023.epw')
# weather_file = WeatherFile(epw_path,
#                            time_steps=CONFIG.ts_per_hour,
#                            azimuth_subdivisions=CONFIG.azimuth_subdivisions,
#                            height_subdivisions=CONFIG.height_subdivisions, )

# #########################################################
# # Measure loading
# measure = pd.read_csv("C:\\Users\\gecky\\OneDrive - Università degli Studi di Padova\\PhD_directory\\AAU material\\DataComfortCooling\\Collected_data\\IC-Meter-QR20F6237B-Indoor-Minutes-01-Jun-2023-29-Oct-2023.csv",
#                       skiprows = 0, header = 1, delimiter = ';', decimal = ',', index_col = 0, parse_dates = True)
# measure.drop(["DATE (EUROPE/COPENHAGEN)", "TIME (EUROPE/COPENHAGEN)"], axis = 1, inplace = True)
# measure_h = measure.resample("1H").mean()
# day_of_the_year = measure_h.loc[date_to_calibrate].index.dayofyear
# T_meas = measure_h.loc[date_to_calibrate]["TEMPERATURE (°C)"].values

# #########################################################
# # Simulation for calibration
# #day_of_the_year = datetime.datetime.strptime(date_to_calibrate + " 00:00", '%m/%d/%Y %H:%M').dayofyear
# day_of_the_year = day_of_the_year[0] - 1
# start_time_step = day_of_the_year*24
# end_time_step = day_of_the_year*24 + 24
# RMSE = simulation(x0,
#                  weather = weather_file,
#                  T_meas = T_meas, 
#                  start_time_step = start_time_step, 
#                  end_time_step = end_time_step)
# start = time.time()
# x_opt = scipy.optimize.least_squares(simulation,
#                                       x0,
#                                       bounds = (x0_lb, x0_ub),
#                                       xtol=1e-2,
#                                       method = 'trf', 
#                                       kwargs = {"weather": weather_file,
#                                                 "T_meas": T_meas,
#                                                 "start_time_step": start_time_step,
#                                                 "end_time_step": end_time_step}
#                                       )

# # Running a final simulation to get results
# x0 = x_opt.x
# RMSE_fin = simulation(x0,
#                   weather = weather_file,
#                   T_meas = T_meas, 
#                   start_time_step = start_time_step, 
#                   end_time_step = end_time_step)


# # Trying another scipy function for optimization
# # bounds = scipy.optimize.Bounds(x0_lb, x0_ub)
# # x_opt = scipy.optimize.minimize(simulation,
# #                                 x0,
# #                                 args=(weather_file, T_meas, start_time_step, end_time_step),
# #                                 bounds = bounds,
# #                                 tol=1e-5,
# #                                 )

# # Trying another scipy function for global optimization
# # args=(weather_file, T_meas, start_time_step, end_time_step)
# # x_opt = scipy.optimize.basinhopping(simulation,
# #                                     x0, 
# #                                     niter=1,
# #                                     minimizer_kwargs={"args": args},
# #                                     )


# stop = time.time()
# print(f'Total calibration time for one day: {(stop-start):.1f} s')

# fig, ax = plt.subplots()
# ax.plot(x_opt.x[1:])

# n_timesteps = end_time_step - start_time_step
# time_interval = np.array([range(n_timesteps)]).T
# fig2, ax2 = plt.subplots()
# ax2.plot(time_interval, T_meas, time_interval, temp_profiles[:,1])



#########################################################
#########################################################
# Calibration process
#########################################################
# Input simulation
# Simulation day
start_date_to_calibrate = dt.datetime(year=2023, month=6, day=20)
# start_date_to_calibrate = "06/20/2023" #m/d/y American format

# Initial vector and boundaries
C_0 = np.array([1])
InternalLoad_0 = np.array([0.42])
Opening_NV_0 = np.array([0]*10 + [0.1]*5 + [0.2]*5 + [0]*4)
x0 = np.hstack([C_0, InternalLoad_0, Opening_NV_0])
#x0 = np.hstack([InternalLoad_0, Opening_NV_0])

C_0_lb = np.array([0.5])  #0.5
InternalLoad_0_lb = np.array([0])
Opening_NV_0_lb = np.array([0]*24)
C_0_ub = np.array([3])  #1.5
InternalLoad_0_ub = np.array([1])
Opening_NV_0_ub = np.array([1]*24)
x0_lb = np.hstack([C_0_lb, InternalLoad_0_lb, Opening_NV_0_lb])
x0_ub = np.hstack([C_0_ub, InternalLoad_0_ub, Opening_NV_0_ub])
#x0_lb = np.hstack([InternalLoad_0_lb, Opening_NV_0_lb])
#x0_ub = np.hstack([InternalLoad_0_ub, Opening_NV_0_ub])

#########################################################
# Epw loading
epw_path = os.path.join('..', 'Ventilation_CPH', 'Weather_CPH_2023.epw')
weather_file = WeatherFile(epw_path,
                           time_steps=CONFIG.ts_per_hour,
                           azimuth_subdivisions=CONFIG.azimuth_subdivisions,
                           height_subdivisions=CONFIG.height_subdivisions, )

#########################################################
# Measure loading
measure = pd.read_csv("C:\\Users\\gecky\\OneDrive - Università degli Studi di Padova\\PhD_directory\\AAU material\\DataComfortCooling\\1_Collected_data\\IC-Meter-QR37D39EA3-Indoor-Minutes-01-Jun-2023-29-Oct-2023.csv",
                      skiprows = 0, header = 1, delimiter = ';', decimal = ',', index_col = 0, parse_dates = True)
measure.drop(["DATE (EUROPE/COPENHAGEN)", "TIME (EUROPE/COPENHAGEN)"], axis = 1, inplace = True)
measure_h = measure.resample("1H").mean().ffill()

#########################################################
# Cycle for the single day optimization over a week
# Initializing vectors of results
C_th = pd.DataFrame()
sched_hg = np.array([])
sched_op = np.array([])
temperatures = np.array([[],[],[]]).T
NV_fr = np.array([])
sim_days = 8 #number of simulated days
# Initializing table of results as dictionary to be saved on csv file
results_cal = {"T_meas [°C]" : np.zeros([sim_days*24*CONFIG.ts_per_hour, 1]),
               "T_sim [°C]" : np.zeros([sim_days*24*CONFIG.ts_per_hour, 1]),
               "T_out [°C]" : np.zeros([sim_days*24*CONFIG.ts_per_hour, 1]),
               "C_th_rel [-]" : np.zeros([sim_days*24*CONFIG.ts_per_hour, 1]),
               "Sched_HG [-]" : np.zeros([sim_days*24*CONFIG.ts_per_hour, 1]),
               "Q_hg [W]" : np.zeros([sim_days*24*CONFIG.ts_per_hour, 1]),
               "Sched_op_w [-]" : np.zeros([sim_days*24*CONFIG.ts_per_hour, 1]),
               "NV_ACH [vol/h]" : np.zeros([sim_days*24*CONFIG.ts_per_hour, 1]),
               "NV_vfr [m3/h]" : np.zeros([sim_days*24*CONFIG.ts_per_hour, 1])
               }
start = time.time()
for day_step in range(sim_days):
    time_delta = dt.timedelta(days=day_step)
    date_to_calibrate = start_date_to_calibrate + time_delta
    date_to_calibrate = date_to_calibrate.strftime("%m/%d/%Y")
    
    day_of_the_year = measure_h.loc[date_to_calibrate].index.dayofyear
    T_meas = measure_h.loc[date_to_calibrate]["TEMPERATURE (°C)"].values
    
    #########################################################
    # Simulation for calibration
    #day_of_the_year = datetime.datetime.strptime(date_to_calibrate + " 00:00", '%m/%d/%Y %H:%M').dayofyear
    day_of_the_year = day_of_the_year[0] - 1
    start_time_step = day_of_the_year*24
    end_time_step = day_of_the_year*24 + 24
    RMSE = simulation(x0,
                      weather_file,
                      T_meas,
                      start_time_step,
                      end_time_step)
    
    # x_opt = scipy.optimize.least_squares(simulation,
    #                                       x0,
    #                                       bounds = (x0_lb, x0_ub),
    #                                       ftol=1e-2,
    #                                       method = 'trf', 
    #                                       kwargs = {"weather": weather_file,
    #                                                 "T_meas": T_meas,
    #                                                 "start_time_step": start_time_step,
    #                                                 "end_time_step": end_time_step}
    #                                       )
    
    # x_opt = scipy.optimize.least_squares(simulation,
    #                                       x0,
    #                                       bounds = (x0_lb, x0_ub),
    #                                       ftol=1e-4,
    #                                       method = 'trf',
    #                                       args = (weather_file,
    #                                               T_meas,
    #                                               start_time_step,
    #                                               end_time_step)
    #                                       )
    
    # # Trying another scipy function for optimization
    bounds = scipy.optimize.Bounds(x0_lb, x0_ub)
    # x_opt = scipy.optimize.minimize(simulation,
    #                                 x0,
    #                                 args=(weather_file, T_meas, start_time_step, end_time_step),
    #                                 bounds = bounds,
    #                                 tol=1e-2,
    #                                 )
    
    # # Trying another scipy function for global optimization
    # # args=(weather_file, T_meas, start_time_step, end_time_step)
    # # x_opt = scipy.optimize.basinhopping(simulation,
    # #                                     x0, 
    # #                                     niter=1,
    # #                                     minimizer_kwargs={"args": args},
    # #                                     )
    
    # Trying another scipy function for global optimization
    x_opt = scipy.optimize.differential_evolution(simulation,
                                                  x0 = x0,
                                                  bounds = bounds,
                                                  args=(weather_file, T_meas, start_time_step, end_time_step),
                                                  tol = 1e-3,
                                                  disp = True
                                                  )
    
    # Trying another scipy function for global optimization
    # x_opt = scipy.optimize.direct(simulation,
    #                               bounds = bounds,
    #                               args=(weather_file, T_meas, start_time_step, end_time_step),
    #                               )
    
    
    C_day = pd.DataFrame([x_opt.x[0]], columns=["Thermal capacity"], index=[date_to_calibrate])
    C_th = pd.concat([C_th, C_day])
    sched_hg = np.append(sched_hg, x_opt.x[1])
    sched_op = np.append(sched_op, x_opt.x[2:])
    # sched_hg = np.append(sched_hg, x_opt.x[0:24])
    # sched_op = np.append(sched_op, x_opt.x[24:])
    
    # Running a final simulation to get results
    x_fin = x_opt.x
    RMSE_fin, temp_profiles, NV_dict, IHG = simulation(x_fin,
                                                       weather_file,
                                                       T_meas,
                                                       start_time_step,
                                                       end_time_step,
                                                       flag = "sim")
    
    temperatures = np.append(temperatures, temp_profiles, axis=0)
    NV_fr = np.append(NV_fr, NV_dict['airflow_rate']['vol/h'][start_time_step:end_time_step], axis=0)
    
    # Updating table of results
    results_cal['T_meas [°C]'][0+day_step*24*CONFIG.ts_per_hour:24+day_step*24*CONFIG.ts_per_hour] = np.array([temp_profiles[:,0]]).T
    results_cal['T_sim [°C]'][0+day_step*24*CONFIG.ts_per_hour:24+day_step*24*CONFIG.ts_per_hour] = np.array([temp_profiles[:,1]]).T
    results_cal['T_out [°C]'][0+day_step*24*CONFIG.ts_per_hour:24+day_step*24*CONFIG.ts_per_hour] = np.array([temp_profiles[:,2]]).T
    results_cal['C_th_rel [-]'][0+day_step*24*CONFIG.ts_per_hour:24+day_step*24*CONFIG.ts_per_hour] = np.array([[x_fin[0]]*24*CONFIG.ts_per_hour]).T
    results_cal['Sched_HG [-]'][0+day_step*24*CONFIG.ts_per_hour:24+day_step*24*CONFIG.ts_per_hour] = np.array([[x_fin[1]]*24*CONFIG.ts_per_hour]).T
    results_cal['Q_hg [W]'][0+day_step*24*CONFIG.ts_per_hour:24+day_step*24*CONFIG.ts_per_hour] = np.array([IHG[start_time_step:end_time_step]]).T
    results_cal['Sched_op_w [-]'][0+day_step*24*CONFIG.ts_per_hour:24+day_step*24*CONFIG.ts_per_hour] = np.array([x_fin[2:]]).T
    results_cal['NV_ACH [vol/h]'][0+day_step*24*CONFIG.ts_per_hour:24+day_step*24*CONFIG.ts_per_hour] = np.array([NV_dict['airflow_rate']['vol/h'][start_time_step:end_time_step]]).T
    results_cal['NV_vfr [m3/h]'][0+day_step*24*CONFIG.ts_per_hour:24+day_step*24*CONFIG.ts_per_hour] = np.array([NV_dict['airflow_rate']['m3/h'][start_time_step:end_time_step]]).T

    print(RMSE_fin)
    
    
stop = time.time()
print(f'Total calibration time for {sim_days} days: {(stop-start):.1f} s')

#########################################################
# Final results and graph printing
# Average thermal capacity over the simulated week
C_th_av = C_th.mean()

# Window opening and heat gain schedules
fig, [ax1, ax2] = plt.subplots(nrows=2)
# ax.plot(x_opt.x[1:])
#ax1.plot(sched_hg)
ax2.plot(sched_op)

# Temperature profiles (outdoor, measured and simulated)
starting_timestep = (int(start_date_to_calibrate.strftime("%j"))-1)*24
n_timesteps = end_time_step - starting_timestep
#n_timesteps = 24
time_interval = np.array([range(n_timesteps)]).T
fig2, ax2 = plt.subplots()
ax2.plot(time_interval, temperatures[:,0], time_interval, temperatures[:,1], time_interval, temperatures[:,2])

# Heat gains schedule and value
#HG = sched_hg * 500
#fig, [ax1, ax2] = plt.subplots(nrows=2)
# ax.plot(x_opt.x[1:])
#ax1.plot(HG)
#ax2.plot(sched_op, 'r')

# Window opening schedule and NV air change rates
fig, [ax1, ax2] = plt.subplots(nrows=2)
# ax.plot(x_opt.x[1:])
ax1.plot(sched_op, 'r')
ax2.plot(NV_fr, 'b')

#########################################################
# Saving table of results on a csv file
output_file_name = "Calibration_results_L09"
result_table = pd.DataFrame(0., index = list(range(24*CONFIG.ts_per_hour))*sim_days, columns = results_cal.keys())
for label in results_cal.keys():
    result_table[label] = results_cal[label]
result_table.to_csv(os.path.join('.', 'Results', f"{output_file_name}.csv"), float_format='%.4f', index = True, sep =";")