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

import datetime
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
        weather=None,
        T_meas=None,
        start_time_step=None,
        end_time_step=None
        ):
        
    C_m = x0[0]
    Loads_schedule = x0[1:25]
    
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
    int_wall = Construction.from_U_value("InternalWall", 0.3, weight_class = "Medium", construction_type = "IntWall")  # I took a reasonable U-value for the internal walls between the unit and the other flats or common areas
    int_partition = Construction.from_U_value("InternalPartition", 1.8, weight_class = "Medium", construction_type = "IntWall")  # I took a reasonable U-value for the internal partitions
    int_floor = Construction.from_U_value("InternalFloor", 0.4, weight_class = "Medium", construction_type = "IntFloor")  # I took a reasonbale U-value for the internal slab
    
    
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
        vertices=((0, 12.5, 13.2), (0, 12.5, 15.9), (4.2, 12.5, 15.9), (4.2, 12.5, 13.2)),  # The third coordinate takes into account the flat floor number (each floor of 3.3 m - external dimensions)
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
        vertices=((4.2, 0, 13.2), (4.2, 12.5, 13.2), (4.2, 12.5, 15.9), (4.2, 0, 15.9)),  # The third coordinate takes into account the flat floor number (each floor of 3.3 m - external dimensions)
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
        vertices=((4.2, 0, 15.9), (4.2, 12.5, 15.9), (0, 12.5, 15.9), (0, 6.5, 15.9), (-2.9, 6.5, 15.9), (-2.9, 2.7, 15.9), (0.9, 2.7, 15.9), (0.9, 0, 15.9)),
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
        construction=int_wall
    )
    
    int_part = SurfaceInternalMass(
        "IntPart",
        area=40,
        surface_type="IntWall",
        construction=int_partition
    )
    
    int_floor = SurfaceInternalMass(
        "IntFloor",
        area=61,  # Floor area based the average value between internal and external dimensions
        surface_type="IntFloor",
        construction=int_floor
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
    
    app_sched = Schedule(          # This schedule actually considers all the sensible loads (including people)
        "AppliancesSchedule",
        "dimensionless",
        np.tile(Loads_schedule, 365)[:delay_ts],
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
        np.array(([.3] * 8 * ts_h + [.5] * 2 * ts_h + [.3] * 4 * ts_h + [.5] * 10 * ts_h) * 365)[:delay_ts],
    )
    
    nv_obj = NaturalVentilation(
        name='nat_vent',
        unit='%',
        nominal_value=0, # Windows are considered as close during the whole period
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
    
    # Create thermal zone
    tz1 = ThermalZone(
        name="Zone 1",
        surface_list=[wall_North, wall_East, roof, int_wall, int_part, int_floor],
        net_floor_area=56, # Net floor area based on internal dimensions
        volume=56*2.7)
    
    tz1._ISO13790_params()
    # tz1._VDI6007_params()
    
    tz1.Cm = C_m*tz1.Cm
    
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
    
    tz1.add_air_handling_unit(ahu, weather)
    
    cooling_1C_peak_load = tz1.design_sensible_cooling_load(weather, model = "1C")
    heating_peak_load = tz1.design_heating_load(-5.)
    
    tz1.add_domestic_hot_water(weather, dhw_1)
    
    bd = Building("Bd1", thermal_zones_list=[tz1], model = "1C")
    bd.set_hvac_system("Traditional Gas Boiler, Centralized, Low Temp Radiator", "A-W chiller, Centralized, Radiant surface")
    bd.set_hvac_system_capacity(weather)
    start = time.time()
    df_res = bd.simulate(weather, t_start = start_time_step, t_stop = end_time_step, preprocessing_ts=50)
    #print(f"1C model: \n\t{8760 * 2 - 1} time steps\n\t{(time.time() - start):.2f} s")
    T_model = df_res["TZ Ta [°C]"]["Zone 1"].values[:len(T_meas)]
    RMSE = np.sqrt(np.sum((T_meas-T_model)**2)/len(T_meas))
    RMSE_np = np.sqrt(np.square(np.subtract(T_meas,T_model)).mean())
    temp_profiles = np.array([T_meas, T_model]).T
    # print(f'RMSE: {RMSE} °C')
    
    return RMSE #, RMSE_np, temp_profiles

    
    
    
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


#########################################################
#########################################################
# Calibration process
#########################################################
# Input simulation
# Simulation day
date_to_calibrate = "06/22/2023" #m/d/y American format

# Initial vector and boundaries
C_0 = np.array([1])
InternalLoad_0 = np.array([0.3]*24)
x0 = np.hstack([C_0, InternalLoad_0])

C_0_lb = np.array([0.5])  #0.5
InternalLoad_0_lb = np.array([0]*24)
C_0_ub = np.array([1.5])  #1.5
InternalLoad_0_ub = np.array([1]*24)
x0_lb = np.hstack([C_0_lb, InternalLoad_0_lb])
x0_ub = np.hstack([C_0_ub, InternalLoad_0_ub])

#########################################################
# Epw loading
epw_path = os.path.join('..', 'Ventilation_CPH', 'Weather_CPH_2023.epw')
weather_file = WeatherFile(epw_path,
                           time_steps=CONFIG.ts_per_hour,
                           azimuth_subdivisions=CONFIG.azimuth_subdivisions,
                           height_subdivisions=CONFIG.height_subdivisions, )

#########################################################
# Measure loading
measure = pd.read_csv("C:\\Users\\gecky\\OneDrive - Università degli Studi di Padova\\PhD_directory\\AAU material\\DataComfortCooling\\Collected_data\\IC-Meter-QR20F6237B-Indoor-Minutes-01-Jun-2023-29-Oct-2023.csv",
                      skiprows = 0, header = 1, delimiter = ';', decimal = ',', index_col = 0, parse_dates = True)
measure.drop(["DATE (EUROPE/COPENHAGEN)", "TIME (EUROPE/COPENHAGEN)"], axis = 1, inplace = True)
measure_h = measure.resample("1H").mean()
day_of_the_year = measure_h.loc[date_to_calibrate].index.dayofyear
T_meas = measure_h.loc[date_to_calibrate]["TEMPERATURE (°C)"].values

#########################################################
# Simulation for calibration
#day_of_the_year = datetime.datetime.strptime(date_to_calibrate + " 00:00", '%m/%d/%Y %H:%M').dayofyear
day_of_the_year = day_of_the_year[0] - 1
start_time_step = day_of_the_year*24
end_time_step = day_of_the_year*24 + 24
RMSE= simulation(x0,
                 weather = weather_file,
                 T_meas = T_meas, 
                 start_time_step = start_time_step, 
                 end_time_step = end_time_step)
start = time.time()
x_opt = scipy.optimize.least_squares(simulation,
                                      x0,
                                      bounds = (x0_lb, x0_ub),
                                      xtol=1e-2,
                                      method = 'trf', 
                                      kwargs = {"weather": weather_file,
                                                "T_meas": T_meas,
                                                "start_time_step": start_time_step,
                                                "end_time_step": end_time_step})

# Trying another scipy function for optimization
# x_opt = 

stop = time.time()
print(f'Total calibration time for one day: {(stop-start):.1f} s')

fig, ax = plt.subplots()
ax.plot(x_opt.x[1:])

fig2, ax2 = plt.subplots()
ax2.plot(T_meas)