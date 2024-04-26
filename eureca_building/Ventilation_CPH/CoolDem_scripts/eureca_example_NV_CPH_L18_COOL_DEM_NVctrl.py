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

from eureca_building.fluids_properties import air_properties, vapour_properties
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
       
    C_m, TotIntGain, weather, start_time_step, end_time_step, preprocessing_ts = args
        
    #########################################################
    # Definition of opaque constructions and windows
    
    # path = os.path.join(
    #     "..",
    #     "example_scripts",
    #     "materials_and_construction_test.xlsx",
    # )
    
    ext_wall_West = Construction.from_U_value("WestExtWall", 0.19, weight_class = "Medium", construction_type = "ExtWall")  # At the moment, I left "medium" weight class
    ext_wall_South = Construction.from_U_value("SouthExtWall", 0.19, weight_class = "Medium", construction_type = "ExtWall")
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
    wall_West = Surface(
        "Wall West",
        vertices=((0., 0., 6.), (0., 0., 8.7), (0., 8.6, 8.7), (0., 8.6, 6.)),  # The third coordinate takes into account the flat floor number (each floor of 3.0 m - external dimensions)
        wwr=0.40,
        surface_type="ExtWall",
        construction=ext_wall_West,
        window=window,
        n_window_layers=1,
        h_window=2.1,
        h_bottom=0
    )
    
    wall_South = Surface(
        "Wall South",
        vertices=((0., 0., 6.), (12.3, 0., 6.), (12.3, 0., 8.7), (0., 0., 8.7)),  # The third coordinate takes into account the flat floor number (each floor of 3.0 m - external dimensions)
        wwr=0.33,
        surface_type="ExtWall",
        construction=ext_wall_South,
        window=window,
        n_window_layers=1,
        h_window=2.1,
        h_bottom=0
    )
            
    int_wall = SurfaceInternalMass(
        "IntWall",
        area=56,
        surface_type="IntWall",
        construction=internal_wall
    )
    
    int_part = SurfaceInternalMass(
        "IntPart",
        area=192,     # Area is doubled since partitions inside the occupied zone present two surfaces
        surface_type="IntWall",
        construction=internal_partition
    )
    
    int_floor = SurfaceInternalMass(
        "IntFloor",
        area=100,  # Floor area based the average value between internal and external dimensions
        surface_type="IntFloor",
        construction=internal_floor
    )
    
    int_ceil = SurfaceInternalMass(
        "IntCeil",
        area=100,  # Ceiling area based the average value between internal and external dimensions
        surface_type="IntCeiling",
        construction=internal_ceiling
    )
    
    
    #########################################################
    # Loads
    
    ts_h = CONFIG.ts_per_hour
    delay_ts = 8760*ts_h+1-ts_h
    
    # A schedule for occupants' presence
    people_sched = Schedule(
        "PeopleOccupancy",
        "dimensionless",
        np.array(([1] * 7 * ts_h + [0.5] * 2 * ts_h + [0.1] * 4 * ts_h + [0.2] * 4 * ts_h + [0.5] * 3 * ts_h + [0.8] * 3 * ts_h + [1] * 1 * ts_h) * 365)[:delay_ts],
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
    
    # A schedule for appliances (plausible values from EN 16798‑1:2019 considering occupancy schedule)
    app_sched = Schedule(
        "AppSched",
        "dimensionless",
        np.array(([0.5] * 7 * ts_h + [0.7] * 2 * ts_h + [0.5] * 8 * ts_h + [0.8] * 6 * ts_h + [0.5] * 1 * ts_h) * 365)[:delay_ts],
    )
    
    # A schedule for lighting (plausible values from EN 16798‑1:2019 considering occupancy schedule)
    light_sched = Schedule(
        "LightSched",
        "dimensionless",
        np.array(([0] * 7 * ts_h + [0.15] * 2 * ts_h + [0] * 8 * ts_h + [0.2] * 6 * ts_h + [0] * 1 * ts_h) * 365)[:delay_ts],
    )
    
    
    # Setting Internal Heat Loads
    people = People(
        name='occupancy_tz',
        unit='px',
        nominal_value=4,
        schedule=people_sched,
        fraction_latent=0.45,
        fraction_radiant=0.3,
        fraction_convective=0.7,
        metabolic_rate=120,
    )
    
    electric_devices = ElectricLoad(
        name='ElectricLoads',
        unit='W/m2',
        nominal_value=3.0,  # 3 W/m2 for apartments from Standard EN 16798-1:2019
        schedule=app_sched,
        fraction_radiant=0.3,
        fraction_convective=0.7,
        fraction_to_zone=1
    )
    
    lights = Lights(
        name="Lights",
        nominal_value=8.0,  # 8 W/m2 for apartments from Standard EN 16798-1:2019
        unit="W/m2",
        schedule=light_sched,
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
        value = 26
    )
    
    # cool_t = Schedule(
    #     name = "t_cool",
    #     schedule_type="temperature",
    #     schedule = np.array(([26] * 7 * ts_h + [50] * 10 * ts_h + [26] * 7 * ts_h) * 365)[:delay_ts],
    # )
    
    heat_h = Schedule.from_constant_value(
        name = "h_heat",
        schedule_type = "dimensionless",
        value = 0
    )
    
    # Activation of dehumidification above 60% RH (constant setpoint)
    cool_h = Schedule.from_constant_value(
        name = "h_cool",
        schedule_type = "dimensionless",
        value = 0.60
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
        schedule = np.array(([0] * 7 * ts_h + [1] * 10 * ts_h + [0] * 7 * ts_h) * 365)[:delay_ts],
        )
    
    inf_obj = Infiltration(
        name='inf',
        unit='Vol/h',
        nominal_value=0.3,
        schedule=inf_sched,
    )
    
    #########################################################
    # Natural ventilation
    
    natural_vent_sched = Schedule(
        "nat_vent_sched",
        "dimensionless",
        np.array(([1] * 7 * ts_h + [0] * 10 * ts_h + [1] * 7 * ts_h) * 365)[:delay_ts],
    )
    
    nv_obj = NaturalVentilation(
        name='nat_vent',
        unit='%',
        nominal_value=20, # A maximum of 20% of the glazed surface is considered openable
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
        nominal_value=0.0580,   # Ventilation flow rate in m3/s
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
        surface_list=[wall_West, wall_South, int_wall, int_part, int_floor, int_ceil],
        net_floor_area=89, # Net floor area based on internal dimensions
        volume=89*2.7,
        NV_ACH_limit=4,
        )
    
    # Calculation of thermal zone RC network parameters
    tz1._ISO13790_params()
    # tz1._VDI6007_params()

    # tz1.Cm = C_m*tz1.Cm
    # Fixed value for thermal capacity found from calibration process
    tz1.Cm = C_m*tz1.Cm
    
    
    # Adding internal heat loads to thermal zone
    tz1.add_internal_load(people, electric_devices, lights)
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
    
    tz1.add_air_handling_unit(ahu, weather)
    
    
    # Calculation of design heating and sensible cooling systems' power
    cooling_1C_peak_load = tz1.design_sensible_cooling_load(weather, model = "1C")
    heating_peak_load = tz1.design_heating_load(-5.)
    # Remove design sensible cooling system power
    tz1.design_sensible_cooling_system_power = -1000000  # [W]
    
    
    # Adding DHW to thermal zone
    tz1.add_domestic_hot_water(weather, dhw_1)
    
    
    # Creating building with thermal zones and HVAC systems and arranging the simulation
    bd = Building("L18", thermal_zones_list=[tz1], model = "1C")
    # bd.set_hvac_system("Traditional Gas Boiler, Centralized, Low Temp Radiator", "A-W chiller, Centralized, Radiant surface")
    bd.set_hvac_system("IdealLoad", "IdealLoad")
    bd.set_hvac_system_capacity(weather)
    path_output_folder = "C:\\Users\\gecky\\Desktop\\EUReCA\\eureca_building\\Ventilation_CPH\\Results\\CoolDemRes"
    
    # Building simulation with fixed input
    # df_res = bd.simulate(weather, t_start=start_time_step, t_stop=end_time_step, preprocessing_ts=100, output_folder=path_output_folder, output_type="csv")
    
    # Building simulation with variable window opening for NV
    for tz in bd._thermal_zones_list:
        tz.reset_init_values(T = 25.)

    results = {
        'TZ Ta [°C]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ To [°C]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ Tmr [°C]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ RH [-]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ sensible load [W]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ latent load [W]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ AHU pre heater load [W]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ AHU post heater load [W]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ DHW volume flow rate [L/s]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ DHW demand [W]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ windows opening [-]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ NV volume flowrate [m3/h]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ NV ACH [vol/h]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ MV volume flowrate [m3/h]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ MV ACH [vol/h]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ inf volume flowrate [m3/h]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),
        'TZ inf ACH [vol/h]' : np.zeros([CONFIG.number_of_time_steps, len(bd._thermal_zones_list)]),


        'Heating system gas consumption [Nm3]' : np.zeros([CONFIG.number_of_time_steps, 1]),
        'Heating system oil consumption [L]' : np.zeros([CONFIG.number_of_time_steps, 1]),
        'Heating system wood consumption [kg]' : np.zeros([CONFIG.number_of_time_steps, 1]),
        'Heating system electric consumption [Wh]' : np.zeros([CONFIG.number_of_time_steps, 1]),
        'Cooling system electric consumption [Wh]': np.zeros([CONFIG.number_of_time_steps, 1]),
        'Appliances electric consumption [Wh]': np.zeros([CONFIG.number_of_time_steps, 1]),
    }

    electric_consumption = np.array([tz.electric_load for tz in bd._thermal_zones_list]).sum(axis=0) / CONFIG.ts_per_hour
    results['Appliances electric consumption [Wh]'][:, 0] = electric_consumption[CONFIG.start_time_step:CONFIG.final_time_step]

    results['TZ DHW volume flow rate [L/s]'] = 1000 * np.array([tz.domestic_hot_water_volume_flow_rate for tz in bd._thermal_zones_list]).T[CONFIG.start_time_step:CONFIG.final_time_step]
    results['TZ DHW demand [W]'] = np.array([tz.domestic_hot_water_demand for tz in bd._thermal_zones_list]).T[CONFIG.start_time_step:CONFIG.final_time_step]

    for t in range(start_time_step - preprocessing_ts, end_time_step):
        # Solving single thermal zone balance (with variable window opening schedule for NV)
        heat_load, cool_load, air_t, air_rh = 0., 0., 0., 0.
        for tz in bd._thermal_zones_list:
            # Changing window opening schedule to provide variable NV
            for IHG in tz.internal_loads_list:
                # Condition 1: windows opened only if people are present
                # if isinstance(IHG, People) and IHG.schedule.schedule[t] == 0:
                if isinstance(IHG, People) and IHG.schedule.schedule[t] < 0.5:
                    tz.natural_ventilation.schedule.schedule[t] = 0
                    tz.natural_ventilation.windows_opening[t] = 0 * tz.natural_ventilation.nominal_value_absolute * tz.natural_ventilation.schedule.schedule[t]
                # elif isinstance(IHG, People) and IHG.schedule.schedule[t] != 0:
                elif isinstance(IHG, People) and IHG.schedule.schedule[t] >= 0.5:
                    # Condition 2: windows closed if outdoor temperature is below 16°C (it is unlikely that tenants open windows if outside is colder)
                    if weather.hourly_data['out_air_db_temperature'][t] < 16.:
                        tz.natural_ventilation.schedule.schedule[t] = 0
                        tz.natural_ventilation.windows_opening[t] = 0 * tz.natural_ventilation.nominal_value_absolute * tz.natural_ventilation.schedule.schedule[t]
                    # Condition 3: if outdoor temperature is above 16°C, different sub-conditions in relation to indoor-outdoor temperatures combinations
                    else:
                        # Sub-condition 3.1: window opening control if indoor temperature is equal to or below 22°C
                        if tz.zone_air_temperature <= 22.:
                            if weather.hourly_data['out_air_db_temperature'][t] <= 22.:
                                tz.natural_ventilation.schedule.schedule[t] = 0
                                tz.natural_ventilation.windows_opening[t] = 0 * tz.natural_ventilation.nominal_value_absolute * tz.natural_ventilation.schedule.schedule[t]
                            elif weather.hourly_data['out_air_db_temperature'][t] > 22.:
                                tz.natural_ventilation.schedule.schedule[t] = 1
                                tz.natural_ventilation.windows_opening[t] = 1 * tz.natural_ventilation.nominal_value_absolute * tz.natural_ventilation.schedule.schedule[t]
                        # Sub-condition 3.2: window opening control if indoor temperature is equal to or above 25°C
                        if tz.zone_air_temperature >= 25.:
                            if weather.hourly_data['out_air_db_temperature'][t] <= 26.:
                                tz.natural_ventilation.schedule.schedule[t] = 1
                                tz.natural_ventilation.windows_opening[t] = 1 * tz.natural_ventilation.nominal_value_absolute * tz.natural_ventilation.schedule.schedule[t]
                            elif weather.hourly_data['out_air_db_temperature'][t] > 26.:
                                tz.natural_ventilation.schedule.schedule[t] = 0
                                tz.natural_ventilation.windows_opening[t] = 0 * tz.natural_ventilation.nominal_value_absolute * tz.natural_ventilation.schedule.schedule[t]
                        # Sub-condition 3.3: window opening control if indoor temperature is between 22 and 25°C, taking into account the window opening state at the previous timestep
                        else:
                            if tz.natural_ventilation.schedule.schedule[t-1] == 0:
                                tz.natural_ventilation.schedule.schedule[t] = 0
                                tz.natural_ventilation.windows_opening[t] = 0 * tz.natural_ventilation.nominal_value_absolute * tz.natural_ventilation.schedule.schedule[t]
                            elif tz.natural_ventilation.schedule.schedule[t-1] == 1 and weather.hourly_data['out_air_db_temperature'][t] <= 26.:
                                tz.natural_ventilation.schedule.schedule[t] = 1
                                tz.natural_ventilation.windows_opening[t] = 1 * tz.natural_ventilation.nominal_value_absolute * tz.natural_ventilation.schedule.schedule[t]
                            else:
                                tz.natural_ventilation.schedule.schedule[t] = 0
                                tz.natural_ventilation.windows_opening[t] = 0 * tz.natural_ventilation.nominal_value_absolute * tz.natural_ventilation.schedule.schedule[t]
                        
                else:
                    pass
            
            # Setting infiltration occurrence complementary to window opening for NV 
            tz.infiltration_air_flow_rate[t] = 0
            for inf in tz.infiltration_list:
                if tz.natural_ventilation.windows_opening[t] == 0:
                    inf.schedule.schedule[t] = 1
                else:
                    inf.schedule.schedule[t] = 0
                tz.infiltration_air_flow_rate[t] += inf.get_flow_rate(weather, area=tz._net_floor_area, volume=tz._volume)[0][t]
        
                
            tz.solve_timestep(t, weather, model = bd._model)
            air_t += tz.zone_air_temperature
            air_rh += tz.zone_air_rel_humidity
            if tz.sensible_zone_load > 0.:
                heat_load += tz.sensible_zone_load
            else:
                cool_load += tz.sensible_zone_load

            if tz.air_handling_unit.preh_deu_Dem > 0.:
                heat_load += tz.air_handling_unit.preh_deu_Dem
            else:
                cool_load += tz.air_handling_unit.preh_deu_Dem
            heat_load += tz.air_handling_unit.posth_Dem
            
            heat_load += tz.domestic_hot_water_demand[t]
        
        air_t /= len(bd._thermal_zones_list)
        air_rh /= len(bd._thermal_zones_list)
        bd.heating_system.solve_system(heat_load, weather, t, air_t, air_rh)
        bd.cooling_system.solve_system(cool_load, weather, t, air_t, air_rh)
        
        # Solving entire building balance (to be used if it consists of one thermal zone)
        # bd.solve_timestep(t, weather)
        results['TZ Ta [°C]'][t - start_time_step,:] = [tz.zone_air_temperature for tz in bd._thermal_zones_list]
        results['TZ To [°C]'][t - start_time_step,:] = [tz.zone_operative_temperature for tz in bd._thermal_zones_list]
        results['TZ Tmr [°C]'][t - start_time_step,:] = [tz.zone_mean_radiant_temperature for tz in bd._thermal_zones_list]
        results['TZ RH [-]'][t - start_time_step,:] = [tz.zone_air_rel_humidity for tz in bd._thermal_zones_list]

        results['TZ sensible load [W]'][t - start_time_step, :] = [tz.sensible_zone_load for tz in bd._thermal_zones_list]
        results['TZ latent load [W]'][t - start_time_step, :] = [tz.latent_zone_load for tz in bd._thermal_zones_list]

        results['TZ AHU pre heater load [W]'][t - start_time_step, :] = [tz.air_handling_unit.preh_deu_Dem for tz in bd._thermal_zones_list]
        results['TZ AHU post heater load [W]'][t - start_time_step, :] = [tz.air_handling_unit.posth_Dem for tz in bd._thermal_zones_list]
        
        results['TZ windows opening [-]'][t - start_time_step, :] = [tz.nat_vent_info['windows_opening']['open_fraction'][t] for tz in bd._thermal_zones_list]
        results['TZ NV volume flowrate [m3/h]'][t - start_time_step, :] = [tz.nat_vent_info['airflow_rate']['m3/h'][t] for tz in bd._thermal_zones_list]
        results['TZ NV ACH [vol/h]'][t - start_time_step, :] = [tz.nat_vent_info['airflow_rate']['vol/h'][t] for tz in bd._thermal_zones_list]
        results['TZ MV volume flowrate [m3/h]'][t - start_time_step, :] = [tz.air_handling_unit.air_flow_rate_kg_S[t]/air_properties["density"]*3600 for tz in bd._thermal_zones_list]
        results['TZ MV ACH [vol/h]'][t - start_time_step, :] = [tz.air_handling_unit.air_flow_rate_kg_S[t]/air_properties["density"]*3600/tz._volume for tz in bd._thermal_zones_list]
        results['TZ inf volume flowrate [m3/h]'][t - start_time_step, :] = [tz.infiltration_air_flow_rate[t]/air_properties["density"]*3600 for tz in bd._thermal_zones_list]
        results['TZ inf ACH [vol/h]'][t - start_time_step, :] = [tz.infiltration_air_flow_rate[t]/air_properties["density"]*3600/tz._volume for tz in bd._thermal_zones_list]

        results['Heating system gas consumption [Nm3]'][t - start_time_step,0] = bd.heating_system.gas_consumption
        results['Heating system oil consumption [L]'][t - start_time_step,0] = bd.heating_system.oil_consumption
        results['Heating system wood consumption [kg]'][t - start_time_step,0] = bd.heating_system.wood_consumption
        results['Heating system electric consumption [Wh]'][t - start_time_step,0] = bd.heating_system.electric_consumption
        results['Cooling system electric consumption [Wh]'][t - start_time_step,0] = bd.cooling_system.electric_consumption

    # Saving results
    tz_labels = [res for res in results.keys() if res.startswith("TZ")]
    bd_labels = [res for res in results.keys() if not res.startswith("TZ")]
    tz_names = [tz.name for tz in bd._thermal_zones_list]
    columns_tz = pd.MultiIndex.from_product([tz_labels,tz_names])
    columns_bd = pd.MultiIndex.from_product([bd_labels,[f"Bd {bd.name}"]])
    tz_res = pd.DataFrame(0., index = range(CONFIG.number_of_time_steps), columns = columns_tz)
    bd_res = pd.DataFrame(0., index = range(CONFIG.number_of_time_steps), columns = columns_bd)
    total = pd.concat([bd_res, tz_res], axis=1)
    for tz_result_label in tz_labels:
        total[tz_result_label] = results[tz_result_label]
    for bd_result_label in bd_labels:
        total[bd_result_label] = results[bd_result_label]
    
    # Creating output file
    total.to_csv(os.path.join(path_output_folder, f"Results {bd.name}.csv"), float_format='%.2f', index = False, sep =';')
    
    
    # Output of the apartment simulation to be plotted (only tz1)
    nat_vent_dict = tz1.nat_vent_info
    infiltrations = tz1.infiltration_air_flow_rate
    phi_int_tot = tz_loads["convective [W]"] + tz_loads["radiative [W]"]
    outdoor_temp_profile = weather.hourly_data["out_air_db_temperature"][start_time_step:end_time_step]
    # T_calc = df_res["TZ Ta [°C]"]["Zone 1"].values[:(end_time_step - start_time_step)]
    T_calc = total["TZ Ta [°C]"]["Zone 1"].values[:(end_time_step - start_time_step)]
    temp_profiles = np.array([T_calc, outdoor_temp_profile]).T
    phi_ia = tz1.phi_ia[start_time_step:end_time_step]
    phi_m = tz1.phi_m[start_time_step:end_time_step]
    phi_st = tz1.phi_st[start_time_step:end_time_step]
      
    return temp_profiles, nat_vent_dict, infiltrations, phi_int_tot, (phi_ia, phi_m, phi_st)

    
#########################################################
#########################################################
# Cooling load calculation (CASE 1 - Infiltration fixed daily schedule with maximum of 5 Vol/h)
#########################################################
# Input simulation
# Simulation period
start_date = dt.datetime(year=2023, month=6, day=1)
start_day = start_date.timetuple().tm_yday - 1
start_time_step = start_day*24*CONFIG.ts_per_hour
# start_date = "06/01/2023" #m/d/y American format
stop_date = dt.datetime(year=2023, month=9, day=30)
stop_day = stop_date.timetuple().tm_yday - 1
end_time_step = stop_day*24*CONFIG.ts_per_hour + 24*CONFIG.ts_per_hour
# Preprocessing timesteps
preprocessing_ts = 100

# Input parameters (thermal capacity and maximum total internal heat gain multipliers)
C_m = 1.5
# C_m = 1.5
# TotIntGain= 0.5
TotIntGain = 0.5
# TotIntGain = 0.3

#########################################################
# Epw loading
epw_path = os.path.join('..', 'Weather_CPH_2023.epw')
weather_file = WeatherFile(epw_path,
                           year=2023,
                           time_steps=CONFIG.ts_per_hour,
                           azimuth_subdivisions=CONFIG.azimuth_subdivisions,
                           height_subdivisions=CONFIG.height_subdivisions,
                           )

#########################################################
# Initializing vectors of results and running apartment simulation
temperatures = np.array([[],[]]).T
NV_fr = np.array([])
inf_mfr = np.array([])

start = time.time()
                                        
temp_profiles, NV_dict, infiltrations, IHG, phi_nodes = simulation(C_m,
                                                                   TotIntGain,
                                                                   weather_file,
                                                                   start_time_step,
                                                                   end_time_step,
                                                                   preprocessing_ts,
                                                                   )
    
    
temperatures = np.append(temperatures, temp_profiles, axis=0)
NV_fr = np.append(NV_fr, NV_dict['airflow_rate']['vol/h'][start_time_step:end_time_step], axis=0)
inf_mfr = np.append(inf_mfr, infiltrations[start_time_step:end_time_step], axis=0)

stop = time.time()
print(f'Total calculation time from {start_date.strftime("%m/%d/%Y") + " 00:00:00"} to {stop_date.strftime("%m/%d/%Y") + " 23:00:00"}: {(stop-start):.1f} s')

#########################################################
# Final results and graph printing

# Temperature profiles (outdoor and simulated)
starting_timestep = (int(start_date.strftime("%j"))-1)*24*CONFIG.ts_per_hour
n_timesteps = end_time_step - starting_timestep
time_interval = np.array([range(n_timesteps)]).T
fig1, ax1 = plt.subplots()
ax1.plot(time_interval, temperatures[:,0], time_interval, temperatures[:,1])

# NV air change rates and infiltration mass flow rate
fig2, ax2 = plt.subplots(nrows=1)
ax2.plot(time_interval, NV_fr, 'b', time_interval, inf_mfr, 'r')

# Internal heat gains
fig3, ax3 = plt.subplots(nrows=1)
ax3.plot(IHG, 'r')

# Zone loads on the RC network nodes
fig4, ax4 = plt.subplots(nrows=1)
ax4.plot(time_interval, phi_nodes[0], 'b', time_interval, phi_nodes[1], 'r', time_interval, phi_nodes[2], 'g')

# Angle of incidence and irradiances on vertical surfaces and horizontal plane
fig5, ax5 = plt.subplots(nrows=1)
ax5.plot(weather_file.hourly_data_irradiances[0][0]['global'][start_time_step:end_time_step])
ax5.plot(weather_file.hourly_data_irradiances[-90][90]['global'][start_time_step:end_time_step])
ax5.plot(weather_file.hourly_data_irradiances[0][90]['global'][start_time_step:end_time_step])
ax5.plot(weather_file.hourly_data_irradiances[90][90]['global'][start_time_step:end_time_step])
ax5.plot(weather_file.hourly_data_irradiances[-180][90]['global'][start_time_step:end_time_step])

fig6, ax6 = plt.subplots(nrows=1)
ax6.plot(weather_file.hourly_data_irradiances[0][0]['direct'][start_time_step:end_time_step])
ax6.plot(weather_file.hourly_data_irradiances[-90][90]['direct'][start_time_step:end_time_step])
ax6.plot(weather_file.hourly_data_irradiances[0][90]['direct'][start_time_step:end_time_step])
ax6.plot(weather_file.hourly_data_irradiances[90][90]['direct'][start_time_step:end_time_step])
ax6.plot(weather_file.hourly_data_irradiances[-180][90]['direct'][start_time_step:end_time_step])

fig7, ax7 = plt.subplots(nrows=1)
ax7.plot(weather_file.hourly_data_irradiances[0][0]['AOI'][start_time_step:end_time_step])
ax7.plot(weather_file.hourly_data_irradiances[-90][90]['AOI'][start_time_step:end_time_step])
ax7.plot(weather_file.hourly_data_irradiances[0][90]['AOI'][start_time_step:end_time_step])
ax7.plot(weather_file.hourly_data_irradiances[90][90]['AOI'][start_time_step:end_time_step])
ax7.plot(weather_file.hourly_data_irradiances[-180][90]['AOI'][start_time_step:end_time_step])