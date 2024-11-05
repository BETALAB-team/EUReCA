"""
This module includes functions to model the building and it includes the Building class
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import copy
import logging
import os

import numpy as np
import pandas as pd

from eureca_building.config import CONFIG
from eureca_building._auxiliary_function_for_monthly_calc import get_monthly_value_from_annual_vector
from eureca_building.thermal_zone import ThermalZone
from eureca_building.pv_system import PV_system
from eureca_building.solar_thermal_system import SolarThermal_Collector
from eureca_building.weather import WeatherFile
from eureca_building.systems import hvac_heating_systems_classes, hvac_cooling_systems_classes, System, Refrigerator
from eureca_building.exceptions import SimulationError
# %% Building class
class Building:
    """This class is a wrapper for ThermalZone objects and HVAC objects
    """

    def __init__(self, name: str, thermal_zones_list:list,  model:str = "2C"):
        """Constructor of the building class. Memorizes the attributes by means of properties setter.
        Checks also the validity of some attributes

        Parameters
        ----------
        name : str
            Name of the building
        thermal_zone : list
            list of ThermalZone ibjects objects
        model : str, default 2C
            model to be used: 1C or 2C
        """
        # self.PV_systems=[]
        self.name = name
        self._thermal_zones_list = thermal_zones_list
        self._model = model
        self._total_area=sum(tz._net_floor_area for tz in thermal_zones_list)

    @property
    def _thermal_zones_list(self) -> list:
        return self.__thermal_zones_list

    @_thermal_zones_list.setter
    def _thermal_zones_list(self, value: list):
        try:
            value = list(value)
        except ValueError:
            raise TypeError(f"Building {self.name}, the thermal_zone_list must be a list or a tuple: {type(value)}")
        for tz in value:
            if not isinstance(tz, ThermalZone):
                raise TypeError(f"Building {self.name}, non ThermalZone object in thermal_zones_list. ")
        self.__thermal_zones_list = value

    @property
    def _model(self) -> str:
        return self.__model
    
    @_model.setter
    def _model(self, value: str):
        try:
            value = str(value)
        except ValueError:
            raise TypeError(f"Building {self.name}, the model must be a str: {type(value)}")
        if not value in ["1C","2C"]:
            raise TypeError(f"Building {self.name}, model must be 1C or 2C. Model = {value}")
        self.__model = value

    @property
    def heating_system(self) -> System:
        return self._heating_system

    @heating_system.setter
    def heating_system(self, value: System):
        if not isinstance(value, System):
            raise TypeError(f"Building {self.name}, the heating system must be a System object: {type(value)}")
        self._heating_system = value

    @property
    def cooling_system(self) -> System:
        return self._cooling_system

    @cooling_system.setter
    def cooling_system(self, value: System):
        if not isinstance(value, System):
            raise TypeError(f"Building {self.name}, the cooling system must be a System object: {type(value)}")
        self._cooling_system = value
        
    @property
    def refrigerator_system(self) -> System:
        return self._refrigerator_system

    @refrigerator_system.setter
    def refrigerator_system(self, value: System):
        if not isinstance(value, System):
            raise TypeError(f"Building {self.name}, the refrigerator system must be a System object: {type(value)}")
        self._refrigerator_system = value
    
    def set_hvac_system(self, heating_system, cooling_system):
        f"""Sets using roperties the heating and cooling system type (strings)

        Available heating systems: {hvac_heating_systems_classes.keys()}
        Available cooling systems: {hvac_cooling_systems_classes.keys()}

        Parameters
        ----------
        heating_system : str
            string to define building heating system
        cooling_system : str
            string to define building cooling system

        Raises
        ------
        KeyError
            if the heating/cooling system is not included in the available list (See above)
        TypeError
            if the heating system does not comply with the Systems metaclass, which is necessary for simulations

        """
        try:
            self.heating_system = hvac_heating_systems_classes[heating_system](heating_system_key = heating_system)
        except KeyError:
            raise KeyError(f"Building {self.name}, heating system not allowed: current heating system {heating_system}. Available heating systems:\n{hvac_heating_systems_classes.keys()}")
        if not isinstance(self.heating_system, System):
            raise TypeError((f"Building {self.name}, heating system does not comply with System class. The heating system class must be created using System interface"))
        try:
            self.cooling_system = hvac_cooling_systems_classes[cooling_system](cooling_system_key = cooling_system)
        except KeyError:
            raise KeyError(f"Building {self.name}, cooling system not allowed: current cooling system {cooling_system}. Available cooling systems:\n{hvac_cooling_systems_classes.keys()}")
        if not isinstance(self.cooling_system, System):
            raise TypeError((f"Building {self.name}, cooling system does not comply with System class. The cooling system class must be created using System interface"))

        # This overrides the convective and radiative fraction to zone
        for tz in self._thermal_zones_list:
            tz.heating_sigma = self.heating_system.sigma
            tz.cooling_sigma = self.cooling_system.sigma
    def add_refrigeration(self, LT_ratio=1):
        print(1)
        self.refrigeration_system=Refrigerator()
    
    def set_refrigerator_capacity(self, EER_mean, refrigeration_electric_ratio=0.29, working_hour_ratio=0.6):
        self.refrigeration_system.set_system_capacity(self._total_area, EER_mean)
        
    def set_hvac_system_capacity(self, weather_object):
        f"""Calls the thermal zone heating and cooling capacity for all themrmal zones (must be run after the calculation of zone loads)

        Parameters
        ----------
        weather_object : eureca_building.weather.WeatherFile
            WeatherFile object to use to simulate 
            

        Raises
        ------
        SimulationError
            if thermal_zone design load calculation has not been carried out yet

        """
        heating_capacity, cooling_capacity = 0. ,0.
        dhw_flow_rate = 0.
        try:
            for tz in self._thermal_zones_list:
                cooling_capacity += tz.design_sensible_cooling_system_power
                heating_capacity += tz.design_heating_system_power
                dhw_flow_rate += tz.domestic_hot_water_volume_flow_rate
        except AttributeError:
            raise SimulationError(f"""
Building {self.name}: set_hvac_system_capacity method can run only after ThermalZones design load is calculated. 
Please run thermal zones design_sensible_cooling_load and design_heating_load
""")
        self.heating_system.set_system_capacity(heating_capacity, weather_object)
        self.cooling_system.set_system_capacity(cooling_capacity, weather_object)
        self.heating_system.set_dhw_design_capacity_tank(dhw_flow_rate, weather_object)

    def add_pv_system(self, weather_obj):

        '''
        PV production
        '''
        building_surface_list=[]
        for tz in self._thermal_zones_list:
            for s in tz._surface_list:
                if s.surface_type=="Roof":
                    building_surface_list.append(s)

        self.pv_system = PV_system(name=f"Bd {self.name} PV system",
                               weatherobject=weather_obj,
                               surface_list=building_surface_list)

    def add_solar_thermal(self, weather_obj):

        dhw_flow_rate = 0.
        try:
            for tz in self._thermal_zones_list:
                dhw_flow_rate += tz.domestic_hot_water_volume_flow_rate.sum()*3600*1000/(CONFIG.ts_per_hour*365)
                
        except AttributeError:
            raise SimulationError(f"""
                                  Building {self.name}: set_hvac_system_capacity method can run only after ThermalZones design load is calculated. 
                                  Please run thermal zones design_sensible_cooling_load and design_heating_load
                                  """)
                      
        building_surface_list=[]
        for tz in self._thermal_zones_list:
            for s in tz._surface_list:
                building_surface_list.append(s)

        # try: 
        self.heating_system.solar_thermal_system=SolarThermal_Collector(name=f"Bd {self.name} ST system",
                                   dhw=dhw_flow_rate,
                                   weatherobject=weather_obj,
                                   surface_list=building_surface_list)
            
        # except AttributeError:
        #     logging.warning(
        #         f"Bd {self.name} : Add solar thermal should be called after a heating system is created. The simulation will neglect the solar thermal")
 
        
    def solve_timestep(self, t: int, weather: WeatherFile):
        """Runs the thermal zone and hvac systems simulation for the timestep t

        Parameters
        ----------
        t : int
            timestep
        weather_object : eureca_building.weather.WeatherFile
            WeatherFile object to use to simulate
        """

        heat_load, dhw_load, cool_load, air_t, air_rh = 0., 0., 0., 0., 0.
        for tz in self._thermal_zones_list:
            tz.solve_timestep(t, weather, model = self._model)
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

            # For the moment not latent
            if tz.latent_zone_load > 0.:
                heat_load += tz.latent_zone_load
            else:
                cool_load += tz.latent_zone_load

            if tz.latent_zone_load > 0.:
                heat_load += tz.latent_zone_load
            else:
                cool_load += tz.latent_zone_load

            # DHW
            dhw_load += tz.domestic_hot_water_demand[t]
            Refrigerated_end_uses=["supermarket"]


        air_t /= len(self._thermal_zones_list)
        air_rh /= len(self._thermal_zones_list)
        if hasattr(self,'refrigeration_system'):          
            self.refrigeration_system.solve_system( weather,t,T_des=30)
        self.heating_system.solve_system(heat_load, dhw_load, weather, t, air_t, air_rh)
        self.cooling_system.solve_system(cool_load, weather, t, air_t, air_rh)
    def simulate(self,
                 weather_object: WeatherFile,
                 t_start: int = CONFIG.start_time_step,
                 t_stop: int = CONFIG.final_time_step,
                 preprocessing_ts: int = 100 * CONFIG.ts_per_hour,
                 output_folder: str = None,
                 output_type: str = "csv",
                 ):
        """Simulate a period and i stores the outputs. Calls solve_timestep method

        Parameters
        ----------
        weather_object : eureca_building.weather.WeatherFile
            WeatherFile object to use to simulate (must be appliad after the calculation of zone loads
        t_start : int (Default first timestep of simulation)
            starting timestep
        t_stop : int (last timestep of simulation)
            stop timestep
        preprocessing_ts : int
            number of preprocessing timesteps
        output_folder : str, default None
            if not None prints building results in the selected folder
        output_type : str, default "csv"
            parquet or csv as output file

        Returns
        ----------
        pandas.DataFrame
            building time step results
        """
        for tz in self._thermal_zones_list:
            tz.reset_init_values()
            
            
        
        results = {
            'TZ Ta [°C]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ To [°C]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ Tmr [°C]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ RH [-]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ sensible load [W]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ latent load [W]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ AHU pre heater load [W]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ AHU post heater load [W]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ AHU electric load [W]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ DHW volume flow rate [L/s]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ DHW demand [W]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),

            'DHW tank charging mode [-]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'DHW tank charge [Wh]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'DHW tank charge [-]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'DHW tank charging rate [W]' : np.zeros([CONFIG.number_of_time_steps, 1]),

            # 'Storage Tank Charge [%]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Solar Thermal Production [Wh]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Non-Renewable DHW [Wh]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system gas consumption [Nm3]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system oil consumption [L]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system gasoline consumption [L]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system coal consumption [kg]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system wood consumption [kg]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system pellet consumption [kg]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system LPG consumption [kg]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system DH consumption [Wh]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system electric consumption [Wh]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Cooling system electric consumption [Wh]': np.zeros([CONFIG.number_of_time_steps, 1]),
            # 'PV Production [W]': np.zeros([CONFIG.number_of_time_steps, 1]),
            'AHU electric consumption [Wh]': np.zeros([CONFIG.number_of_time_steps, 1]),
            'Appliances electric consumption [Wh]': np.zeros([CONFIG.number_of_time_steps, 1]),
            'Electric consumption [Wh]':np.zeros([CONFIG.number_of_time_steps, 1]),
            'Refrigerator Heat Absorbed [Wh]':np.zeros([CONFIG.number_of_time_steps, 1])
        }
        
        
        
        # # Associate solar thermal to the building
        # self.add_solar_thermal(weather_object)

        
        electric_consumption = np.array([tz.electric_load for tz in self._thermal_zones_list]).sum(axis=0) / CONFIG.ts_per_hour
        results['Appliances electric consumption [Wh]'][:, 0] = electric_consumption[CONFIG.start_time_step:CONFIG.final_time_step]

        results['TZ DHW volume flow rate [L/s]'] = 1000 * np.array([tz.domestic_hot_water_volume_flow_rate for tz in self._thermal_zones_list]).T[CONFIG.start_time_step:CONFIG.final_time_step]
        results['TZ DHW demand [W]'] = np.array([tz.domestic_hot_water_demand for tz in self._thermal_zones_list]).T[CONFIG.start_time_step:CONFIG.final_time_step]

        for t in range(t_start - preprocessing_ts, t_stop):
            self.solve_timestep(t, weather_object)

                 
            results['TZ Ta [°C]'][t - t_start,:] = [tz.zone_air_temperature for tz in self._thermal_zones_list]
            results['TZ To [°C]'][t - t_start,:] = [tz.zone_operative_temperature for tz in self._thermal_zones_list]
            results['TZ Tmr [°C]'][t - t_start,:] = [tz.zone_mean_radiant_temperature for tz in self._thermal_zones_list]
            results['TZ RH [-]'][t - t_start,:] = [tz.zone_air_rel_humidity for tz in self._thermal_zones_list]

            results['TZ sensible load [W]'][t - t_start, :] = [tz.sensible_zone_load for tz in self._thermal_zones_list]
            results['TZ latent load [W]'][t - t_start, :] = [tz.latent_zone_load for tz in self._thermal_zones_list]

            results['TZ AHU pre heater load [W]'][t - t_start, :] = [tz.air_handling_unit.preh_deu_Dem for tz in self._thermal_zones_list]
            results['TZ AHU post heater load [W]'][t - t_start, :] = [tz.air_handling_unit.posth_Dem for tz in self._thermal_zones_list]
            results['TZ AHU electric load [W]'][t - t_start, :] = [tz.AHU_electric_consumption for tz in
                                                                      self._thermal_zones_list]

            results['DHW tank charging mode [-]'][t - t_start, 0] = self.heating_system.charging_mode
            results['DHW tank charge [-]'][t - t_start, 0] = self.heating_system.dhw_tank_current_charge_perc
            results['DHW tank charge [Wh]'][t - t_start, 0] = self.heating_system.dhw_tank_current_charge
            results['Non-Renewable DHW [Wh]'][t - t_start,0] = self.heating_system.dhw_capacity_to_tank
            try:
                results['Solar Thermal Production [Wh]'][t - t_start,0] = self.heating_system.solar_gain_out
            except AttributeError:
                results['Solar Thermal Production [Wh]'][t - t_start, 0] = 0


            results['Heating system gas consumption [Nm3]'][t - t_start,0] = self.heating_system.gas_consumption
            results['Heating system oil consumption [L]'][t - t_start,0] = self.heating_system.oil_consumption
            results['Heating system gasoline consumption [L]'][t - t_start,0] = self.heating_system.gasoline_consumption
            results['Heating system LPG consumption [kg]'][t - t_start,0] = self.heating_system.lpg_consumption
            results['Heating system coal consumption [kg]'][t - t_start,0] = self.heating_system.coal_consumption
            results['Heating system wood consumption [kg]'][t - t_start,0] = self.heating_system.wood_consumption
            results['Heating system pellet consumption [kg]'][t - t_start,0] = self.heating_system.pellet_consumption
            results['Heating system DH consumption [Wh]'][t - t_start,0] = self.heating_system.DH_consumption
            results['Heating system electric consumption [Wh]'][t - t_start,0] = self.heating_system.electric_consumption
            results['Cooling system electric consumption [Wh]'][t - t_start,0] = self.cooling_system.electric_consumption
            results['AHU electric consumption [Wh]'][t - t_start,0] = results['TZ AHU electric load [W]'][t - t_start, :].sum() / CONFIG.ts_per_hour
            if hasattr(self,'refrigeration_system'):          
                results['Refrigerator Heat Absorbed [Wh]'][t - t_start,0]=self.refrigeration_system.heat_absorbed
        # results[ 'Solar Thermal PRoduction [Wh]'] = np.array(self.heating_system.solar_gain)
        # print((np.max(results['Solar Thermal Production [Wh]'])))
        # Saving results

        tz_labels = [res for res in results.keys() if res.startswith("TZ")]
        bd_labels = [res for res in results.keys() if not res.startswith("TZ")]
        tz_names = [tz.name for tz in self._thermal_zones_list]
        columns_tz = pd.MultiIndex.from_product([tz_labels,tz_names])
        columns_bd = pd.MultiIndex.from_product([bd_labels,[f"Bd {self.name}"]])
        Time_index = pd.date_range(start = CONFIG.start_date,periods = CONFIG.number_of_time_steps, freq = f"{CONFIG.time_step}s")
        tz = pd.DataFrame(0., index = range(CONFIG.number_of_time_steps), columns = columns_tz)
        bd = pd.DataFrame(0., index = range(CONFIG.number_of_time_steps), columns = columns_bd)
        total = pd.concat([bd, tz], axis=1)
        for tz_result_label in tz_labels:
            total[tz_result_label] = results[tz_result_label]
        for bd_result_label in bd_labels:
            total[bd_result_label] = results[bd_result_label]
        total.index=Time_index    
        total['Electric consumption [Wh]'] += total["Heating system electric consumption [Wh]"]\
                + total["Cooling system electric consumption [Wh]"] \
                + total["Appliances electric consumption [Wh]"] \
                + total['AHU electric consumption [Wh]']

        # Associate PV to the building
        if hasattr(self, 'pv_system'):
            pv_production=self.pv_system.pv_production()
            [BatteryState , tobattery, frombattery, togrid, fromgrid, directsolar]=self.pv_system.Battery_charge(electricity=total['Electric consumption [Wh]'].iloc[:, 0].values,pv_prod=pv_production)
        else:
            pv_production = 0.
            [BatteryState, tobattery, frombattery, togrid, fromgrid, directsolar] = [np.nan]*6
            fromgrid = total['Electric consumption [Wh]'].iloc[:, 0].values
            togrid = 0.

        total["PV production [Wh]",f"Bd {self.name}"]=pv_production
        total["Battery State [%]",f"Bd {self.name}"]=BatteryState
        total["Given to Batteries [Wh]",f"Bd {self.name}"]=tobattery
        
        total["Taken from the Batteries [Wh]",f"Bd {self.name}"]=frombattery
        total["Given to Grid [Wh]",f"Bd {self.name}"]=togrid
        total["Taken from the Gird [Wh]",f"Bd {self.name}"]=fromgrid
        total["directly from the PV [Wh]",f"Bd {self.name}"]=directsolar
        total["PV System self consumption",f"Bd {self.name}"]=(frombattery+directsolar)/(fromgrid+frombattery+directsolar)

         
        #total = pd.concat([total, pv_production], axis=1)
        #pv_production=tz.pv_production.interpolate(method="time")
        if output_folder != None:
            if not os.path.isdir(output_folder):
                os.mkdir(output_folder)
            if output_type == 'csv':
                total.to_csv(os.path.join(output_folder, f"Results {self.name}.csv"), float_format='%.2f', index = False, sep =";")
            elif output_type == 'parquet':
                total.to_parquet(os.path.join(output_folder, f"Results {self.name}.parquet.snappy"), engine="pyarrow", compression = "snappy")
            else:
                raise KeyError(f"Building simulation: output file type can be either 'csv' or 'parquet'. Current output type: {output_type}")
        return total

    def simulate_quasi_steady_state(self,
                 weather_object: WeatherFile,
                 output_folder: str = None,
                 output_type: str = "csv",
                 ):
        """Simulate a period and i stores the outputs. Calls solve_timestep method

        Parameters
        ----------
        weather_object : eureca_building.weather.WeatherFile
            WeatherFile object to use to simulate (must be appliad after the calculation of zone loads
        output_folder : str, default None
            if not None prints building results in the selected folder
        output_type : str, default "csv"
            parquet or csv as output file

        Returns
        ----------
        pandas.DataFrame
            building time step results
        """
        for tz in self._thermal_zones_list:
            tz.reset_init_values()

        results = {}

        electric_consumption = np.array([tz.electric_load for tz in self._thermal_zones_list]).sum(
            axis=0) / CONFIG.ts_per_hour # Wh

        results['Appliances electric consumption [Wh]'] = get_monthly_value_from_annual_vector(electric_consumption,
        method='sum')

        DHW_Demand = np.array([tz.domestic_hot_water_demand for tz in self._thermal_zones_list]).T.sum(axis = 1) / CONFIG.ts_per_hour
        DHW_Demand = get_monthly_value_from_annual_vector(DHW_Demand,method='sum')

        heat_demand = np.array([0]*12)
        cool_demand = np.array([0]*12)
        for tz in self._thermal_zones_list:
            tz.solve_quasisteadystate_method(weather_object)
            shd = np.clip(tz.heat_sensible_zone_demand_qss_method, 0, None)*1000 # Wh
            lhd = np.clip(tz.heat_latent_zone_demand_qss_method, 0, None)*1000 # Wh
            sad = np.clip(tz.sensible_AHU_demand_qss_method, 0, None)*1000 # Wh
            lad = np.clip(tz.latent_AHU_demand_qss_method, 0, None)*1000 # Wh
            heat_demand = heat_demand + shd + lhd + sad + lad

            shd = np.clip(tz.cool_sensible_zone_demand_qss_method, None, 0)*1000 # Wh
            lhd = np.clip(tz.cool_latent_zone_demand_qss_method, None, 0)*1000 # Wh
            sad = np.clip(tz.sensible_AHU_demand_qss_method, None, 0)*1000 # Wh
            lad = np.clip(tz.latent_AHU_demand_qss_method, None, 0)*1000 # Wh
            cool_demand = cool_demand + shd + lhd + sad + lad

        try:
            self.heating_system.solve_quasi_steady_state(heat_demand, DHW_Demand)
        except AttributeError:
            raise AttributeError("Heating system not allowed. If solving with quasi steady state, heating system must have a quasi steady state method solution... ")
        try:
            self.cooling_system.solve_quasi_steady_state(cool_demand)
        except AttributeError:
            raise AttributeError("Cooling system not allowed. If solving with quasi steady state, cooling system must have a quasi steady state method solution... ")

        results['TZ heating demand [Wh]'] = heat_demand
        results['TZ cooling demand [Wh]'] = cool_demand
        results['TZ DHW demand [Wh]'] = DHW_Demand

        results['Heating system gas consumption [Nm3]'] = self.heating_system.gas_consumption
        results['Heating system oil consumption [L]'] = self.heating_system.oil_consumption
        results['Heating system coal consumption [kg]'] = self.heating_system.coal_consumption
        results['Heating system wood consumption [kg]'] = self.heating_system.wood_consumption
        results['Heating system pellet consumption [kg]'] = self.heating_system.pellet_consumption
        results['Heating system DH consumption [Wh]'] = self.heating_system.DH_consumption
        results['Heating system electric consumption [Wh]'] = self.heating_system.electric_consumption
        results['Cooling system electric consumption [Wh]'] = self.cooling_system.electric_consumption

        results = pd.DataFrame(results)

        # total = pd.concat([total, pv_production], axis=1)
        # pv_production=tz.pv_production.interpolate(method="time")
        if output_folder != None:
            if not os.path.isdir(output_folder):
                os.mkdir(output_folder)
            if output_type == 'csv':
                results.to_csv(os.path.join(output_folder, f"Results {self.name}.csv"), float_format='%.2f', index=False,
                             sep=";")
            elif output_type == 'parquet':
                results.to_parquet(os.path.join(output_folder, f"Results {self.name}.parquet.snappy"), engine="pyarrow",
                                 compression="snappy")
            else:
                raise KeyError(
                    f"Building simulation: output file type can be either 'csv' or 'parquet'. Current output type: {output_type}")

        return results

    def get_geojson_feature_parser(self):
        """Function to get the json dictionary of building properties to stamp the output geojson

        Returns
        ----------
        dict
            dict with some info of the building
        """

        floors = []

        for s in self._thermal_zones_list[0]._surface_list:
            if s.surface_type == "GroundFloor":
                vtxs = s._vertices
                vtxs = [[vtx[0], vtx[1]] for vtx in vtxs]
                floors.append(vtxs)

        return {
            "type": "Feature",
            "properties": {
                "id": self.name,
                "new_id": self.name,
                "Name": self.name,
                },
            "geometry":{
                "type": "MultiPolygon",
                "coordinates": [floors]
            }
        }
    


