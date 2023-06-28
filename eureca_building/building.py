"""
This module includes functions to model the building
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
from eureca_building.thermal_zone import ThermalZone
from eureca_building.weather import WeatherFile
from eureca_building.systems import hvac_heating_systems_classes, hvac_cooling_systems_classes, System
from eureca_building.exceptions import SimulationError
# %% Building class
class Building:
    """
    This class is a wrapper for ThermalZone objects and HVAC objects
    """

    def __init__(self, name: str, thermal_zones_list:list,  model:str = "2C"):
        """
        Args:
            name: str
                Name of the building
            thermal_zone: list
                list of ThermalZone ibjects objects
            model: str (default 2C)
                model to be used: 1C or 2C
        """

        self.name = name
        self._thermal_zones_list = thermal_zones_list
        self._model = model

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

    def set_hvac_system(self, heating_system, cooling_system):
        f"""

        Args:
            heating_system: str
                string to define building heating system
            cooling_system: str
                string to define building cooling system

        Available heating systems: {hvac_heating_systems_classes.keys()}
        Available cooling systems: {hvac_cooling_systems_classes.keys()}

        Returns:
            None

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

    def set_hvac_system_capacity(self, weather_object):
        f"""

        Args:
            weather_object: WeatherFile
                WeatherFile object to use to simulate (must be appliad after the calculation of zone loads
            

        Returns:
            None

        """
        heating_capacity, cooling_capacity = 0. ,0.
        try:
            for tz in self._thermal_zones_list:
                cooling_capacity += tz.design_sensible_cooling_system_power
                heating_capacity += tz.design_heating_system_power
        except AttributeError:
            raise SimulationError(f"""
Building {self.name}: set_hvac_system_capacity method can run only after ThermalZones design load is calculated. 
Please run thermal zones design_sensible_cooling_load and design_heating_load
""")
        self.heating_system.set_system_capacity(heating_capacity, weather_object)
        self.cooling_system.set_system_capacity(cooling_capacity, weather_object)

    def solve_timestep(self, t: int, weather: WeatherFile):
        """
        Args:
            t: int
                timestep
            weather_object: WeatherFile
                WeatherFile object to use to simulate (must be appliad after the calculation of zone loads
        """
        heat_load, cool_load, air_t, air_rh = 0., 0., 0., 0.
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
            # if tz.latent_zone_load > 0.:
            #     heat_load += tz.latent_zone_load
            # else:
            #     cool_load += tz.latent_zone_load

            # if tz.latent_zone_load > 0.:
            #     heat_load += tz.latent_zone_load
            # else:
            #     cool_load += tz.latent_zone_load

            # DHW
            heat_load += tz.domestic_hot_water_demand[t]

        air_t /= len(self._thermal_zones_list)
        air_rh /= len(self._thermal_zones_list)

        self.heating_system.solve_system(heat_load, weather, t, air_t, air_rh)
        self.cooling_system.solve_system(cool_load, weather, t, air_t, air_rh)

    def simulate(self,
                 weather_object: WeatherFile,
                 t_start: int = CONFIG.start_time_step,
                 t_stop: int = CONFIG.final_time_step,
                 preprocessing_ts: int = 100 * CONFIG.ts_per_hour,
                 output_folder: str = None
                 ):
        """
        Args:
            weather_object: WeatherFile
                WeatherFile object to use to simulate (must be appliad after the calculation of zone loads
            t_start: int (Default first timestep of simulation)
                starting timestep
            t_stop: int (last timestep of simulation)
                stop timestep
            preprocessing_ts: int
                number of preprocessing timesteps
            output_folder: str (default = None)
                if not None prints building results in the selected folder
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
            'TZ DHW volume flow rate [L/s]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),
            'TZ DHW demand [W]' : np.zeros([CONFIG.number_of_time_steps, len(self._thermal_zones_list)]),

            'Heating system gas consumption [Nm3]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system oil consumption [L]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system wood consumption [kg]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Heating system electric consumption [Wh]' : np.zeros([CONFIG.number_of_time_steps, 1]),
            'Cooling system electric consumption [Wh]': np.zeros([CONFIG.number_of_time_steps, 1]),
            'Appliances electric consumption [Wh]': np.zeros([CONFIG.number_of_time_steps, 1]),
        }

        electric_consumption = np.array([tz.electric_load for tz in self._thermal_zones_list]).sum(axis=0) / CONFIG.ts_per_hour
        results['Appliances electric consumption [Wh]'][:, 0] = electric_consumption[CONFIG.start_time_step:CONFIG.final_time_step]

        results['TZ DHW volume flow rate [L/s]'] = 1000 * np.array([tz.domestic_hot_water_volume_flow_rate for tz in self._thermal_zones_list]).T[CONFIG.start_time_step:CONFIG.final_time_step]
        results['TZ DHW demand [W]'] = np.array([tz.domestic_hot_water_demand for tz in self._thermal_zones_list]).T[CONFIG.start_time_step:CONFIG.final_time_step]

        for t in range(t_start - preprocessing_ts, t_stop):
            self.solve_timestep(t, weather_object)
            results['TZ Ta [°C]'][t,:] = [tz.zone_air_temperature for tz in self._thermal_zones_list]
            results['TZ To [°C]'][t,:] = [tz.zone_operative_temperature for tz in self._thermal_zones_list]
            results['TZ Tmr [°C]'][t,:] = [tz.zone_mean_radiant_temperature for tz in self._thermal_zones_list]
            results['TZ RH [-]'][t,:] = [tz.zone_air_rel_humidity for tz in self._thermal_zones_list]

            results['TZ sensible load [W]'][t, :] = [tz.sensible_zone_load for tz in self._thermal_zones_list]
            results['TZ latent load [W]'][t, :] = [tz.latent_zone_load for tz in self._thermal_zones_list]

            results['TZ AHU pre heater load [W]'][t, :] = [tz.air_handling_unit.preh_deu_Dem for tz in self._thermal_zones_list]
            results['TZ AHU post heater load [W]'][t, :] = [tz.air_handling_unit.posth_Dem for tz in self._thermal_zones_list]


            results['Heating system gas consumption [Nm3]'][t,0] = self.heating_system.gas_consumption
            results['Heating system oil consumption [L]'][t,0] = self.heating_system.oil_consumption
            results['Heating system wood consumption [kg]'][t,0] = self.heating_system.wood_consumption
            results['Heating system electric consumption [Wh]'][t,0] = self.heating_system.electric_consumption
            results['Cooling system electric consumption [Wh]'][t,0] = self.cooling_system.electric_consumption

        # Saving results

        tz_labels = [res for res in results.keys() if res.startswith("TZ")]
        bd_labels = [res for res in results.keys() if not res.startswith("TZ")]
        tz_names = [tz.name for tz in self._thermal_zones_list]
        columns_tz = pd.MultiIndex.from_product([tz_labels,tz_names])
        columns_bd = pd.MultiIndex.from_product([bd_labels,[f"Bd {self.name}"]])
        tz = pd.DataFrame(0., index = range(CONFIG.number_of_time_steps), columns = columns_tz)
        bd = pd.DataFrame(0., index = range(CONFIG.number_of_time_steps), columns = columns_bd)
        total = pd.concat([bd, tz], axis=1)
        for tz_result_label in tz_labels:
            total[tz_result_label] = results[tz_result_label]
        for bd_result_label in bd_labels:
            total[bd_result_label] = results[bd_result_label]

        if output_folder != None:
            if not os.path.isdir(output_folder):
                os.mkdir(output_folder)
            total.to_csv(os.path.join(output_folder, f"Results {self.name}.csv"), float_format='%.2f', index = False)

        return total

    def get_geojson_feature_parser(self):

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


