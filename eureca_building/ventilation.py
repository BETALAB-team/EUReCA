"""
This module includes functions to model natural ventilation and infiltration
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging

import numpy as np
from scipy.optimize import fsolve

from eureca_building.schedule_properties import ventilation_prop
from eureca_building.schedule import Schedule
from eureca_building.fluids_properties import air_properties, gravitational_acceleration
from eureca_building.weather import WeatherFile
from eureca_building.surface import Surface

from eureca_building.exceptions import (
    InvalidScheduleType
)


def calc_neutral_plane_nat_vent(x, *data):
    a,b,c,h_t,h_b = data
    y = 0.
    for i in range(len(a)):
        y += (-2*a[i]*(np.abs(b*(x - h_t[i]) + c[i]))**(3/2)) / (3*b)
        y += (2*a[i]*(np.abs(b*(x - h_b[i]) + c[i]))**(3/2)) / (3*b)
    return y

class Ventilation:
    """
    Ventilation
    """

    def __init__(
            self,
            name: str,
            unit: str,
            nominal_value: float,
            schedule: Schedule,
            tag: str = None,
    ):
        """
        VentilationObject

        Args:
            name: str
                name
            unit: str
                value of the unit: ["Vol/h", "kg/s", "kg/(m2 s)", "m3/s", "m3/(m2 s)"]
            nominal_value: float
                the value to be multiplied by the schedule
            schedule: Schedule
                Schedule object
            tag: str
                a tag to define the type of internal load
        """
        self.name = name
        self.unit = unit
        self.nominal_value = nominal_value
        self.schedule = schedule
        self.tag = tag

    @property
    def nominal_value(self):
        return self._nominal_value

    @nominal_value.setter
    def nominal_value(self, value):
        try:
            value = float(value)
        except ValueError:
            raise ValueError(f"Ventilation object {self.name}, nominal value not float: {value}")
        self._nominal_value = value

    @property
    def schedule(self):
        return self._schedule

    @schedule.setter
    def schedule(self, value):
        if not isinstance(value, Schedule):
            raise ValueError(f"Ventilation object {self.name}, schedule type not Schedule: {type(value)}")
        if value.schedule_type not in ["dimensionless", ]:
            raise InvalidScheduleType(
                f"Ventilation object  {self.name}, schedule type must be 'Dimensionless': {value.schedule_type}"
            )
        if np.any(np.less(value.schedule, 0.)):
            raise ValueError(
                f"Ventilation object  {self.name}, schedule type has some negative values"
            )

        self._schedule = value

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        if not isinstance(value, str):
            raise ValueError(f"Ventilation object {self.name}, unit is not a string: {type(value)}")
        if value not in ventilation_prop["infiltration"]["unit"]:
            raise ValueError(
                f"Ventilation object  {self.name}, unit must be chosen from: {ventilation_prop['infiltration']['unit']}"
            )
        self._unit = value

    def _get_absolute_value_nominal(self, area=None, volume=None):
        """
        Returns ventilation nominal value in kg/s
        Args:
            area: None
                [m2]: must be provided if the load is area specific
            volume: None
                [m3]: must be provided if the load is volume specific

        Returns:
            float
        """
        air_density = air_properties['density']
        try:
            self.nominal_value_absolute = {
                "Vol/h": self.nominal_value * volume * air_density / 3600,
                "kg/s": self.nominal_value,
                "kg/(m2 s)": self.nominal_value * area,
                "m3/s": self.nominal_value * air_density,
                "m3/(m2 s)": self.nominal_value * area * air_density,
            }[self.unit]
        except TypeError:
            raise AttributeError(
                f"Ventilation object  {self.name}, to calculate ventilation mass flow rate with specific unit you have to provide a volume or an area"
            )
        except KeyError:
            raise ValueError(
                f"Ventilation object  {self.name}, unit must be chosen from: {ventilation_prop['infiltration']['unit']}"
            )

    def get_air_flow_rate(self, area=None, volume=None) -> np.array:
        """
        Calc the air mass flow rate in kg/s
        Args:
            area: float
                area [m2]
            volume: float
                volume [m3]

        Returns:
            np.array

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area, volume=volume)
        return self.nominal_value_absolute * self.schedule.schedule

    def get_vapour_flow_rate(self, weather, area=None, volume=None) -> np.array:
        """
        Calc the vapour mass flow rate in kg/s
        Args:
            weather: Weather
                Weather class object
            area: float
                area [m2]
            volume: float
                volume [m3]

        Returns:
            np.array

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area, volume=volume)
        return self.nominal_value_absolute * self.schedule.schedule * weather.hourly_data['out_air_specific_humidity']

    def get_flow_rate(self, weather, *args, **kwargs) -> list:
        """
        Return the air and vapour flow rate from natural ventilation.
        weather object must be passed
        Args:
            weather: Weather
                weather object used for outdoor specific humidity

        Returns:
            [np.array, np.array, np.array]:
                the schedules: air flow rate [kg/s], vapour [kg/s]

        """
        if "area" not in kwargs.keys():
            area = None
        else:
            area = kwargs['area']
        if "volume" not in kwargs.keys():
            volume = None
        else:
            volume = kwargs['volume']
        vapuor = self.get_vapour_flow_rate(weather, area=area, volume=volume)
        air = self.get_air_flow_rate(area=area, volume=volume)
        return air, vapuor

class Infiltration(Ventilation):
    pass


class NaturalVentilation(Ventilation):
    
    def __init__(
            self,
            name: str,
            unit: str,
            nominal_value: float,
            schedule: Schedule,
            tag: str = None,
            surfaces_with_opening: list = None,
            weather: WeatherFile = None,
    ):
        super().__init__(name,unit,nominal_value,schedule,tag)
        self._get_windows_opening()
        if (weather != None) and (surfaces_with_opening != None):
            self.define_pressure_coef(weather, surfaces_with_opening)
                
    @property
    def schedule(self):
        return self._schedule

    @schedule.setter
    def schedule(self, value):
        if not isinstance(value, Schedule):
            raise ValueError(f"Natural ventilation object {self.name}, schedule type not Schedule: {type(value)}")
        if value.schedule_type not in ["dimensionless",]:
            raise InvalidScheduleType(
                f"Natural ventilation object  {self.name}, schedule type must be 'Dimensionless': {value.schedule_type}"
            )
        if np.any(np.less(value.schedule, 0.)):
            raise ValueError(
                f"Natural ventilation object  {self.name}, schedule type has some negative values"
            )

        if np.any(np.greater(value.schedule, 1.)):
            raise ValueError(
                f"Natural ventilation object  {self.name}, opening schedule type has some values above 1"
            )

        self._schedule = value

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        if not isinstance(value, str):
            raise ValueError(f"Ventilation object {self.name}, unit is not a string: {type(value)}")
        if value not in ventilation_prop["natural"]["unit"]:
            raise ValueError(
                f"Ventilation object  {self.name}, unit must be chosen from: {ventilation_prop['natural']['unit']}"
            )
        self._unit = value

    def _get_absolute_value_nominal(self):
        """
        Returns natural ventilation opening nominal value in a value from 0 to 1
        Args:
            area: None
                [m2]: must be provided if the load is area specific
            volume: None
                [m3]: must be provided if the load is volume specific

        Returns:
            float
        """
        air_density = air_properties['density']
        try:
            self.nominal_value_absolute = {
                "%": self.nominal_value / 100,
                "-": self.nominal_value,
            }[self.unit]
        except TypeError:
            raise AttributeError(
                f"Natural Ventilation object  {self.name}, to calculate natural ventilation the unit must be - or %"
            )
        except KeyError:
            raise ValueError(
                f"Natural Ventilation object  {self.name}, unit must be chosen from: {ventilation_prop['natural']['unit']}"
            )

    def get_air_flow_rate(self):
        """
        Not implemented for Natural Ventilation
        """
        raise NotImplementedError(f"Class Natural Ventilation: get_air_flow_rate method not implmented")

    def get_vapour_flow_rate(self):
        """
        Not implemented for Natural Ventilation
        """
        raise NotImplementedError(f"Class Natural Ventilation: get_vapour_flow_rate method not implmented")

    def get_flow_rate(self):
        """
        Not implemented for Natural Ventilation
        """
        raise NotImplementedError(f"Class Natural Ventilation: get_flow_rate method not implmented")

    def _get_windows_opening(self) -> np.array:
        """
        Calc the windows opening from the input schedule in 0-1 range
        Returns:
            np.array

        """
        try:
            self.nominal_value_absolute
        except AttributeError:
            self._get_absolute_value_nominal()
        self.windows_opening = self.nominal_value_absolute * self.schedule.schedule
        return self.windows_opening

    def define_pressure_coef(self, weather, surfaces_with_opening):
        if not isinstance(weather, WeatherFile):
            raise ValueError(f"Natural Ventilation object: weather parameter is not of WeatherFile class: {type(weather)}")
        
        # TODO: check that number of layers is the same for each surface
        
        self.surfaces_with_opening = surfaces_with_opening
        
        windw_dir = weather.hourly_data["wind_direction"]
        surfaces_with_opening # Type list
        
        aspect_ratio = 1
        pressure_coeff_fun = lambda alfa: (0.603*np.log(1.248 - 0.703*np.sin(alfa/2) -\
                                            1.175*(np.sin(alfa))**2 + \
                                            0.131*(np.sin(2*np.log(aspect_ratio)*alfa))**3 + \
                                            0.769*np.cos(alfa/2)  +\
                                            0.07*np.log(aspect_ratio)**2*(np.sin(alfa/2)**2) + \
                                            0.717*np.cos(alfa/2)**2))

        # def pressure_coeff_fun(alfa, aspect_ratio): 
        #     return (0.603*np.log(1.248 - 0.703*np.sin(alfa/2) -\
        #                                     1.175*(np.sin(alfa))**2 + \
        #                                     0.131*(np.sin(2*np.log(aspect_ratio)*alfa))**3 + \
        #                                     0.769*np.cos(alfa/2)  +\
        #                                     0.07*np.log(aspect_ratio)**2*(np.sin(alfa/2)**2) + \
        #                                     0.717*np.cos(alfa/2)**2))

        for s in surfaces_with_opening:
            if not isinstance(s,Surface):
                raise TypeError(f"Natural Ventilation object: the list of surfaces with opening must be full of Surfaces objects. Surface type: {type(s)}")
            # Calcolo deì coefficient della supeficie
            azimuth_mod = s._azimuth + 180.      
            angle_of_incidence = windw_dir - azimuth_mod - 180
            angle_of_incidence[angle_of_incidence > 360] = angle_of_incidence[angle_of_incidence > 360] - 360     
            angle_of_incidence[angle_of_incidence > 360] = angle_of_incidence[angle_of_incidence > 360] - 360     
            angle_of_incidence[angle_of_incidence < 0] = angle_of_incidence[angle_of_incidence < 0] + 360
            angle_of_incidence[angle_of_incidence < 0] = angle_of_incidence[angle_of_incidence < 0] + 360
            angle_of_incidence[angle_of_incidence > 180] = 360 - angle_of_incidence[angle_of_incidence > 180]
            pressure_coeff = pressure_coeff_fun(angle_of_incidence/180*np.pi)
            s.wind_pressure_coeff = pressure_coeff
            s.angle_of_incidence = angle_of_incidence
            s._c_coeff = s.wind_pressure_coeff*weather.hourly_data["wind_speed"]**2   # Coeff c = wind pressure coeff * (wind speed)^2 for calculation of natural ventilation flow rate

    def get_timestep_ventilation_mass_flow(self, ts, t_zone, weather):
        # To be completed and commented
        
        a_coeff = [s._a_coeff for s in self.surfaces_with_opening]
        b_coeff = (2*(t_zone - weather.hourly_data["out_air_db_temperature"][ts])*gravitational_acceleration)/(t_zone +273.15)   # Coeff b = 2*(t_zona - t_esterna)*g/t_zona for calculation of natural ventilation flow rate
        c_coeff = [s._c_coeff[ts] for s in self.surfaces_with_opening]
        
        z_n = np.zeros(len(self.surfaces_with_opening[0]._h_bottom_windows))
        for p in range(len(self.surfaces_with_opening[0]._h_bottom_windows)):
            h_top = [s._h_top_windows[p] for s in self.surfaces_with_opening]
            h_bottom = [s._h_bottom_windows[p] for s in self.surfaces_with_opening]
            res = fsolve(calc_neutral_plane_nat_vent, 5., args = (a_coeff, b_coeff, c_coeff, h_top, h_bottom))
            # TODO: fix when T_ext and T_int are similar (only wind effect)
            z_n[p] = res[0]
            
        return z_n



class MechanicalVentilation(Ventilation):
    """
    Tha same as Natural/Ventilation object but different check on units
    """

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        if not isinstance(value, str):
            raise ValueError(f"Ventilation object {self.name}, unit is not a string: {type(value)}")
        if value not in ventilation_prop["mechanical"]["unit"]:
            raise ValueError(
                f"Ventilation object  {self.name}, unit must be chosen from: {ventilation_prop['mechanical']['unit']}"
            )
        self._unit = value

    def _get_absolute_value_nominal(self, area=None, volume=None):
        """
        Returns ventilation nominal value in kg/s
        Args:
            area: None
                [m2]: must be provided if the load is area specific
            volume: None
                [m3]: must be provided if the load is volume specific

        Returns:
            float
        """
        air_density = air_properties['density']
        try:
            self.nominal_value_absolute = {
                "Vol/h": self.nominal_value * volume * air_density / 3600,
                "kg/s": self.nominal_value,
                "kg/(m2 s)": self.nominal_value * area,
                "m3/s": self.nominal_value * air_density,
                "m3/(m2 s)": self.nominal_value * area * air_density,
            }[self.unit]
        except TypeError:
            raise AttributeError(
                f"Ventilation object  {self.name}, to calculate ventilation mass flow rate with specific unit you have to provide a volume or an area"
            )
        except KeyError:
            raise ValueError(
                f"Ventilation object  {self.name}, unit must be chosen from: {ventilation_prop['mechanical']['unit']}"
            )