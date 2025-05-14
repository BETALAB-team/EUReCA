"""This module includes functions to model natural ventilation and infiltration
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
    """TODO : Giacomo per favore riempi questa

    Parameters
    ----------
    x : float

    data : tuple

    """
    a,b,c,h_t,h_b = data
    y = 0.
    for i in range(len(a)):
        y += (-2*a[i]*(np.abs(b*(x - h_t[i]) + c[i]))**(3/2)) / (3*b)
        y += (2*a[i]*(np.abs(b*(x - h_b[i]) + c[i]))**(3/2)) / (3*b)
    return y

class Ventilation:
    """Ventilation class
    """

    def __init__(
            self,
            name: str,
            unit: str,
            nominal_value: float,
            schedule: Schedule,
            tag: str = None,
    ):
        """VentilationObject creation

        Parameters
        ----------
        name : str
            name
        unit : str
            value of the unit: ["Vol/h", "kg/s", "kg/(m2 s)", "m3/s", "m3/(m2 s)"]
        nominal_value : float
            the value to be multiplied by the schedule
        schedule : eureca_building.schedule.Schedule
            Schedule object with a fractional schedule
        tag : str, default None
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
        """Returns ventilation nominal value in kg/s. This method might be overridden by a child class

        Parameters
        ----------
        area : float, default None
            [m2]: must be provided if the load is area specific
        volume : float, default None
            [m3]: must be provided if the load is volume specific
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
        """Returns ventilation air flow rate in kg/s. This method might be overridden by a child class

        Parameters
        ----------
        area : float, default None
            [m2]: must be provided if the load is area specific
        volume : float, default None
            [m3]: must be provided if the load is volume specific

        Returns
        ----------
        numpy.array
            air flow rate in kg/s

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area, volume=volume)
        return self.nominal_value_absolute * self.schedule.schedule

    def get_vapour_flow_rate(self, weather, area=None, volume=None) -> np.array:
        """Calc the vapour mass flow rate in kg/s

        Parameters
        ----------
        weather : eureca_building.weather.WeatherFile
            Weather object
        area : float, default None
            [m2]: must be provided if the load is area specific
        volume : float, default None
            [m3]: must be provided if the load is volume specific

        Returns
        ----------
        numpy.array
            vapour flow rate in kg/s
        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area, volume=volume)
        return self.nominal_value_absolute * self.schedule.schedule * weather.hourly_data['out_air_specific_humidity']

    def get_flow_rate(self, weather, *args, **kwargs) -> list:
        """Return the air and vapour flow rate from natural ventilation.
        weather object must be passed

        Parameters
        ----------
        weather : eureca_building.weather.WeatherFile
            Weather object
        args : list
            additional args
        kwargs : dict
            additional kwargs

        Returns
        ----------
        tuple
            tuple of two numpy.array (air and vapour flow rates)

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
    """This is just an inherited version of the Ventilation class, without any change
    """
    pass


class NaturalVentilation(Ventilation):
    """Inheritaded from the Ventilation class
    """
    def __init__(
            self,
            name: str,
            unit: str,
            nominal_value: float,
            schedule: Schedule,
            tag: str = None,
            surfaces_with_opening: list = None,
            weather: WeatherFile = None,
            vol_flow_limit: float = None
    ):
        """Init method. Call the Ventilation super() init and then store few more input

        Parameters
        ----------
        name : str
            name
        unit : str
            value of the unit: ["Vol/h", "kg/s", "kg/(m2 s)", "m3/s", "m3/(m2 s)"]
        nominal_value : float
            the value to be multiplied by the schedule
        schedule : eureca_building.schedule.Schedule
            Schedule object with a fractional schedule
        tag : str, default None
            a tag to define the type of internal load
        surfaces_with_opening : list
            list of eureca_building.surface.Surface objects (those considered for the natural ventilation purposes
        weather : eureca_building.weather.WeatherFile
            Weather object
        vol_flow_limit : float
            Limit to natural vent in m3/s
        """
        super().__init__(name,unit,nominal_value,schedule,tag)
        self._get_windows_opening()
        if (weather != None) and (surfaces_with_opening != None):
            self.define_pressure_coef(weather, surfaces_with_opening)
        
        if vol_flow_limit is not None:
            try:
                vol_flow_limit = float(vol_flow_limit)
            except TypeError:
                raise TypeError(f"Ventilation object {self.name}, maximum air flow rate not number or numeric string: {vol_flow_limit}")
            except ValueError:
                raise ValueError(f"Ventilation object {self.name}, maximum air flow rate not float: {vol_flow_limit}")
            self.vol_flow_limit = vol_flow_limit
        else:
            self.vol_flow_limit = 1e15
                
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
        """Returns natural ventilation nominal value in kg/s. Overrides the parent class method

        Parameters
        ----------
        args
        kwargs
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
        """Not implemented for Natural Ventilation
        """
        raise NotImplementedError(f"Class Natural Ventilation: get_air_flow_rate method not implmented")

    def get_vapour_flow_rate(self):
        """Not implemented for Natural Ventilation
        """
        raise NotImplementedError(f"Class Natural Ventilation: get_vapour_flow_rate method not implmented")

    def get_flow_rate(self):
        """Not implemented for Natural Ventilation
        """
        raise NotImplementedError(f"Class Natural Ventilation: get_flow_rate method not implmented")

    def _get_windows_opening(self) -> np.array:
        """Calc the windows opening from the input schedule in 0-1 range

        Returns
        ----------
        numpy.array
            Wiondow opening schedule [0-1]
        """
        try:
            self.nominal_value_absolute
        except AttributeError:
            self._get_absolute_value_nominal()
        self.windows_opening = self.nominal_value_absolute * self.schedule.schedule
        return self.windows_opening

    def define_pressure_coef(self, weather, surfaces_with_opening):
        """TODO : Per Giacomo compila la documentazione

        Parameters
        ----------
        weather : eureca_building.weather.WeatherFile
            WeatherFile object
        surfaces_with_opening : list
            list of eureca_building.surface.Surface objects (those considered for the natural ventilation purposes

        """
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
        """TODO : Per Giacomo compila la documentazione

        Parameters
        ----------
        ts : int
            time step of simulation
        t_zone : float
            zone temperature [°C]
        weather : eureca_building.weather.WeatherFile
            WeatherFile object
        vol_flow_limit : float
            Upper limit for the NV flow rate [m3/s]
            
        Returns
        ----------
        float
            NV air flow rate [m3/s]
            
        """
        # To be completed and commented
        
        opening_percentage = self.windows_opening[ts]
        # Avoid singularity in neutral plane calc
        if abs(weather.hourly_data["out_air_db_temperature"][ts] - t_zone) <= 0.2:
            __t_zone__ = weather.hourly_data["out_air_db_temperature"][ts] - 0.2
        else:
            __t_zone__ = t_zone
        a_coeff = [s._a_coeff * opening_percentage for s in self.surfaces_with_opening]
        b_coeff = (2*(__t_zone__ - weather.hourly_data["out_air_db_temperature"][ts])*gravitational_acceleration)/(__t_zone__ +273.15)   # Coeff b = 2*(t_zona - t_esterna)*g/t_zona for calculation of natural ventilation flow rate
        c_coeff = [s._c_coeff[ts] for s in self.surfaces_with_opening]
        
        z_n = np.zeros(len(self.surfaces_with_opening[0]._h_bottom_windows))
        vol_flow = np.zeros([len(self.surfaces_with_opening[0]._h_bottom_windows), len(self.surfaces_with_opening)])
        vol_flow_sopra = np.zeros(len(self.surfaces_with_opening[0]._h_bottom_windows))
        vol_inflow = np.zeros(len(self.surfaces_with_opening[0]._h_bottom_windows))
        vol_outflow = np.zeros(len(self.surfaces_with_opening[0]._h_bottom_windows))
        for floor in range(len(self.surfaces_with_opening[0]._h_bottom_windows)):
            h_top = [s._h_top_windows[floor] for s in self.surfaces_with_opening]
            h_bottom = [s._h_bottom_windows[floor] for s in self.surfaces_with_opening]
            res = fsolve(calc_neutral_plane_nat_vent, 5., args = (a_coeff, b_coeff, c_coeff, h_top, h_bottom))
            z_n[floor] = res[0]
            # else:
            # TODO: fix when T_ext and T_int are similar (only wind effect)
            # z_n[floor] = np.nan
            

            # i = 0
            # for s in self.surfaces_with_opening:
            #     #if s._h_bottom_windows[floor] < z_n[floor]:

            #     vol_flow[floor][i] += (2*s._a_coeff * opening_percentage*(np.abs(b_coeff*(z_n[floor] - s._h_bottom_windows[floor]) + s._c_coeff[ts]))**(3/2)) / (3*b_coeff)
            #     #if s._h_top_windows[floor] < z_n[floor]:
            #     vol_flow[floor][i] -= (2*s._a_coeff * opening_percentage*(np.abs(b_coeff*(z_n[floor] - s._h_top_windows[floor]) + s._c_coeff[ts]))**(3/2)) / (3*b_coeff)
            #     i+=1
            
            i = 0
            for s in self.surfaces_with_opening:
                vol_flow[floor][i] += (2*s._a_coeff * opening_percentage*(np.abs(b_coeff*(z_n[floor] - s._h_bottom_windows[floor]) + s._c_coeff[ts]))**(3/2)) / (3*b_coeff)
                vol_flow[floor][i] -= (2*s._a_coeff * opening_percentage*(np.abs(b_coeff*(z_n[floor] - s._h_top_windows[floor]) + s._c_coeff[ts]))**(3/2)) / (3*b_coeff)
                if s._h_bottom_windows[floor] < z_n[floor] < s._h_top_windows[floor]:
                    vol_inflow[floor] += np.abs((2*s._a_coeff * opening_percentage*(np.abs(b_coeff*(z_n[floor] - s._h_bottom_windows[floor]) + s._c_coeff[ts]))**(3/2)) / (3*b_coeff))
                    vol_outflow[floor] += np.abs(2*s._a_coeff * opening_percentage*(np.abs(b_coeff*(z_n[floor] - s._h_top_windows[floor]) + s._c_coeff[ts]))**(3/2)) / (3*b_coeff)
                elif z_n[floor] < s._h_bottom_windows[floor]:
                    if vol_flow[floor][i] > 0:
                        vol_inflow[floor] += vol_flow[floor][i]
                    else:
                        vol_outflow[floor] += vol_flow[floor][i]
                else:
                    if vol_flow[floor][i] > 0:
                        vol_inflow[floor] += vol_flow[floor][i]
                    else:
                        vol_outflow[floor] += vol_flow[floor][i]
                i += 1

            # vol_flow_sopra[floor] = 0
            # for s in self.surfaces_with_opening:
            #     if s._h_bottom_windows[floor] > z_n[floor]:
            #         vol_flow_sopra[floor] -= (2*s._a_coeff * opening_percentage*(np.abs(b_coeff*(z_n[floor] - s._h_bottom_windows[floor]) + s._c_coeff[ts]))**(3/2)) / (3*b_coeff)
            #     if s._h_top_windows[floor] > z_n[floor]:
            #         vol_flow_sopra[floor] += (2*s._a_coeff * opening_percentage*(np.abs(b_coeff*(z_n[floor] - s._h_top_windows[floor]) + s._c_coeff[ts]))**(3/2)) / (3*b_coeff)
        # vol_flow_tot = vol_flow.sum()
        vol_inflow_tot = vol_inflow.sum()
        if vol_inflow_tot > self.vol_flow_limit:
            vol_inflow_tot = self.vol_flow_limit
        # vol_outflow_tot = vol_outflow.sum()
        # mass_flow_tot = vol_flow_tot * air_properties["density"] # kg/s
        # vapour_flow_tot = mass_flow_tot * weather.hourly_data["out_air_specific_humidity"][ts] # kg_vap/s
        # return mass_flow_tot, z_n, vol_inflow_tot, vol_outflow_tot, vol_flow #, vapour_flow_tot
        return vol_inflow_tot

class MechanicalVentilation(Ventilation):
    """The same as Natural/Ventilation object but different check on units
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

    def _get_absolute_value_nominal(self, area=0, volume=0):
        """Calcs ventilation nominal value in kg/s. This method overrides the parent class method

        Parameters
        ----------
        area : float, default None
            [m2]: must be provided if the load is area specific
        volume : float, default None
            [m3]: must be provided if the load is volume specific
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