"""
This module includes functions and classes to model internal heat gains
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging

import numpy as np

from eureca_building.schedule_properties import internal_loads_prop
from eureca_building.schedule import Schedule
from eureca_building.fluids_properties import vapour_properties

from eureca_building.exceptions import (
    ConvectiveRadiantFractionError,
    InvalidHeatGainUnit,
    InvalidScheduleType,
    AreaNotProvided,
    PeopleNotProvided,
)


# TODO: implement interface for methods to be implemented
# https://realpython.com/python-interface/

class InternalLoad:
    """Internal Gain Class: parent class to set some common things
    """

    def __init__(
            self,
            name: str,
            nominal_value: float,
            schedule: Schedule,
            tag: str = None,
    ):
        """Parent class for some inherited InternalLoads.
        It load the input and checks them throughout properties setter

        Parameters
        ----------
        name : str
            name
        nominal_value : float
            the value to be multiplied by the schedule
        schedule : eureca_building.Schedule
            Schedule object
        tag : str
            a tag to define the type of internal load

        """
        self.name = name
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
            raise ValueError(f"Internal Heat Gain {self.name}, nominal value not float: {value}")
        self._nominal_value = value

    @property
    def schedule(self):
        return self._schedule

    @schedule.setter
    def schedule(self, value):
        if not isinstance(value, Schedule):
            raise ValueError(f"Internal Heat Gain {self.name}, schedule type not Schedule: {type(value)}")
        if value.schedule_type not in ["dimensionless", "percent", ]:
            raise InvalidScheduleType(
                f"Internal Heat Gain {self.name}, schedule type must be 'Dimensionless' or 'Percent': {value.schedule_type}"
            )
        self._schedule = value

    @property
    def fraction_to_zone(self):
        return self._fraction_to_zone

    @fraction_to_zone.setter
    def fraction_to_zone(self, value):
        if value < 0 or value > 1.:
            raise ValueError(f"Internal Heat Gain {self.name}, fraction to zone outside range [0.,1.]: {value}")
        self._fraction_to_zone = float(value)

    @property
    def fraction_latent(self):
        return self._fraction_latent

    @fraction_latent.setter
    def fraction_latent(self, value):
        if value < 0 or value > 1.:
            raise ValueError(f"Internal Heat Gain {self.name}, fraction latent outside range [0.,1.]: {value}")
        self._fraction_latent = float(value)

    @property
    def fraction_radiant(self):
        return self._fraction_radiant

    @fraction_radiant.setter
    def fraction_radiant(self, value):
        if value < 0 or value > 1.:
            raise ValueError(f"Internal Heat Gain {self.name}, fraction radiant outside range [0.,1.]: {value}")
        try:
            if abs(value + self._fraction_convective - 1.) > 1e-5:
                raise ConvectiveRadiantFractionError(
                    f"Internal Heat Gain {self.name}, radiant/convective fraction sum not 1. Radiant = {value}, convective = {self.fraction_convective}"
                )
        except AttributeError:
            # This is just to avoid the check if self.fraction_radiant doesn't exist
            pass

        self._fraction_radiant = float(value)

    @property
    def fraction_convective(self):
        return self._fraction_convective

    @fraction_convective.setter
    def fraction_convective(self, value):
        if value < 0 or value > 1.:
            raise ValueError(f"Internal Heat Gain {self.name}, fraction convective outside range [0.,1.]: {value}")
        try:
            if abs(value + self.fraction_radiant - 1.) > 1e-5:
                raise ConvectiveRadiantFractionError(
                    f"Internal Heat Gain {self.name}, radiant/convective fraction sum not 1. Convective = {value}, radiant = {self.fraction_radiant}"
                )
        except AttributeError:
            # This is just to avoid the check if self.fraction_radiant doesn't exist
            pass
        self._fraction_convective = float(value)

    def get_convective_load(self, *args, **kwarg) -> np.array:
        """Just an empty method to raise an NotImplementedError Exception. This way any inherited class implements it

        Parameters
        ----------
        args
        kwarg

        """
        raise NotImplementedError(
            f"""
You must override the get_convective_load method for each class inherited from InternalLoad
Return value must be a np.array
"""
        )

    def get_radiant_load(self, *args, **kwarg) -> np.array:
        """Just an empty method to raise an NotImplementedError Exception. This way any inherited class implements it

        Parameters
        ----------
        args
        kwarg

        """
        raise NotImplementedError(
            f"""
You must override the get_radiant_load method for each class inherited from InternalLoad
Return value must be a np.array
"""
        )

    def get_latent_load(self, *args, **kwarg) -> np.array:
        """Just an empty method to raise an NotImplementedError Exception. This way any inherited class implements it

        Parameters
        ----------
        args
        kwarg

        """
        raise NotImplementedError(
            f"""
You must override the get_latent_load method for each class inherited from InternalLoad
Return value must be a np.array
"""
        )

    def get_loads(self, *args, **kwargs) -> list:
        """Return the convective, radiant, latent, electric load (numpy.array)
        If the calculation method is specific (W/m2 or px/m2) the area must be passed as kwarg (example area=12.5)

        Parameters
        ----------
        area : float
            Area in m2. pass it as kwarg: load_obj.get_loads(area = 12.5)

        Parameters
        ----------
        tuple
            [numpy.array, numpy.array, numpy.array, numpy.array]
            the schedules: convective [W], radiant [W], vapour [kg_vap/s], electric [W]

        """
        if "area" not in kwargs.keys():
            area = None
        else:
            area = kwargs['area']
        conv = self.get_convective_load(area=area)
        rad = self.get_radiant_load(area=area)
        lat = self.get_latent_load(area=area)
        el = self.get_electric_load(area=area)
        return conv, rad, lat, el


class People(InternalLoad):
    def __init__(
            self,
            name: str,
            nominal_value: float,
            unit: str,
            schedule: Schedule,
            fraction_latent: float = 0.55,
            fraction_radiant: float = 0.3,
            fraction_convective: float = 0.7,
            metabolic_rate: float = 110,
            tag: str = None,
    ):
        f"""Inherited from InternalLoad class. CHecks the values throughout propeties setter methods
        The sum of radiant and convective fraction must be 1
        
        Parameters
        ----------
        name : str
            name
        nominal_value : float
            the value to be multiplied by the schedule
        unit : str
            define the unit from the list {internal_loads_prop["people"]["unit"]}
        schedule : eureca_building.Schedule
            Schedule object
        fraction_latent : float, default 0.55
            latent fraction (between 0 and 1)
        fraction_radiant : float, default 0.3
            radiant fraction of the sensible part (between 0 and 1)
        fraction_convective : float, default 0.7
            convective fraction of the sensible part (between 0 and 1)
        metabolic_rate : float, default 110
            Metabolic rate [W/px]
        tag : str, default None
            Tag string to define the type of load
        """
        super().__init__(name, nominal_value, schedule, tag=tag)
        self.unit = unit
        self.fraction_latent = fraction_latent
        self.fraction_radiant = fraction_radiant
        self.fraction_convective = fraction_convective
        self.metabolic_rate = metabolic_rate

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        if not isinstance(value, str):
            raise TypeError(f"People load {self.name}, type is not a str: {value}")
        if value not in internal_loads_prop["people"]["unit"]:
            raise InvalidHeatGainUnit(
                f"People load {self.name}, unit not in: {internal_loads_prop['people']['unit']}\n{value}"
            )
        if value in ["W/m2", "px/m2", ]:
            self._calculation_method = "floor_area"
        elif value in ["W", "px", ]:
            self._calculation_method = "absolute"
        self._unit = value

    @property
    def metabolic_rate(self):
        return self._metabolic_rate

    @metabolic_rate.setter
    def metabolic_rate(self, value):
        try:
            value = float(value)
        except ValueError:
            raise ValueError(f"People load {self.name}, metabolic_rate is not a float: {value}")
        if value < 0.:
            raise ValueError(
                f"People load {self.name}, negative metabolic rate: {value}"
            )
        if value > 250.:
            logging.warning(
                f"People load {self.name}, metabolic rate over 250 W/px: {value}"
            )
        self._metabolic_rate = value

    def _get_absolute_value_nominal(self, area=None):
        """Memorizes the occupancy nominal value in W

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provided if the load is specific

        """
        if self._calculation_method == "floor_area":
            if area == None:
                raise AreaNotProvided(
                    f"Internal Heat Gain {self.name}, specific load but area not provided."
                )
            area = float(area)
        else:
            area = 1.
        if self.unit in ["px", "px/m2", ]:
            px_w_converter = self.metabolic_rate
        else:
            px_w_converter = 1.

        # This value is in W
        self.nominal_value_absolute = area * px_w_converter * self.nominal_value

    def get_convective_load(self, area=None):
        """Returns the convective load numpy.array
        If the calculation method is specific (W/m2 or px/m2) the area must be passed

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provided if the load is specific

        Returns
        ----------
        numpy.array
            the schedule [W]

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area)
        convective_fraction = (1 - self.fraction_latent) * self.fraction_convective
        return convective_fraction * self.nominal_value_absolute * self.schedule.schedule

    def get_radiant_load(self, area=None):
        """Returns the radiant load numpy.array
        If the calculation method is specific (W/m2 or px/m2) the area must be passed

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provide if the load is specific

        Returns
        ----------
        numpy.array
            the schedule [W]

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area)
        radiant_fraction = (1 - self.fraction_latent) * self.fraction_radiant
        return radiant_fraction * self.nominal_value_absolute * self.schedule.schedule

    def get_latent_load(self, area=None):
        """Return the latent load numpy.array
        If the calculation method is specific (W/m2 or px/m2) the area must be passed

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provided if the load is specific

        Returns
        ----------
        numpy.array
            the schedule [W]

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area)
        vapour_nominal_value_kg_s = self.fraction_latent * self.nominal_value_absolute / vapour_properties[
            'latent_heat']
        return self.schedule.schedule * vapour_nominal_value_kg_s

    def get_electric_load(self, area=None):
        """Return the electric consumption load numpy.array
        If the calculation method is specific (W/m2 or px/m2) the area must be passed

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provided if the load is specific

        Returns
        ----------
        numpy.array
            the schedule [W]

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area)
        return self.nominal_value_absolute * self.schedule.schedule * 0


class ElectricLoad(InternalLoad):
    def __init__(
            self,
            name: str,
            nominal_value: float,
            unit: str,
            schedule: Schedule,
            fraction_to_zone: float = 1.,
            fraction_radiant: float = 0.3,
            fraction_convective: float = 0.7,
            number_of_people: float = None,
            tag: str = None,
    ):
        f"""Inherited from InternalLoad class. Uses properties set methods to check types 
        The sum of radiant and convective fraction must be 1

        Parameters
        ----------
        name : str
            name
        nominal_value : float
            the value to be multiplied by the schedule
        unit : str
            define the unit from the list {internal_loads_prop["people"]["unit"]}
        schedule : eureca_building.Schedule
            Schedule object
        fraction_to_zone : float, default 1.
            fraction that is actually counted as Heat Load for the zone (between 0 and 1)
        fraction_radiant : float, default 0.3
            radiant fraction of the sensible part (between 0 and 1)
        fraction_convective : float, default 0.7
            convective fraction of the sensible part (between 0 and 1)
        number_of_people : float, default None
            if the unit is W/px then this number must be passed
        tag : str, default None
            Tag string to define the type of load
            
        """
        super().__init__(name, nominal_value, schedule, tag=tag)
        self.unit = unit
        self.fraction_to_zone = fraction_to_zone
        self.fraction_radiant = fraction_radiant
        self.fraction_convective = fraction_convective
        self.number_of_people = number_of_people

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        if not isinstance(value, str):
            raise TypeError(f"ElectricLoad load {self.name}, type is not a str: {value}")
        if value not in internal_loads_prop["electric"]["unit"]:
            raise InvalidHeatGainUnit(
                f"ElectricLoad load {self.name}, unit not in: {internal_loads_prop['people']['unit']}\n{value}"
            )
        if value in ["W/m2"]:
            self._calculation_method = "floor_area"
        elif value in ["W/px"]:
            self._calculation_method = "people"
        elif value in ["W"]:
            self._calculation_method = "absolute"
        self._unit = value

    @property
    def number_of_people(self):
        return self._number_of_people

    @number_of_people.setter
    def number_of_people(self, value):
        if value is not None:
            try:
                value = value
            except ValueError:
                raise ValueError(f"ElectricLoad load {self.name}, number_of_people is not a float: {value}")
            if value < 0.:
                raise ValueError(
                    f"ElectricLoad load {self.name}, negative number_of_people: {value}"
                )
        self._number_of_people = value

    def _get_absolute_value_nominal(self, area=None):
        """Memorizes the electric load nominal value in W

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provided if the load is specific

        """
        if self._calculation_method == "floor_area":
            if area == None:
                raise AreaNotProvided(
                    f"Internal Heat Gain {self.name}, specific load but area not provided."
                )
            area = float(area)
        else:
            area = 1.
        if self._calculation_method == "people":
            if self.number_of_people == None:
                raise PeopleNotProvided(
                    f"Internal Heat Gain {self.name}, people calculation but number of people not provided."
                )
            number_of_people = float(self.number_of_people)
        else:
            number_of_people = 1.

        # print(area)
        # print(number_of_people)
        # This value is in W
        self.nominal_value_absolute = number_of_people * area * self.nominal_value

    def get_convective_load(self, area=None):
        """Returns the convective load numpy.array
        If the calculation method is specific (W/m2 or px/m2) the area must be passed

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provided if the load is specific

        Returns
        ----------
        numpy.array
            the schedule [W]

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area)
        convective_fraction = self.fraction_to_zone * self.fraction_convective
        return convective_fraction * self.nominal_value_absolute * self.schedule.schedule

    def get_radiant_load(self, area=None):
        """Returns the radiant load numpy.array
        If the calculation method is specific (W/m2 or px/m2) the area must be passed

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provide if the load is specific

        Returns
        ----------
        numpy.array
            the schedule [W]

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area)
        radiant_fraction = self.fraction_to_zone * self.fraction_radiant
        return radiant_fraction * self.nominal_value_absolute * self.schedule.schedule

    def get_latent_load(self, area=None):
        """Return the latent load numpy.array
        If the calculation method is specific (W/m2 or px/m2) the area must be passed

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provided if the load is specific

        Returns
        ----------
        numpy.array
            the schedule [W]

        """
        return self.schedule.schedule * 0

    def get_electric_load(self, area=None):
        """Return the electric consumption load numpy.array
        If the calculation method is specific (W/m2 or px/m2) the area must be passed

        Parameters
        ----------
        area : float, default None
            Area of the zone in [m2]: must be provided if the load is specific

        Returns
        ----------
        numpy.array
            the schedule [W]

        """
        # try:
        #     self.nominal_value_absolute
        # except AttributeError:
        self._get_absolute_value_nominal(area=area)
        return self.nominal_value_absolute * self.schedule.schedule


class Lights(ElectricLoad):
    def __init__(
            self,
            name: str,
            nominal_value: float,
            unit: str,
            schedule: Schedule,
            fraction_to_zone: float = 1.,
            fraction_radiant: float = 0.3,
            fraction_convective: float = 0.7,
            number_of_people: float = None,
            tag: str = None,
    ):
        f"""Inherited from ElectricLoad class. Uses properties set methods to check types 
        This is just a wrapper for an EletricLoad object as they are similar

        Parameters
        ----------
        name : str
            name
        nominal_value : float
            the value to be multiplied by the schedule
        unit : str
            define the unit from the list {internal_loads_prop["people"]["unit"]}
        schedule : Schedule
            Schedule object
        fraction_to_zone : float, default 1.
            fraction that is actually counted as Heat Load for the zone (between 0 and 1)
        fraction_radiant : float, default 0.3
            radiant fraction of the sensible part (between 0 and 1)
        fraction_convective : float, default 0.7
            convective fraction of the sensible part (between 0 and 1)
        number_of_people : float, default None
            if the unit is W/px then this number must be passed
        tag : str, default None
            Tag string to define the type of load
        """
        super().__init__(
            name,
            nominal_value,
            unit,
            schedule,
            fraction_to_zone=fraction_to_zone,
            fraction_radiant=fraction_radiant,
            fraction_convective=fraction_convective,
            number_of_people=number_of_people,
            tag=tag,
        )
