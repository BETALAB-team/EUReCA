"""
This module includes functions to model setpoint of the thermal zone
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging

import numpy as np

from eureca_building.schedule_properties import setpoint_prop
from eureca_building.schedule import Schedule
from eureca_building.fluids_properties import vapour_properties

from eureca_building.exceptions import (
    SetpointTypeNotAllowed,
)


# TODO: implement interface for methods to be implemented
# https://realpython.com/python-interface/

class Setpoint:
    """
    Internal Gain Class
    """

    def __init__(
            self,
            name: str,
            setpoint_type: str,
            tag: str = None,
    ):
        f"""
        Parent class for some inherited SetpointObjects

        Args:
            name: str
                name
            setpoint_type: float
                type from : {setpoint_prop.keys()}
            tag: str
                a tag to define the type of internal load
        """
        self.name = name
        self.setpoint_type = setpoint_type
        self.tag = tag

    @property
    def setpoint_type(self):
        return self._setpoint_type

    @setpoint_type.setter
    def setpoint_type(self, value):
        try:
            value = str(value)
        except ValueError:
            raise ValueError(f"Setpoint {self.name}, schedule_type is not a str: {value}")
        if value not in setpoint_prop.keys():
            raise SetpointTypeNotAllowed(
                f"Setpoint {self.name}, {value} schedule_type not allowed. Chose from: {setpoint_prop.keys()}")
        self._setpoint_type = value

    def get_convective_load(self, *args, **kwarg) -> np.array:
        raise NotImplementedError(
            f"""
You must override the get_convective_load method for each class inherited from InternalLoad
Return value must be a np.array
"""
        )


class SetpointDualBand(Setpoint):
    def __init__(self,
                 name: str,
                 setpoint_type: str,
                 schedule_lower: Schedule,
                 schedule_upper: Schedule,
                 tag: str = None, ):
        super().__init__(name, setpoint_type, tag)
        self.schedule_lower = schedule_lower
        self.schedule_upper = schedule_upper
        # Check if setpoints schedules are not intersecting
        if np.any(np.greater_equal(self.schedule_lower.schedule, self.schedule_upper.schedule)):
            raise ValueError(f"Setpoint {self.name}, lower schedule higher than upper schedule.")

    @property
    def schedule_lower(self):
        return self._schedule_lower

    @schedule_lower.setter
    def schedule_lower(self, value):
        if not isinstance(value, Schedule):
            raise ValueError(f"Setpoint {self.name}, lower schedule type not Schedule: {type(value)}")
        if np.any(np.less(value.schedule, setpoint_prop[self.setpoint_type]['limit'][0])):
            logging.warning(
                f"Setpoint {self.name}, lower schedule goes below {setpoint_prop[self.setpoint_type]['limit'][0]} {setpoint_prop[self.setpoint_type]['unit']}")
        self._schedule_lower = value

    @property
    def schedule_upper(self):
        return self._schedule_upper

    @schedule_upper.setter
    def schedule_upper(self, value):
        if not isinstance(value, Schedule):
            raise ValueError(f"Setpoint {self.name}, lower schedule type not Schedule: {type(value)}")
        if np.any(np.greater(value.schedule, setpoint_prop[self.setpoint_type]['limit'][1])):
            logging.warning(
                f"Setpoint {self.name}, lower schedule goes above {setpoint_prop[self.setpoint_type]['limit'][1]} {setpoint_prop[self.setpoint_type]['unit']}")
        self._schedule_upper = value
