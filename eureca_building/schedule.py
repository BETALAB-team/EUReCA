"""
This module includes functions to model any schedule
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging

import numpy as np
# import matplotlib.pyplot as plt

from eureca_building.schedule_properties import schedule_types
from eureca_building.config import CONFIG
from eureca_building.exceptions import (
    InvalidScheduleType,
    ScheduleOutsideBoundaryCondition,
    InvalidScheduleDimension,
    ScheduleLengthNotConsistent,
)


class Schedule:
    """
    Represents time-varying schedules for energy modeling (e.g., temperature, availability).
    
    Supports validation of bounds, unit types, and generation from constant or daily profiles.
    """

    def __init__(
            self,
            name: str,
            schedule_type: str,
            schedule: np.array,
            lower_limit=None,
            upper_limit=None,
    ):
        f"""Schedule Constructor and check the input values and types
        
        Parameters
        ----------
        name : str
            name
        schedule_type : str
            type of the {schedule_types["unit_type"]}
        schedule : numpy.array
            the schedule array, length equal to 8760 time the number of time steps per hour
        upper_limit : float, default None
            upper limit to check schedule validity
        lower_limit: float, default None
            upper limit to check schedule validity
        """
        self.name = str(name)
        self.schedule_type = schedule_type
        self._lower_limit = lower_limit
        self._upper_limit = upper_limit
        self.schedule = schedule

    @property
    def schedule_type(self):
        return self._schedule_type

    @schedule_type.setter
    def schedule_type(self, value):
        if not isinstance(value, str):
            raise TypeError(f"Schedule {self.name}, type is not a str: {value}")
        if value not in schedule_types["unit_type"]:
            raise InvalidScheduleType(
                f"Schedule {self.name}, type not in: {schedule_types['unit_type']}\n{value}"
            )
        self._schedule_type = value

    @property
    def _lower_limit(self):
        return self.__lower_limit

    @_lower_limit.setter
    def _lower_limit(self, value):
        if value is not None:
            if not isinstance(value, float) and not isinstance(value, int):
                raise TypeError(f"Schedule {self.name}, lower limit is not a number: {value}")
            self.__lower_limit = value
        else:
            self.__lower_limit = -1e20
        if self.schedule_type == "Percent":
            if value is not None:
                logging.warning(
                    f"""
Schedule {self.name}, the schedule is a percentage schedule but a lower limit was set.
Lower limit set to 0."""
                )
            self.__lower_limit = 0.

    @property
    def _upper_limit(self):
        return self.__upper_limit

    @_upper_limit.setter
    def _upper_limit(self, value):
        if value is not None:
            if not isinstance(value, float) and not isinstance(value, int):
                raise TypeError(f"Schedule {self.name}, upper limit is not a number: {value}")
            self.__upper_limit = value
        else:
            self.__upper_limit = 1e20
        if self.schedule_type == "Percent":
            if value is not None:
                logging.warning(
                    f"""
Schedule {self.name}, the schedule is a percentage schedule but a upper limit was set.
Lower limit set to 1."""
                )
            self.__upper_limit = 1.

    @property
    def schedule(self):
        return self._schedule

    @schedule.setter
    def schedule(self, _value):
        try:
            value = np.array(_value, dtype=float)
        except ValueError:
            raise ValueError(f"Schedule {self.name}, non-numeric values in the schedule")
        if value.ndim > 1:
            raise InvalidScheduleDimension(f"Schedule {self.name}, schedule dimension higher than 1: {value.ndim}")
        if np.any(np.greater(value, self._upper_limit)):
            raise ScheduleOutsideBoundaryCondition(
                f"Schedule {self.name}, there is a value above the upper limit: upper limit {self._upper_limit}"
            )
        if np.any(np.less(value, self._lower_limit)):
            raise ScheduleOutsideBoundaryCondition(
                f"Schedule {self.name}, there is a value below the lower limit: lower limit {self._lower_limit}"
            )
        if len(value) != CONFIG.number_of_time_steps_year:
            raise ScheduleLengthNotConsistent(
                f"""
Schedule {self.name}: the length of the schedule is not consistent 
with the number of time steps provided. 
Schedule length : {len(value)}
Number of time steps: {CONFIG.number_of_time_steps_year}
                """
            )

        self._schedule = value

    # def plot(self):
        # plt.plot(self.schedule)
        # plt.title(f'Schedule: {self.name}')

    @classmethod
    def from_daily_schedule(
            cls,
            name: str,
            schedule_type: str,
            schedule_week_day: np.array,
            schedule_saturday: np.array,
            schedule_sunday: np.array,
            schedule_holiday: np.array,
            lower_limit=None,
            upper_limit=None,
            starting_day: int = 0,
            holidays: tuple = (),
    ):
        f"""Class method. This method allows to create a simulation schedule using daily profiles.

        Parameters
        ----------
        name : str
            name
        schedule_type :
            type of the {schedule_types["unit_type"]}
        schedule_week_day : numpy.array
            week_day schedule (length 24 * n_ts)
        schedule_saturday : numpy.array
            saturday  schedule (length 24 * n_ts)
        schedule_sunday : numpy.array
            sunday  schedule (length 24 * n_ts)
        schedule_holiday : numpy.array
            holiday schedule (length 24 * n_ts)
        upper_limit : float, default None
            upper limit to check schedule validity
        lower_limit : float, default None
            upper limit to check schedule validity
        holidays : tuple
            tuple of holidays (with int from 0 to 364)
        starting_day : int
            day to start the year (0 monday, 1 tuesday, ... 6 sunday)

        Returns
        ----------
        eureca_building.schedule.Schedule
        """
        try:
            holidays = tuple(holidays)
        except ValueError:
            raise TypeError(
                f"Schedule {name}, holidays is not a tuple: holidays = {holidays}"
            )
        for i in holidays:
            if not isinstance(i, int) or i > 364:
                raise TypeError(
                    f"Schedule {name}, holidays list must contain only int less from 0 to 364: holidays = {holidays}"
                )

        if not isinstance(starting_day, int) or starting_day > 6:
            raise TypeError(
                f"Schedule {name}, starting day must be an int between 0 and 6: starting day = {starting_day}"
            )


        week = np.hstack([np.tile(schedule_week_day, 5),schedule_saturday,schedule_sunday])
        year = np.tile(week,54)
        year_net = year[starting_day * 24 * CONFIG.ts_per_hour : (starting_day * 24 + 8760) * CONFIG.ts_per_hour]
        if CONFIG.ts_per_hour > 1:
            year_net = year_net[:(1-CONFIG.ts_per_hour)]
        for day in holidays:
            year_net[day * 24 * CONFIG.ts_per_hour: (day + 1) * 24 * CONFIG.ts_per_hour] = schedule_holiday

        schedule = cls(
            name = name,
            schedule_type = schedule_type,
            schedule = year_net,
            upper_limit = upper_limit,
            lower_limit = lower_limit,
        )

        return schedule


    @classmethod
    def from_constant_value(
            cls,
            name: str,
            schedule_type: str,
            value: float,
            lower_limit=None,
            upper_limit=None,
    ):
        f"""Class method. This method allows to create a simulation schedule from a constant value.

        Parameters
        ----------
        name : str
            Name
        schedule_type : str
            type of the {schedule_types["unit_type"]}
        value : float
            the schedule value (might be of different units)
        upper_limit : float, default None
            upper limit to check schedule validity
        lower_limit: float, default None
            upper limit to check schedule validity

        Returns
        -------
        eureca_building.schedule.Schedule
        """

        sched = np.array([value] * 24 * 365 * CONFIG.ts_per_hour)
        if CONFIG.ts_per_hour > 1:
            sched = sched[:(1-CONFIG.ts_per_hour)]
        schedule = cls(
            name=name,
            schedule_type=schedule_type,
            schedule=sched,
            upper_limit=upper_limit,
            lower_limit=lower_limit,
        )

        return schedule
