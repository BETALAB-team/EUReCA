"""
This module includes functions to model Domestic Hot Water consumptions
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging
import math

import numpy as np
import pandas as pd

from eureca_building.config import CONFIG
from eureca_building.schedule_properties import domestic_hot_water_prop
from eureca_building.fluids_properties import water_properties
from eureca_building.schedule import Schedule
from eureca_building.exceptions import InvalidScheduleType

def _event_distribution(number_of_daily_events,daily_vector_distibution, number_of_days):
    """
    Parameters
    ----------
    n
    x
    pdf vettore distribuzione

    Returns
    -------

    """
    y_guess = np.random.rand(number_of_days, int(number_of_daily_events))
    cdf = np.cumsum(daily_vector_distibution)
    x_event = np.interp(y_guess,cdf,np.arange(len(cdf)))
    # get randomly som time events from the cdf and y_guess selected randomly
    return np.round(x_event).astype(int)

def dhw_calc_calculation(volume_unit, numunits):
    """

    Parameters
    ----------
    volume_unit
    numunits
    time_step

    Returns
    -------

    """
    # The DHW calc is always run with a 5 min time step. At the end it is resapled to the CONFIG time step
    volume_unit = volume_unit * 1000 # m3 to l
    time_step = 5 # min
    time_steps_hour = int(60 / time_step)
    time_steps_day = int(24 * time_steps_hour)

    total_consumed_volume_dw = np.zeros([365*24*time_steps_hour, numunits])

    for unit in range(numunits):


        total_days = domestic_hot_water_prop["total_days"]

        consumed_volume = np.zeros([total_days, time_steps_day, len(domestic_hot_water_prop["DHWcalc_uses"].items())])

        use_number = 0
        for use, values in domestic_hot_water_prop["DHWcalc_uses"].items():

            average_time_drawoff = values["average_time"]
            # l/d split by drawoff type
            daily_vol_total_drawoff_use = volume_unit * values["percentage_to_total_consumption"]
            daily_drawoff_on_event = np.round(daily_vol_total_drawoff_use / values["mean"]).astype(int)
            # Maximum number of on events in a timestep for each drawoff
            # n_max_time_step_on_events = np.ceil(time_step / average_time_drawoff)

            # This does a resampling with a linear interpolation method  to keep 1 as sum of each column
            temporal_dist = np.interp(np.arange(0, 24, 1 / time_steps_hour), np.arange(0, 24),
                                      values["temporal_distribution"]) / time_steps_hour

            # array with the time steps when the use is on [number of daily on events, number of days]
            temporal_dist_year = _event_distribution(daily_drawoff_on_event, temporal_dist, total_days)
            if values["volume_pdf"] == "lognormal":
                m = values["mean"]
                v = values["std"]
                mu = np.log((m ** 2) / np.sqrt(v + m ** 2))
                sigma = np.sqrt(np.log(v / (m ** 2) + 1))
                flow_rate = np.random.lognormal(mu, sigma, [total_days, daily_drawoff_on_event])  # l/min
            elif values["volume_pdf"] == "normal":
                flow_rate = np.random.normal(values["mean"], values["std"], [total_days, daily_drawoff_on_event])
            else:
                raise ValueError(
                    f'DHW calculation. The volume flow rate probabilty distribution function is not allowed. PDF: {values["volume_pdf"]}. Aloowed PDFs: [lognormal, normal]')
            volume_use = np.abs(flow_rate)  # l ad accensione
            # To avoid negative consumptions (possible with normal dist)
            volume_use[volume_use < 0] = 0.

            for ts in range(daily_drawoff_on_event):
                consumed_volume[np.arange(total_days), temporal_dist_year[:, ts], use_number] = volume_use[:, ts]

            # The rescale needed because otherwise the rounding process provoke an underestimation
            consumed_volume[:, :, use_number] = consumed_volume[:, :, use_number]/(consumed_volume[:, :, use_number].sum() / 365) * daily_vol_total_drawoff_use

            use_number += 1

        total_consumed_volume_dw[:,unit] = consumed_volume.sum(axis=2).reshape(365*24*time_steps_hour)

    total_consumed_volume = total_consumed_volume_dw.sum(axis=1)
    total_consumed_volume_rs = pd.Series(total_consumed_volume, index = pd.date_range(start='1/1/2018 00:00', periods=8760*time_steps_hour, freq = f"{time_step}min")).resample(f"{CONFIG.time_step}S").sum().values

    return total_consumed_volume_rs # [CONFIG.start_time_step:CONFIG.final_time_step]

class DomesticHotWater:
    """DomesticHotWater object
    Class to manage all the calculations involved in the Domestic Hot Water consumption
    """

    def __init__(
            self,
            name: str,
            calculation_method: str,
            unit = None,
            schedule = None,
            n_of_occupants = None,
    ):
        f"""Constructor for DomesticHotWater. Memorizes the attributes anc checks them through properties setter

        Parameters
        ----------
        name : str
            name of the object
        calculation_method : str
            Calculation method, choose from {domestic_hot_water_prop['calculation_method']}
        unit : str
            Unit of the schedule, choose from {domestic_hot_water_prop['unit']}
        schedule : eureca_building.schedule.Schedule
            Schedule object, to be used in case the method is 'schedule'
        """


        self.name = name
        self.calculation_method = calculation_method
        self.unit = unit
        self.schedule = schedule
        self.n_of_people = n_of_occupants
        if self.calculation_method == "Number of occupants" and self.n_of_people == None:
            raise TypeError("DHW object: If 'Number of occupants' is selected as calculation method, a numeric n_of_occupants arg must be introduced")

    @property
    def calculation_method(self):
        return self._calculation_method

    @calculation_method.setter
    def calculation_method(self, value):
        if not isinstance(value, str):
            raise ValueError(f"Domestic hot water {self.name}, calculation method must be a str: {value}")
        if value not in domestic_hot_water_prop['calculation_method']:
            raise ValueError(f"Domestic hot water {self.name}, calculation method not valid: {value}: Choose from {domestic_hot_water_prop['calculatio_method']}")
        self._calculation_method = value

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        if self._calculation_method == "Schedule":
            if not isinstance(value, str):
                raise TypeError(f"Domestic hot water {self.name}, unit is not a str: {value}")
            if value not in domestic_hot_water_prop["unit"]:
                raise TypeError(
                    f"Domestic Hot Water {self.name}, unit not in: {domestic_hot_water_prop['unit']}\n{value}"
                )
        self._unit = value

    @property
    def schedule(self):
        return self._schedule

    @schedule.setter
    def schedule(self, value):
        if self._calculation_method == "Schedule":
            if not isinstance(value, Schedule):
                raise ValueError(f"Domestic Hot Water {self.name}, schedule type not Schedule: {type(value)}.\nIf you chose Schedule calculation method you must provide a mass flow rate schedule")
            if value.schedule_type not in ["mass_flow_rate",]:
                raise InvalidScheduleType(
                    f"Domestic Hot Water {self.name}, schedule type must be 'mass_flow_rate': {value.schedule_type}"
                )
        self._schedule = value

    def get_dhw_yearly_mass_flow_rate(self, area, number_of_units, weather, n_of_people = None):
        """This function calculates the water and mass flow rate consumption, given the area of the building and the number of units (to be used when unit and/or method need them)

        Parameters
        ----------
        area : float
            Area of the building [m2]
        number_of_units : int
            Number of dwellings (for residential calculation done with UNI-TS 11300
        weather : eureca_building.weather.WeatherFile
            WeatherFile object

        Returns
        -------
        tuple
            tuple of numpy.arrays
            volume flow rate [m3/s], dhw heating demand [W]
        """
        if self.calculation_method == "Schedule":
            schedule = self.schedule.schedule
            if self.unit == "L/s":
                volume = schedule / 1000 # to converto to m3/s
            if self.unit == "L/(m2 h)":
                volume = schedule * area / 3600 / 1000 # to converto to m3/s

        else:
            # Calculation of demand and volume with UNI-TS 11300
            # Vw [lt/day] = a [lt/(m2 day)] * Af [m2] + b [lt/day]
            Af = area / number_of_units
            if Af <= 35:
                a = 0.0; b = 50.0
            elif 35 < Af <= 50:
                a = 2.667; b = -43.33
            elif 50 < Af <= 200:
                a = 1.067; b = 36.67
            elif Af > 200:
                a = 0.0; b = 250.0

            # Water Need [m3/day]
            Vw_single = (a * Af + b) / 1000

            if self.calculation_method == "Number of occupants":
                try:
                    self.n_of_people = float(self.n_of_people)
                except ValueError:
                    raise TypeError(f"DHW object: number of people is not numeric : {self.n_of_people}")
                Vw_single = self.n_of_people * 0.04 / number_of_units #  Water Need [m3/day] considering 40 L/px

            if self.calculation_method in ["UNI-TS 11300-2","Number of occupants"]:
                sched = pd.Series([0.500, 0.502, 0.504, 0.957, 0.984, 1.042, 1.102, 1.120, 1.126, 1.131, 1.133, 1.132,
                         1.133, 1.136, 1.133, 1.135, 1.134, 1.134, 1.135, 1.133, 1.122, 1.102, 0.972, 0.498]*365,
                                  index = pd.date_range(start="00:00 01/01/2023", freq='1h', periods = 365*24)).resample(f"{int(CONFIG.time_step/60)}min").ffill().values
                # sched = np.ones(CONFIG.number_of_time_steps_year)

                volume = sched * Vw_single * number_of_units * 365 / CONFIG.number_of_time_steps_year / CONFIG.time_step # To convert from m3/ts to m3/s

            elif self.calculation_method == "DHW calc":
                # TODO: fix stochastic calculation
                volume = dhw_calc_calculation(Vw_single, number_of_units)[:CONFIG.number_of_time_steps_year] / 1000 / CONFIG.time_step # to m3/ts to m3/s

        Cw = water_properties["specific_heat"]      # [J/(kg K)]
        DTw = domestic_hot_water_prop["target temperature [°C]"] - weather.general_data["average_out_air_db_temperature"]       # [°C]
        rho = water_properties["density"] # [kg/m3]

        demand = volume *  rho * Cw * DTw # W

        return volume, demand







