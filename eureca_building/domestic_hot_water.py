"""
This module includes functions to model natural ventilation and infiltration
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging
import math

import numpy as np

from eureca_building.config import CONFIG
from eureca_building.schedule_properties import domestic_hot_water_prop
from eureca_building.fluids_properties import water_properties
from eureca_building.schedule import Schedule
from eureca_building.exceptions import InvalidScheduleType

def distrEventi(n,x,pdf):
    y_guess=np.random.rand(int(n))
    cdf=np.cumsum(pdf)
    [cdf,index]=np.unique(cdf,1)

    x_event=np.interp(y_guess,cdf,x[index])
    time_event=np.round(x_event)
    return time_event

def dhw_calc_calculation(volume_unit, numunits, time_step):

    total_days = domestic_hot_water_prop["total_days"]
    nuses = domestic_hot_water_prop["nuses"]
    vol_aver_drawoff = domestic_hot_water_prop["vol_aver_drawoff"]
    vol_desv_drawoff = domestic_hot_water_prop["vol_desv_drawoff"]
    time_aver_drawoff = domestic_hot_water_prop["time_aver_drawoff"]
    proportion_event = domestic_hot_water_prop["proportion_event"]
    dist = domestic_hot_water_prop["dist"]

    vol_total_drawoff = vol_aver_drawoff * time_aver_drawoff

    vol_total_drawoff_use = np.zeros((nuses))

    for j in range(nuses):
        vol_total_drawoff_use[j] = volume_unit * proportion_event[j]

    draw_offs = np.zeros((nuses))
    draw_offs_dec = np.zeros((nuses))

    for j in range(nuses):
        draw_offs_dec[j] = np.abs(vol_total_drawoff_use[j] / vol_aver_drawoff[j] - int(
            vol_total_drawoff_use[j] / vol_aver_drawoff[j]))
        if draw_offs_dec[j] >= 0.5:
            draw_offs[j] = math.ceil(vol_total_drawoff_use[j] / vol_aver_drawoff[j])
        else:
            draw_offs[j] = np.round(vol_total_drawoff_use[j] / vol_aver_drawoff[j])

    time_steps_hour = 60 / time_step
    time_steps_day = 24 * time_steps_hour

    array_time_step = np.zeros(4)
    for i in range(4):
        array_time_step[i] = time_step

    n_max = array_time_step / time_aver_drawoff
    for i in range(len(n_max)):
        n_max[i] = math.ceil(n_max[i])

    Volume_use_daily_array = np.zeros((total_days, int(time_steps_day)))
    Volume_use_arrayb0 = np.zeros((int(total_days * time_steps_day)))
    dist_use = np.zeros((24, 12))
    Volume_use_array = np.zeros((int(time_steps_day * total_days), nuses))
    Volume_use_sum = np.zeros((1, nuses))
    Volume_aver_drawoff_final = np.zeros(nuses)
    Volume_desv_drawoff_final = np.zeros((1, nuses))
    Volume_use_time = np.zeros(np.int(time_steps_day))

    Volume_use_unit = np.zeros((int(total_days * time_steps_day), numunits))

    for units0 in range(numunits):

        for use in range(nuses):

            dist_use = np.repeat(dist[:, use], 12, axis=0)
            dist_use_t = np.reshape(dist_use, (288, 1)) / time_steps_hour

            vol_aver_drawoff_use = vol_aver_drawoff[use]
            vol_desv_drawoff_use = vol_desv_drawoff[use]
            time_aver_drawoff_use = time_aver_drawoff[use]
            draw_offs_use = int(draw_offs[use])
            n_max1 = n_max[use]

            for day in range(total_days):

                time_event = distrEventi(draw_offs_use, np.arange(0, time_steps_day), dist_use_t)

                if use == 0:
                    m = vol_aver_drawoff_use
                    v = vol_desv_drawoff_use
                    mu = np.log((m ** 2) / np.sqrt(v + m ** 2))
                    sigma = np.sqrt(np.log(v / (m ** 2) + 1))
                    Flow_rated = np.random.lognormal(mu, sigma, int(draw_offs_use))

                else:
                    Flow_rated = np.random.normal(vol_aver_drawoff_use, vol_desv_drawoff_use, int(draw_offs_use))

                Volume_use = np.array(np.double(np.abs(Flow_rated * time_aver_drawoff_use)))

                for i in range(int(time_steps_day)):

                    index = np.where(time_event == i)

                    if len(index) > n_max1:
                        index1 = index[0, int(n_max1)]
                    else:
                        index1 = index

                    if draw_offs_use == 1:
                        if index1 == 0:
                            Volume_use_time = Volume_use
                    else:
                        Volume_use_time = Volume_use[index1]

                    Volume_use_daily_array[(day, i)] = np.sum(Volume_use_time) / time_aver_drawoff_use

            Volume_use_resh1 = np.reshape(Volume_use_daily_array, (1, int(total_days * time_steps_day)))
            Volume_use_array[:, use] = Volume_use_resh1

            Volume_use_sum[0, use] = np.sum(Volume_use_array[:, use])

        Volume_use_unit[:, units0] = np.sum(Volume_use_array, axis=1)

        # if numunits>1:
        Volume_use_arrayb0 = np.sum(Volume_use_unit, axis=1)
        # else:
        #     Volume_use_arrayb0=Volume_use_unit[:,0]

    Volume_totalb1 = np.sum(Volume_use_arrayb0)

    number_units = numunits

    volume_profile = Volume_use_arrayb0
    total_volume = Volume_totalb1

    Volume_meanb1 = ((Volume_totalb1) / numunits) / total_days
    return volume_profile

class DomesticHotWater:

    def __init__(
            self,
            name: str,
            calculation_method: str,
            unit = None,
            schedule = None,
    ):
        self.name = name
        self.calculation_method = calculation_method
        self.unit = unit
        self.schedule = schedule

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

    def get_dhw_yearly_mass_flow_rate(self, area, number_of_units, weather):

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

            if self.calculation_method == "UNI-TS 11300-2":
                volume = np.ones(CONFIG.number_of_time_steps_year) * Vw_single * number_of_units * 365 / CONFIG.number_of_time_steps_year / CONFIG.time_step # To convert from m3/ts to m3/s

            if self.calculation_method == "DHW calc":
                # Broken
                # TODO: fix stochastic calculation
                volume = dhw_calc_calculation(Vw_single, number_of_units, CONFIG.time_step / 60)[:CONFIG.number_of_time_steps_year] / 1000 / CONFIG.time_step # to m3/ts to m3/s
                # L or m3?

        Cw = water_properties["specific_heat"]      # [J/(kg K)]
        DTw = domestic_hot_water_prop["target temperature [°C]"] - weather.general_data["average_out_air_db_temperature"]       # [°C]
        rho = water_properties["density"] # [kg/m3]

        demand = volume *  rho * Cw * DTw # W

        return volume, demand







