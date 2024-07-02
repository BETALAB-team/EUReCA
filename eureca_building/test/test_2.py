"""
Tests for schedules
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import os
from eureca_building.config import load_config

config_path = os.path.join('.', 'eureca_building', 'test', 'config.json')
load_config(config_path)
from eureca_building.config import CONFIG

import pytest
import numpy as np

from eureca_building.schedule import Schedule
from eureca_building.internal_load import InternalLoad, People, ElectricLoad, Lights
from eureca_building.exceptions import (
    InvalidScheduleType,
    ScheduleOutsideBoundaryCondition,
    InvalidScheduleDimension,
)


class TestSchedule:
    """
    This is a test class for the pytest module.
    It tests Schedules class and Internal Heat Gains
    """

    def test_schedule(self):
        # Standard schedule
        Schedule(
            "Temperature1",
            "temperature",
            np.random.rand(CONFIG.number_of_time_steps_year) + 10,
        )

    def test_schedule_2(self):
        # Standard schedule
        with pytest.raises(ScheduleOutsideBoundaryCondition):
            Schedule(
                "Temperature1",
                "temperature",
                np.random.rand(8760*4-3) + 10,
                upper_limit=5,
                lower_limit=-10,
            )

    def test_schedule_3(self):
        # Standard schedule
        with pytest.raises(InvalidScheduleDimension):
            Schedule(
                "Temperature1",
                "temperature",
                np.array([[3, 3, 3, 5]]) + 10,
                upper_limit=5,
                lower_limit=-10,
            )


class TestInternalHeatGains:

    def test_IHG_creation(self):
        # Standard IHG
        sched = Schedule(
            "Percent1",
            "percent",
            np.array([0.1, .2, .3, .5] * int(CONFIG.number_of_time_steps_year / 4+ 1))[:-CONFIG.ts_per_hour + 1],
            upper_limit=1.,
            lower_limit=0.,
        )

        InternalLoad(
            name='test_IHG',
            nominal_value=10.,
            schedule=sched,
        )

        People(
            name='test_IHG',
            unit='W/m2',
            nominal_value=10.,
            schedule=sched,
            fraction_latent=0.55,
            fraction_radiant=0.3,
            fraction_convective=0.7,
            metabolic_rate=110,
        )

    def test_people_schedules(self):
        # Standard IHG
        sched = Schedule(
            "Percent1",
            "percent",
            np.array([0.1, .2, .3, .5] * int(CONFIG.number_of_time_steps_year / 4+ 1))[:-CONFIG.ts_per_hour + 1],
        )

        people1 = People(
            name='test_IHG',
            unit='W/m2',
            nominal_value=10.,
            schedule=sched,
            fraction_latent=0.45,
            fraction_radiant=0.3,
            fraction_convective=0.7,
        )

        conv, rad, lat, el = people1.get_loads(area=2.)
        conv = conv[:4]
        rad = rad[:4]
        lat = lat[:4]
        assert (
                np.linalg.norm(conv - np.array([0.385, 0.77, 1.155, 1.925]) * 2) < 1e-5 and
                np.linalg.norm(rad - np.array([0.165, 0.33, 0.495, 0.825]) * 2) < 1e-5 and
                np.linalg.norm(lat - np.array([1.79928029e-07,
                                               3.59856058e-07,
                                               5.39784086e-07,
                                               8.99640144e-07]) * 2) < 1e-5
        )

    def test_people_schedules_2(self):
        # Standard IHG
        sched = Schedule(
            "Percent1",
            "percent",
            np.array([0.1, .2, .3, .5] * int(CONFIG.number_of_time_steps_year / 4+ 1))[:-CONFIG.ts_per_hour + 1],
        )

        people1 = People(
            name='test_IHG',
            unit='px/m2',
            nominal_value=5.,
            schedule=sched,
            fraction_latent=0.45,
            fraction_radiant=0.3,
            fraction_convective=0.7,
            metabolic_rate=150,
        )

        conv, rad, lat, el = people1.get_loads(area=2.)
        conv = conv[:4]
        rad = rad[:4]
        lat = lat[:4]
        assert (
                np.linalg.norm(conv - np.array([57.75,
                                                115.5,
                                                173.25,
                                                288.75,
                                                ])) < 1e-5 and
                np.linalg.norm(rad - np.array([24.75,
                                               49.5,
                                               74.25,
                                               123.75,
                                               ])) < 1e-5 and
                np.linalg.norm(lat - np.array([2.69892E-05,
                                               5.39784E-05,
                                               8.09676E-05,
                                               0.000134946,
                                               ])) < 1e-3
        )

    def test_people_schedules_3(self):
        # Standard IHG
        sched = Schedule(
            "Percent1",
            "percent",
            np.array([0.1, .2, .3, .5] * int(CONFIG.number_of_time_steps_year / 4 + 1))[:-CONFIG.ts_per_hour + 1],
        )

        el1 = ElectricLoad(
            name='test_IHG',
            unit='W',
            nominal_value=100.,
            schedule=sched,
            fraction_to_zone=.9,
            fraction_radiant=0.45,
            fraction_convective=0.55,
        )

        conv, rad, lat, el = el1.get_loads()
        conv = conv[:4]
        rad = rad[:4]
        lat = lat[:4]
        assert (
                np.linalg.norm(conv - np.array([4.95,
                                                9.9,
                                                14.85,
                                                24.75,
                                                ])) < 1e-5 and
                np.linalg.norm(rad - np.array([4.05,
                                               8.1,
                                               12.15,
                                               20.25,
                                               ])) < 1e-5 and
                np.linalg.norm(lat - np.array([0,
                                               0,
                                               0,
                                               0,
                                               ])) < 1e-5
        )

        el2 = ElectricLoad(
            name='test_IHG',
            unit='W/m2',
            nominal_value=10.,
            schedule=sched,
            fraction_to_zone=.9,
            fraction_radiant=0.45,
            fraction_convective=0.55,
        )

        conv, rad, lat, el = el2.get_loads(area=2)
        conv = conv[:4]
        rad = rad[:4]
        lat = lat[:4]
        assert (
                np.linalg.norm(conv - np.array([0.99,
                                                1.98,
                                                2.97,
                                                4.95,
                                                ])) < 1e-5 and
                np.linalg.norm(rad - np.array([0.81,
                                               1.62,
                                               2.43,
                                               4.05,
                                               ])) < 1e-5 and
                np.linalg.norm(lat - np.array([0,
                                               0,
                                               0,
                                               0,
                                               ])) < 1e-5
        )

        el3 = ElectricLoad(
            name='test_IHG',
            unit='W/px',
            nominal_value=10.,
            schedule=sched,
            fraction_to_zone=.9,
            fraction_radiant=0.45,
            fraction_convective=0.55,
            number_of_people=3,
        )

        conv, rad, lat , el= el3.get_loads()
        conv = conv[:4]
        rad = rad[:4]
        lat = lat[:4]
        assert (
                np.linalg.norm(conv - np.array([1.485,
                                                2.97,
                                                4.455,
                                                7.425,

                                                ])) < 1e-5 and
                np.linalg.norm(rad - np.array([1.215,
                                               2.43,
                                               3.645,
                                               6.075,
                                               ])) < 1e-5 and
                np.linalg.norm(lat - np.array([0,
                                               0,
                                               0,
                                               0,
                                               ])) < 1e-5
        )

        l1 = Lights(
            name='test_IHG',
            unit='W/m2',
            nominal_value=10.,
            schedule=sched,
            fraction_to_zone=.9,
            fraction_radiant=0.45,
            fraction_convective=0.55,
        )

        conv, rad, lat, el = l1.get_loads(area=2)
        conv = conv[:4]
        rad = rad[:4]
        lat = lat[:4]
        assert (
                np.linalg.norm(conv - np.array([0.99,
                                                1.98,
                                                2.97,
                                                4.95,
                                                ])) < 1e-5 and
                np.linalg.norm(rad - np.array([0.81,
                                               1.62,
                                               2.43,
                                               4.05,
                                               ])) < 1e-5 and
                np.linalg.norm(lat - np.array([0,
                                               0,
                                               0,
                                               0,
                                               ])) < 1e-5
        )
