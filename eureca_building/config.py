"""
This module includes classes and functions to manage weather file
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import json
import os
import configparser
import logging
from datetime import datetime, timedelta
import eureca_building.logs

global DEFAULT_CONFIG_FILE
DEFAULT_CONFIG_FILE = os.path.join(".", "eureca_building", "default_config.ini")


def load_config(file: str = None):
    global CONFIG
    try:
        if file.endswith('ini'):
            CONFIG = Config()
            CONFIG.read(file)
        elif file.endswith('json'):
            CONFIG = Config.from_json(file)
        else:
            raise ValueError(f'Config file must be .json or .ini')
    except (FileNotFoundError, AttributeError):
        message = "Config file not found: Trying to load default config"
        print(message)
        logging.warning(message)
        CONFIG = Config()
        CONFIG.read(DEFAULT_CONFIG_FILE)
        message = "Default config loaded"
        print(message)

    globals().update(CONFIG)
    return CONFIG


# %% ---------------------------------------------------------------------------------------------------
class Config(configparser.ConfigParser):
    """
    This class is a container for config settings.

    Methods:
        to_json
        from_json
    """

    @property
    def ts_per_hour(self) -> int:
        return self._ts_per_hour

    @ts_per_hour.setter
    def ts_per_hour(self, value: int):
        try:
            value = int(value)
        except ValueError:
            raise TypeError(f"Config, time_step_per_hour is not an int: {value}")
        if 60 % value != 0:
            raise ValueError(f"Config, time_step_per_hour must be a divider of 60")
        self._ts_per_hour = value

    def read(self, file):
        super().read(file)
        # Generic config settings
        self.ts_per_hour = int(self['simulation settings']['time steps per hour'])
        self.start_date = datetime.strptime(self['simulation settings']['simulation reference year'] + ' ' + self['simulation settings']['start date'], "%Y %m-%d %H:%M")
        self.final_date = datetime.strptime(self['simulation settings']['simulation reference year'] + ' ' +self['simulation settings']['final date'], "%Y %m-%d %H:%M")
        self.heating_start_date = datetime.strptime(self['simulation settings']['simulation reference year'] + ' ' +self['simulation settings']['heating season start'], "%Y %m-%d %H:%M")
        self.heating_final_date = datetime.strptime(self['simulation settings']['simulation reference year'] + ' ' +self['simulation settings']['heating season end'], "%Y %m-%d %H:%M")
        self.cooling_start_date = datetime.strptime(self['simulation settings']['simulation reference year'] + ' ' +self['simulation settings']['cooling season start'], "%Y %m-%d %H:%M")
        self.cooling_final_date = datetime.strptime(self['simulation settings']['simulation reference year'] + ' ' +self['simulation settings']['cooling season end'], "%Y %m-%d %H:%M")
        self.time_step = int(3600 / self.ts_per_hour)  # s
        self.number_of_time_steps = int((self.final_date - self.start_date) / timedelta(
            minutes=self.time_step / 60)) + 1
        start_time_step = int(
            (self.start_date - datetime(self.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=self.time_step / 60))
        self.start_time_step = start_time_step
        self.final_time_step = start_time_step + self.number_of_time_steps

        # Heating/cooling season time steps
        self.heating_season_start_time_step = int(
            (self.heating_start_date - datetime(self.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=self.time_step / 60))
        self.heating_season_end_time_step = int(
            (self.heating_final_date - datetime(self.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=self.time_step / 60))

        self.cooling_season_start_time_step = int(
            (self.cooling_start_date - datetime(self.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=self.time_step / 60))
        self.cooling_season_end_time_step = int(
            (self.cooling_final_date - datetime(self.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=self.time_step / 60))

        self.number_of_time_steps_year = int(8760 * 60 / (self.time_step / 60))
        if self.ts_per_hour > 1:
            self.number_of_time_steps_year -= (self.ts_per_hour - 1)

        self.simulation_reference_year = int(self['simulation settings']['simulation reference year'])

        # Radiation
        self.azimuth_subdivisions = int(self['solar radiation settings']["azimuth subdivisions"])
        self.height_subdivisions = int(self['solar radiation settings']["height subdivisions"])

        self.do_solar_radiation_calculation = False if str.lower(self['solar radiation settings']["do solar radiation calculation"]) == "false" else True
        self.urban_shading_tolerances = [float(i) for i in
                                                self['solar radiation settings']["urban shading tolerances"].split(",")]

    @classmethod
    def from_json(cls, file_path):
        config_dict = cls()
        try:

            with open(file_path, "r") as json_data_file:
                config_dict.read_dict(json.load(json_data_file))
        except FileNotFoundError:
            raise FileNotFoundError(f"Config file {file_path} not found")
        # Generic config settings
        config_dict.ts_per_hour = int(config_dict['simulation settings']['time steps per hour'])
        config_dict.start_date = datetime.strptime(config_dict['simulation settings']['simulation reference year'] + ' ' +config_dict['simulation settings']['start date'], "%Y %m-%d %H:%M")
        config_dict.final_date = datetime.strptime(config_dict['simulation settings']['simulation reference year'] + ' ' +config_dict['simulation settings']['final date'], "%Y %m-%d %H:%M")
        config_dict.heating_start_date = datetime.strptime(config_dict['simulation settings']['simulation reference year'] + ' ' +config_dict['simulation settings']['heating season start'], "%Y %m-%d %H:%M")
        config_dict.heating_final_date = datetime.strptime(config_dict['simulation settings']['simulation reference year'] + ' ' +config_dict['simulation settings']['heating season end'], "%Y %m-%d %H:%M")
        config_dict.cooling_start_date = datetime.strptime(config_dict['simulation settings']['simulation reference year'] + ' ' +config_dict['simulation settings']['cooling season start'], "%Y %m-%d %H:%M")
        config_dict.cooling_final_date = datetime.strptime(config_dict['simulation settings']['simulation reference year'] + ' ' +config_dict['simulation settings']['cooling season end'], "%Y %m-%d %H:%M")
        config_dict.time_step = int(3600 / config_dict.ts_per_hour)  # s
        config_dict.number_of_time_steps = int((config_dict.final_date - config_dict.start_date) / timedelta(
            minutes=config_dict.time_step / 60)) + 1
        start_time_step = int(
            (config_dict.start_date - datetime(config_dict.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=config_dict.time_step / 60))
        config_dict.start_time_step = start_time_step
        config_dict.final_time_step = start_time_step + config_dict.number_of_time_steps

        # Heating/cooling season time steps
        config_dict.heating_season_start_time_step = int(
            (config_dict.heating_start_date - datetime(config_dict.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=config_dict.time_step / 60))
        config_dict.heating_season_end_time_step = int(
            (config_dict.heating_final_date - datetime(config_dict.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=config_dict.time_step / 60))

        config_dict.cooling_season_start_time_step = int(
            (config_dict.cooling_start_date - datetime(config_dict.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=config_dict.time_step / 60))
        config_dict.cooling_season_end_time_step = int(
            (config_dict.cooling_final_date - datetime(config_dict.start_date.year, 1, 1, 00, 00, 00)) / timedelta(
                minutes=config_dict.time_step / 60))

        config_dict.number_of_time_steps_year = int(8760 * 60 / (config_dict.time_step / 60))
        if config_dict.ts_per_hour > 1:
            config_dict.number_of_time_steps_year -= (config_dict.ts_per_hour - 1)

        config_dict.simulation_reference_year = int(config_dict['simulation settings']['simulation reference year'])

        # Radiation
        config_dict.azimuth_subdivisions = int(config_dict['solar radiation settings']["azimuth subdivisions"])
        config_dict.height_subdivisions = int(config_dict['solar radiation settings']["height subdivisions"])

        config_dict.do_solar_radiation_calculation = False if str.lower(
            config_dict['solar radiation settings']["do solar radiation calculation"]) == "false" else True
        config_dict.urban_shading_tolerances = [float(i) for i in
                                                config_dict['solar radiation settings']["urban shading tolerances"].split(",")]

        return config_dict

    def to_json(self, file_path):
        with open(file_path, "w") as outfile:
            json.dump(self, outfile, indent=4, )

# try:
#     CONFIG
# except NameError:
#     load_config()