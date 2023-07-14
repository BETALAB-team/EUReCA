"""This module includes a container class for schedule end-uses
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

'''IMPORTING MODULES'''

import os
import logging

import pandas as pd
import numpy as np

from eureca_building.config import CONFIG
from eureca_building.schedule import Schedule
from eureca_building.internal_load import People, Lights, ElectricLoad, InternalLoad
from eureca_building.setpoints import SetpointDualBand
from eureca_building.ventilation import Infiltration, MechanicalVentilation
from eureca_building.domestic_hot_water import DomesticHotWater

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions to create the schedule EndUses

def load_schedules(path):
    '''Schedule type loading in case you use daily schedules
    This function takes the path to the spreadsheet consisting of the schedules and loads; it
    Works for the daily schedule (see file ScheduleSemp.xlsx in eureca_ubem/Input/)
    
    Parameters
    ----------
    path : str
        Path to the spreadsheet to read file_schedule.xlsx

    Returns
    -------
    dict
        dictionary with EndUse_key/ eureca_ubem.end_uses.EndUse objects
    '''
    
    # Check input data type  

    # read archetype names from the excel sheet
    end_uses_sheet_dict = pd.read_excel(path,sheet_name=None,header=[0],index_col=[0], skiprows = 0)
    general_data = end_uses_sheet_dict["GeneralData"][['System Convective Fraction', 'AHU humidity control',
                                   'AHU sensible heat recovery', 'AHU latent heat recovery',
                                   'Outdoor Air Ratio','DomesticHotWater calculation']]

    holidays = [int(i) for i in list(end_uses_sheet_dict["GeneralData"]["Holidays from 0 to 364"].iloc[0].split(','))]

    end_uses_dict = {}

    for k_use, use_df in end_uses_sheet_dict.items():
        if k_use != "GeneralData":
            use_df = use_df.drop(use_df.index[0])
            end_uses_dict[k_use] = EndUse.load_daily_sched(k_use,use_df,general_data,holidays)

    return end_uses_dict

#%%--------------------------------------------------------------------------------------------------- 
#%% EndUse class

class EndUse:
    """This class builds up the loads and schedules for different type of buildings
    """

    def __init__(self,name):
        '''Creates the object with many dictionaries to store various loads and schedules
        
        Parameters
        ----------
        name : str
            name of the archetype

        '''    
        
        # Check input data type
        
        if not isinstance(name, str):
            raise TypeError(f'ERROR EndUse initialization, name must be a string: name {name}')
        
        # Inizialization
        
        self.name = name

        # Schedules in files are set using hourly timestep
        self.heat_gains = {
            'appliances':None,
            'lighting':None,
            'people':None,
        }
        self.domestic_hot_water = {
            'domestic_hot_water':None
        }
        self.infiltration = {
            'infiltration': None,
        }
        self.zone_system = {
            'temperature_setpoint': None,
            'humidity_setpoint': None,
            'convective_fraction':None,
        }
        self.air_handling_unit_system = {
            'ventilation_flow_rate': None,
            'ahu_availability': None,
            'ahu_humidity_control': None,
            'ahu_supply_temperature': None,
            'ahu_supply_humidity': None,
            "ahu_sensible_heat_recovery": None,
            "ahu_latent_heat_recovery": None,
            "outdoor_air_ratio": None,
        }
        '''
        Take care of units!!!
        '''

    @property
    def domestic_hot_water(self):
        return self._domestic_hot_water

    @domestic_hot_water.setter
    def domestic_hot_water(self, value):
        if not isinstance(value,dict):
            raise ValueError(f"EndUse object {self.name}, domestic_hot_water must be a dict: {value}")

        if (not isinstance(value["domestic_hot_water"], DomesticHotWater)) and value["domestic_hot_water"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting domestic hot water. Domestic hot water not a DomesticHotWater type")

        self._domestic_hot_water = value

    @property
    def heat_gains(self):
        return self._heat_gains

    @heat_gains.setter
    def heat_gains(self, value):
        if not isinstance(value,dict):
            raise ValueError(f"EndUse object {self.name}, heat gains must be a dict: {value}")

        if (not isinstance(value["appliances"], InternalLoad)) and value["appliances"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting internal loads. Appliances not a InternalLoad type")
        if (not isinstance(value["lighting"], InternalLoad)) and value["lighting"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting internal loads. Lights not a InternalLoad type")
        if (not isinstance(value["people"], People)) and value["people"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting internal loads. People not a People type")

        self._heat_gains = value

    @property
    def infiltration(self):
        return self._infiltration

    @infiltration.setter
    def infiltration(self, value):
        if not isinstance(value,dict):
            raise ValueError(f"EndUse object {self.name}, infiltration must be a dict: {value}")

        if (not isinstance(value["infiltration"], Infiltration)) and value["infiltration"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting infiltration. Infiltration not a Infiltration type")

        self._infiltration = value

    @property
    def zone_system(self):
        return self._zone_system

    @zone_system.setter
    def zone_system(self, value):
        if not isinstance(value,dict):
            raise ValueError(f"EndUse object {self.name}, zone_system must be a dict: {value}")

        if (not isinstance(value["temperature_setpoint"], Infiltration)) and value["temperature_setpoint"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting zone_system. temperature_setpoint not a SetpointDualBand type")
        if (not isinstance(value["humidity_setpoint"], Infiltration)) and value["humidity_setpoint"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting zone_system. humidity_setpoint not a SetpointDualBand type")
        if (not isinstance(value["convective_fraction"], float)) and value["convective_fraction"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting zone_system. humidity_setpoint not a SetpointDualBand type")

        self._zone_system = value


    @property
    def air_handling_unit_system(self):
        return self._air_handling_unit_system

    @air_handling_unit_system.setter
    def air_handling_unit_system(self, value):
        if not isinstance(value,dict):
            raise ValueError(f"EndUse object {self.name}, air_handling_unit_system must be a dict: {value}")

        if (not isinstance(value["ventilation_flow_rate"], MechanicalVentilation)) and value["ventilation_flow_rate"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting air_handling_unit_system. ventilation_flow_rate not a MechanicalVentilation type")
        if (not isinstance(value["ahu_availability"], Schedule)) and value["ahu_availability"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting air_handling_unit_system. ahu_availability not a Schedule type")

        if (not isinstance(value["ahu_supply_temperature"], Schedule)) and value["ahu_supply_temperature"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting air_handling_unit_system. ahu_supply_temperature not a Schedule type")
        if (not isinstance(value["ahu_supply_humidity"], Schedule)) and value["ahu_supply_humidity"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting air_handling_unit_system. ahu_supply_humidity not a Schedule type")

        if (not isinstance(value["ahu_humidity_control"], bool)) and value["ahu_humidity_control"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting zone_system. ahu_humidity_control not a boolean")
        if (not isinstance(value["ahu_sensible_heat_recovery"], float)) and value["ahu_sensible_heat_recovery"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting zone_system. ahu_sensible_heat_recovery not a float")
        if (not isinstance(value["ahu_latent_heat_recovery"], float)) and value["ahu_latent_heat_recovery"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting zone_system. ahu_latent_heat_recovery not a float")
        if (not isinstance(value["outdoor_air_ratio"], float)) and value["outdoor_air_ratio"] != None:
            raise TypeError(f"EndUse object {self.name}: error in setting zone_system. outdoor_air_ratio not a float")

        self._air_handling_unit_system = value

    @classmethod
    def load_daily_sched(cls, name, daily_df_from_excel, scalar_df_from_excel, holidays):
        '''Class method to create the EndUse object from the spreadsheet page

        Parameters
        ----------
        name : str
            name
        daily_df_from_excel : pandas.DataFrame
            DataFrame containing schedules from the end_use
        scalar_df_from_excel : pandas.DataFrame
            This series includes some additional data about the archetype (from the GeneralInfo page in the spreadsheet)
                                                                            (Sensible and
                                                                           Latent AHU recovery,
                                                                           Convective fraction of internal gains)

        Returns
        -------
        eureca_ubem.end_uses.EndUse
        '''

        # Each schedule is set to a different attribute of the class

        schedules_list = [
                'appliances',
                'lighting',
                'people',
                'heatingTSP',
                'coolingTSP',
                'HeatingRHSP',
                'CoolingRHSP',
                'ventFlowRate',
                'infFlowRate',
                'ahu_supply_temperature_heating',
                'ahu_supply_specific_humidity_heating',
                'ahu_supply_temperature_cooling',
                'ahu_supply_specific_humidity_cooling',
            ]

        sched_df = pd.DataFrame(columns = schedules_list)
        sched_df['appliances'] = daily_df_from_excel['Appliances']
        sched_df['lighting'] = daily_df_from_excel['Lighting']
        sched_df['people'] = daily_df_from_excel['Occupancy (Total)']
        sched_df['heatingTSP'] = daily_df_from_excel['HeatSP']
        sched_df['coolingTSP'] = daily_df_from_excel['CoolSP']
        sched_df['HeatingRHSP'] = daily_df_from_excel['HumSP'] / 100
        sched_df['CoolingRHSP'] = daily_df_from_excel['DehumSP'] / 100
        sched_df['ventFlowRate'] = daily_df_from_excel['Ventilation FlowRate']
        sched_df['infFlowRate'] = daily_df_from_excel['Infiltration FlowRate']
        sched_df['ahu_supply_temperature_heating'] = daily_df_from_excel['Ventilation Supply Temperature heating mode']
        sched_df['ahu_supply_specific_humidity_heating'] = daily_df_from_excel['Ventilation Supply specific humidity heating mode']
        sched_df['ahu_supply_temperature_cooling'] = daily_df_from_excel['Ventilation Supply Temperature cooling mode']
        sched_df['ahu_supply_specific_humidity_cooling'] = daily_df_from_excel['Ventilation Supply specific humidity cooling mode']
        sched_df['dhw'] = daily_df_from_excel['Domestic Hot Water FlowRate']

        sched_df = cls.rescale_df(CONFIG.ts_per_hour, sched_df)

        scalar_data = {}
        try:
            scalar_data['conFrac'] = float(scalar_df_from_excel.loc[name]['System Convective Fraction'])
            scalar_data['AHUHUM'] = bool(scalar_df_from_excel.loc[name]['AHU humidity control'])
            scalar_data['sensRec'] = float(scalar_df_from_excel.loc[name]['AHU sensible heat recovery'])
            scalar_data['latRec'] = float(scalar_df_from_excel.loc[name]['AHU latent heat recovery'])
            scalar_data['outdoorAirRatio'] = float(scalar_df_from_excel.loc[name]['Outdoor Air Ratio'])
            scalar_data['DomesticHotWater calculation'] = str(scalar_df_from_excel.loc[name]['DomesticHotWater calculation'])
        except KeyError:
            raise KeyError(
                f"ERROR Loading end use {name}. GeneralData does not have the correct columns names: ConvFrac, AHUHum, SensRec, LatRec, OutAirRatio")
        except ValueError:
            raise ValueError(f"""ERROR 
                             Loading end use {name}. GeneralData
                             I'm not able to parse the General data. 
                                 ConvFrac should be a float {scalar_df_from_excel.loc[name]['ConvFrac']}
                                 AHUHum should be a boolean {scalar_df_from_excel.loc[name]['AHUHum']}
                                 SensRec should be a float {scalar_df_from_excel.loc[name]['SensRec']}
                                 LatRec  should be a float {scalar_df_from_excel.loc[name]['LatRec']}
                                 OutAirRatio   should be a float {scalar_df_from_excel.loc[name]['OutAirRatio']}
                                 DomesticHotWater calculation   should be a str {scalar_df_from_excel.loc[name]['DomesticHotWater calculation']}
                             """)

        # Check the quality of input data
        if not 0. <= scalar_data['conFrac'] <= 1.:
            logging.warning(f"WARNING Loading end use {name}. Convective fraction of the heat gain outside boundary condition [0-1]: ConvFrac {scalar_data['conFrac']}")
        if not 0. <= scalar_data['sensRec'] <= 1.:
            logging.warning(f"WARNING Loading end use {name}. Sensible recovery of the AHU outside boundary condition [0-1]: sensRec {scalar_data['sensRec']}")
        if not 0. <= scalar_data['latRec'] <= 1.:
            logging.warning(f"WARNING Loading end use {name}. Latent recovery of the AHU outside boundary condition [0-1]: sensRec {scalar_data['latRec']}")
        if not 0. <= scalar_data['outdoorAirRatio'] <= 1.:
            logging.warning(f"WARNING Loading end use {name}. Outdoor air ratio of the AHU outside boundary condition [0-1]: outdoorAirRatio {scalar_data['outdoorAirRatio']}")

        # Creating the end use object

        people_sched = Schedule.from_daily_schedule(
            f"People sched {name}",
            "dimensionless",
            schedule_week_day = sched_df["people"].iloc[:24 * CONFIG.ts_per_hour].values,
            schedule_saturday = sched_df["people"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour].values,
            schedule_sunday = sched_df["people"].iloc[48 * CONFIG.ts_per_hour:24*3 * CONFIG.ts_per_hour].values,
            schedule_holiday = sched_df["people"].iloc[48 * CONFIG.ts_per_hour:24*3 * CONFIG.ts_per_hour].values,
            starting_day = 0,
            holidays = holidays
        )

        people = People(
            # The schedule is already in W/m2
            name=f'People load {name}',
            unit='W/m2',
            nominal_value=1.,
            schedule=people_sched,
            fraction_latent=0.43,
            fraction_radiant=0.3,
            fraction_convective=0.7,
            metabolic_rate=1,
        )

        app_sched = Schedule.from_daily_schedule(
            f"Appliances sched {name}",
            "dimensionless",
            schedule_week_day=sched_df["appliances"].iloc[:24 * CONFIG.ts_per_hour ],
            schedule_saturday=sched_df["appliances"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["appliances"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["appliances"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        app = ElectricLoad(
            name=f'Appliances load {name}',
            unit='W/m2',
            nominal_value=1.,
            schedule=app_sched,
            fraction_radiant=0.5,
            fraction_convective=0.5,
        )

        lights_sched = Schedule.from_daily_schedule(
            f"Lights sched {name}",
            "dimensionless",
            schedule_week_day=sched_df["lighting"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["lighting"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["lighting"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["lighting"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        lights = Lights(
            name=f"Lights load {name}",
            unit='W/m2',
            nominal_value=1.,
            schedule=lights_sched,
            fraction_radiant=0.7,
            fraction_convective=0.3,
        )

        heat_sp_sched = Schedule.from_daily_schedule(
            f"Heating setpoint {name}",
            "temperature",
            schedule_week_day=sched_df["heatingTSP"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["heatingTSP"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["heatingTSP"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["heatingTSP"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        cool_sp_sched = Schedule.from_daily_schedule(
            f"Cooling setpoint {name}",
            "temperature",
            schedule_week_day=sched_df["coolingTSP"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["coolingTSP"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["coolingTSP"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["coolingTSP"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        hum_sp_sched = Schedule.from_daily_schedule(
            f"Humidifying setpoint {name}",
            "dimensionless",
            schedule_week_day=sched_df["HeatingRHSP"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["HeatingRHSP"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["HeatingRHSP"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["HeatingRHSP"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        dehum_sp_sched = Schedule.from_daily_schedule(
            f"Dehumidifying setpoint {name}",
            "dimensionless",
            schedule_week_day=sched_df["CoolingRHSP"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["CoolingRHSP"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["CoolingRHSP"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["CoolingRHSP"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        heat_sp_sched.schedule[CONFIG.heating_season_end_time_step: CONFIG.heating_season_start_time_step] = -100
        hum_sp_sched.schedule[CONFIG.heating_season_end_time_step: CONFIG.heating_season_start_time_step] = -100

        cool_sp_sched.schedule[:CONFIG.cooling_season_start_time_step] = 150
        cool_sp_sched.schedule[CONFIG.cooling_season_end_time_step:] = 150

        dehum_sp_sched.schedule[:CONFIG.cooling_season_start_time_step] = 1.5
        dehum_sp_sched.schedule[CONFIG.cooling_season_end_time_step:] = 1.5

        temp_sp = SetpointDualBand(
            "t_sp",
            "temperature",
            schedule_lower=heat_sp_sched,
            schedule_upper=cool_sp_sched,
        )
        heat_sp = SetpointDualBand(
            "h_sp",
            "relative_humidity",
            schedule_lower=hum_sp_sched,
            schedule_upper=dehum_sp_sched,
        )

        infiltration = Schedule.from_daily_schedule(
            f"Infiltration sched {name}",
            "dimensionless",
            schedule_week_day=sched_df["infFlowRate"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["infFlowRate"].iloc[24 * CONFIG.ts_per_hour :48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["infFlowRate"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["infFlowRate"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        inf_obj = Infiltration(
            name='inf_obj',
            unit='Vol/h',
            nominal_value=1.,
            schedule=infiltration,
        )

        vent_sched = Schedule.from_daily_schedule(
            f"Ventilation mass flow rate sched {name}",
            "dimensionless",
            schedule_week_day=sched_df["ventFlowRate"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["ventFlowRate"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["ventFlowRate"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["ventFlowRate"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        vent_obj = MechanicalVentilation(
            name='vent_obj',
            unit='m3/(m2 s)',
            nominal_value=1,
            schedule=vent_sched,
        )

        supply_t_heating = Schedule.from_daily_schedule(
            f"Ventilation temperature supply sched {name}",
            "dimensionless",
            schedule_week_day=sched_df["ahu_supply_temperature_heating"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["ahu_supply_temperature_heating"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["ahu_supply_temperature_heating"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["ahu_supply_temperature_heating"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        supply_h_heating = Schedule.from_daily_schedule(
            f"Ventilation humidity supply sched {name}",
            "dimensionless",
            schedule_week_day=sched_df["ahu_supply_specific_humidity_heating"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["ahu_supply_specific_humidity_heating"].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["ahu_supply_specific_humidity_heating"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["ahu_supply_specific_humidity_heating"].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        supply_t_cooling = Schedule.from_daily_schedule(
            f"Ventilation temperature supply sched {name}",
            "dimensionless",
            schedule_week_day=sched_df["ahu_supply_temperature_cooling"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["ahu_supply_temperature_cooling"].iloc[
                              24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["ahu_supply_temperature_cooling"].iloc[
                            48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["ahu_supply_temperature_cooling"].iloc[
                             48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        supply_h_cooling = Schedule.from_daily_schedule(
            f"Ventilation humidity supply sched {name}",
            "dimensionless",
            schedule_week_day=sched_df["ahu_supply_specific_humidity_cooling"].iloc[:24 * CONFIG.ts_per_hour],
            schedule_saturday=sched_df["ahu_supply_specific_humidity_cooling"].iloc[
                              24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour],
            schedule_sunday=sched_df["ahu_supply_specific_humidity_cooling"].iloc[
                            48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            schedule_holiday=sched_df["ahu_supply_specific_humidity_cooling"].iloc[
                             48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour],
            starting_day=0,
            holidays=holidays
        )

        # Change heating cooling season
        supply_t_heating.schedule[CONFIG.cooling_season_start_time_step: CONFIG.cooling_season_end_time_step] = supply_t_cooling.schedule[CONFIG.cooling_season_start_time_step: CONFIG.cooling_season_end_time_step]
        supply_h_heating.schedule[CONFIG.cooling_season_start_time_step: CONFIG.cooling_season_end_time_step] = supply_h_cooling.schedule[CONFIG.cooling_season_start_time_step: CONFIG.cooling_season_end_time_step]

        ahu_availability_sched = Schedule.from_constant_value(
            "ahu_availability_sched",
            "availability",
            1,
        )

        ahu_availability_sched.schedule[CONFIG.cooling_season_start_time_step: CONFIG.cooling_season_end_time_step] = -1
        ahu_availability_sched.schedule[CONFIG.heating_season_end_time_step: CONFIG.cooling_season_start_time_step] = 0
        ahu_availability_sched.schedule[CONFIG.cooling_season_end_time_step: CONFIG.heating_season_start_time_step] = 0

        domestic_hot_water_sched = Schedule.from_daily_schedule(
            f"DHW sched {name}",
            "mass_flow_rate",
            schedule_week_day=sched_df['dhw'].iloc[:24 * CONFIG.ts_per_hour].values,
            schedule_saturday=sched_df['dhw'].iloc[24 * CONFIG.ts_per_hour:48 * CONFIG.ts_per_hour].values,
            schedule_sunday=sched_df['dhw'].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour].values,
            schedule_holiday=sched_df['dhw'].iloc[48 * CONFIG.ts_per_hour:24 * 3 * CONFIG.ts_per_hour].values,
            starting_day=0,
            holidays=holidays
        )

        dhw = DomesticHotWater(
            f'DomesticHotWater {name}',
            calculation_method=scalar_data['DomesticHotWater calculation'],
            unit="L/(m2 h)",
            schedule=domestic_hot_water_sched,
        )

        end_use_obj = cls(name)

        end_use_obj.heat_gains['appliances'] = app # InternalLoad type
        end_use_obj.heat_gains['lighting'] = lights # InternalLoad type
        end_use_obj.heat_gains['people'] = people # People Type

        end_use_obj.domestic_hot_water['domestic_hot_water'] = dhw

        end_use_obj.infiltration['infiltration'] = inf_obj # Infiltration Type

        end_use_obj.zone_system['temperature_setpoint'] = temp_sp # SetpointDualBand type
        end_use_obj.zone_system['humidity_setpoint'] = heat_sp # SetpointDualBand type
        end_use_obj.zone_system['convective_fraction'] = scalar_data['conFrac']

        end_use_obj.air_handling_unit_system['ventilation_flow_rate'] = vent_obj # MechanicalVentilation
        end_use_obj.air_handling_unit_system['ahu_availability'] = ahu_availability_sched # Schedule
        end_use_obj.air_handling_unit_system['ahu_humidity_control'] = scalar_data['AHUHUM'] # boolean
        end_use_obj.air_handling_unit_system['ahu_supply_temperature'] = supply_t_heating # Schedule
        end_use_obj.air_handling_unit_system['ahu_supply_humidity'] = supply_h_heating # Schedule
        end_use_obj.air_handling_unit_system['ahu_sensible_heat_recovery'] = scalar_data['sensRec'] # float
        end_use_obj.air_handling_unit_system['ahu_latent_heat_recovery'] = scalar_data['latRec'] # float
        end_use_obj.air_handling_unit_system['outdoor_air_ratio'] = scalar_data['outdoorAirRatio'] # float

        return end_use_obj


    # @classmethod
    # def load_daily_sched(self,daily_df_from_excel, scale_df_from_excel):
    #
    #     '''
    #     Used for ScheduleComp.xlsx Excel file
    #
    #     Parameters
    #         ----------
    #         arch : pandas dataframe
    #             This must include all the schedules' keys
    #
    #         sched : pandas dataframe
    #             This dataframe includes all the yearly schedules
    #
    #     Returns
    #     -------
    #     None.
    #
    #     '''
    #
    #     # Each schedule is set to a different attribute of the class
    #
    #     try:
    #         self.sched_df['appliances'] = sched[arch['Appliances']]
    #         self.sched_df['lighting'] = sched[arch['Lighting']]
    #         self.sched_df['people'] = sched[arch['People (Sensible)']]
    #         self.sched_df['vapour'] = sched[arch['Vapour']]
    #         self.sched_df['heatingTSP'] = sched[arch['HeatTSP']]
    #         self.sched_df['coolingTSP'] = sched[arch['CoolTSP']]
    #         self.sched_df['HeatingRHSP'] = sched[arch['HeatRHSP']]
    #         self.sched_df['CoolingRHSP'] = sched[arch['CoolRHSP']]
    #         self.sched_df['ventFlowRate'] = sched[arch['VentFlowRate']]
    #         self.sched_df['infFlowRate'] = sched[arch['InfFlowRate']]
    #         self.sched_df['plantOnOffSens'] = sched[arch['PlantONOFFSens']]
    #         self.sched_df['plantOnOffLat'] = sched[arch['PlantONOFFLat']]
    #         self.sched_df['AHUOnOff'] = sched[arch['AHUONOFF']]
    #         self.sched_df['AHUHUM'] = sched[arch['AHUHUM']]
    #         self.sched_df['AHUTSupp'] = sched[arch['AHUTSupp']]
    #         self.sched_df['AHUxSupp'] = sched[arch['AHUxSupp']]
    #
    #         self.scalar_data['conFrac'] = float(arch['ConvFrac'])
    #         self.scalar_data['AHUHUM'] = bool(arch['AHUHum'])
    #         self.scalar_data['sensRec'] = float(arch['SensRec'])
    #         self.scalar_data['latRec'] = float(arch['LatRec'])
    #         self.scalar_data['outdoorAirRatio'] = float(arch['OutAirRatio'])
    #     except KeyError:
    #         raise KeyError(f'ERROR Archetype object {self.name}: can not find all schedules')

    @staticmethod
    def rescale_df(ts, sched_df):
        '''Static method to rescale the schedule dataframe from the hour to the simulation time step
        
        Parameters
        ----------
        ts : int
            Number of time steps per hour
        sched_df : pandas.DataFrame
            DataFrame to rescale
                
        Returns
        -------
        pandas.DataFrame
            The rescaled df
        '''   
        
        # Check input data type
        
        if not isinstance(ts, int):
            raise TypeError(f'ERROR input ts is not an integer: ts {ts}')         
    
        # Rescale 
        
        m = str(60/ts) + 'min'
        sched_df = pd.concat([sched_df,sched_df.iloc[-1]])
        time = pd.date_range('2020-01-01', periods=len(sched_df.index), freq='1h')
        sched_df.set_index(time, inplace=True)
        sched_df = sched_df.resample(m).ffill()                               # Steps interpolation Resampling
        return sched_df
        #Boundary = Boundary_0.resample(str(ts)+'S').interpolate(method='linear')        # Linear interpolation Resampling
        # There are several upsample methods: pad(), bfill(), mean(), interpolate(), apply custum function...
        # For info look:
        #   https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.resample.html
        #   https://machinelearningmastery.com/resample-interpolate-time-series-data-python/
        
    # def create_np(self):
    #
    #     '''
    #     Creates the np.array attributes from the self.sched_df dataframe
    #
    #     Parameters
    #         ----------
    #
    #     Returns
    #     -------
    #     None.
    #     '''
    #
    #     self.appliances = self.sched_df['appliances'].to_numpy(dtype = np.float_)
    #     self.lighting = self.sched_df['lighting'].to_numpy(dtype = np.float_)
    #     self.people = self.sched_df['people'].to_numpy(dtype = np.float_)
    #     self.vapour = self.sched_df['vapour'].to_numpy(dtype = np.float_)
    #     self.heatingTSP = self.sched_df['heatingTSP'].to_numpy(dtype = np.float_)
    #     self.coolingTSP = self.sched_df['coolingTSP'].to_numpy(dtype = np.float_)
    #     self.HeatingRHSP = self.sched_df['HeatingRHSP'].to_numpy(dtype = np.float_)
    #     self.CoolingRHSP = self.sched_df['CoolingRHSP'].to_numpy(dtype = np.float_)
    #     self.ventFlowRate = self.sched_df['ventFlowRate'].to_numpy(dtype = np.float_)
    #     self.infFlowRate = self.sched_df['infFlowRate'].to_numpy(dtype = np.float_)
    #     self.plantOnOffSens = self.sched_df['plantOnOffSens'].to_numpy(dtype = np.int_)
    #     self.plantOnOffLat = self.sched_df['plantOnOffLat'].to_numpy(dtype = np.int_)
    #     self.AHUOnOff = self.sched_df['AHUOnOff'].to_numpy(dtype = np.int_)
    #     #self.AHUHUM = self.sched_df['AHUHUM'].to_numpy(dtype = np.bool_)
    #     self.AHUTSupp = self.sched_df['AHUTSupp'].to_numpy(dtype = np.float_)
    #     self.AHUxSupp = self.sched_df['AHUxSupp'].to_numpy(dtype = np.float_)
    #     #self.conFrac = self.sched_df['conFrac'].to_numpy(dtype = np.float_)
    #     #self.AHUHumidistat = self.sched_df['AHUHum']
    #     #self.sensRec = self.sched_df['sensRec'].to_numpy(dtype = np.float_)
    #     #self.latRec = self.sched_df['latRec'].to_numpy(dtype = np.float_)
    #     #self.outdoorAirRatio = self.sched_df['outdoorAirRatio'].to_numpy(dtype = np.float_)

