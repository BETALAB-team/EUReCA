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
from tqdm import tqdm
import pandas as pd
import numpy as np
import datetime as dt

from eureca_building.config import CONFIG
from eureca_building.schedule import Schedule
from eureca_building.internal_load import People, Lights, ElectricLoad, InternalLoad
from eureca_building.setpoints import SetpointDualBand
from eureca_building.ventilation import Infiltration, MechanicalVentilation
from eureca_building.domestic_hot_water import DomesticHotWater

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions to create the schedule EndUses

def load_schedules(City, path):
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
    Case = "Classic EUReCA"
    if path.lower().endswith(".csv"):
        end_uses_csv = pd.read_csv(path, sep=";")
        end_uses_csv_tag = end_uses_csv["tag"].iloc[0]
        if end_uses_csv_tag == "Markov_UU":
            Case = "Markov_UU"
    # read archetype names from the excel sheet
    if Case == "Classic EUReCA":
        
        end_uses_sheet_dict = pd.read_excel(path,sheet_name=None,header=[0],index_col=[0], skiprows = 0)
        general_data = end_uses_sheet_dict["GeneralData"][['System Convective Fraction', 'AHU humidity control',
                                       'AHU sensible heat recovery', 'AHU latent heat recovery',
                                       'Outdoor Air Ratio','DomesticHotWater calculation','Appliances calculation']]
    
        holidays = [int(i) for i in list(end_uses_sheet_dict["GeneralData"]["Holidays from 0 to 364"].iloc[0].split(','))]
    
        end_uses_dict = {}
        end_uses_dict["tag"] = "Classic"
        for k_use, use_df in end_uses_sheet_dict.items():
            if k_use != "GeneralData":
                use_df = use_df.drop(use_df.index[0])
                end_uses_dict[k_use] = EndUse.load_daily_sched(k_use,use_df,general_data,holidays)
                
            
        return end_uses_dict
    if Case == "Markov_UU":
        weekdays_matrix, holidays_matrix, holidays_list, setpoints, settings = load_markov_csv(path)
        end_uses_dict = {}
        end_uses_dict["tag"] = "Markov UU"
        for building_id, building_object in tqdm(City.buildings_objects.items(), desc="Occupancy Simulation"):
            n_people = City.buildings_info[building_id]["Number of occupants"]
            
            individual_schedule = simulate_occupancy(
                wd_matrix = weekdays_matrix,
                we_matrix =  holidays_matrix,
                holidays = holidays_list,
                n_people = n_people,
                random_state = 42)
            individual_enduse_object = load_daily_sched_individual(str(building_id), individual_schedule, holidays_list,  settings, setpoints)
            end_uses_dict[building_id]=individual_enduse_object
        return end_uses_dict


#%%--------------------------------------------------------------------------------------------------- 
#%% Markov 
def load_markov_csv (markovian_csv_file, n_hours = 24, sep = ";"):
    """
    Read a CSV of hourly transition probabilities and build 3D transition matrices.

    Expected columns:
        'hour'        : int in [0, 23]
        'fromstate'   : state id (e.g. 1, 2, 3, ...)
        'wd_to_stateX': weekday prob from 'fromstate' to state X
        'we_to_stateX': weekend/off-day prob from 'fromstate' to state X

    Returns
    -------
    wd_matrix : np.ndarray
        Shape (n_hours, n_states, n_states)
    we_matrix : np.ndarray
        Shape (n_hours, n_states, n_states)
    """
    df = pd.read_csv(markovian_csv_file, sep = sep)

    wd_cols = [c for c in df.columns if c.startswith('wd_to_state')]
    we_cols = [c for c in df.columns if c.startswith('we_to_state')]
    if len(wd_cols) != len(we_cols):
        raise ValueError("Number of weekday and weekend 'to_state' columns does not match.")
    n_states = len(wd_cols)

    from_states = sorted(df['fromstate'].unique())
    state_index = {s: i for i, s in enumerate(from_states)}

    wd_matrix = np.zeros((n_hours, n_states, n_states), dtype=float)
    we_matrix = np.zeros((n_hours, n_states, n_states), dtype=float)

    for _, row in df.iterrows():
        h = int(row['hour'])
        fs = state_index[row['fromstate']]
        if not (0 <= h < n_hours):
            raise ValueError(f"Hour {h} out of range [0, {n_hours-1}].")

        wd_probs = np.array([row[c] for c in wd_cols], dtype=float)
        we_probs = np.array([row[c] for c in we_cols], dtype=float)

        if wd_probs.size != n_states or we_probs.size != n_states:
            raise ValueError("Mismatch between number of states and columns.")

        wd_matrix[h, fs, :] = wd_probs
        we_matrix[h, fs, :] = we_probs
        
    if 'holidays' in df.columns:
        holidays = sorted(df['holidays'].dropna().astype(int).unique().tolist())
    else:
        holidays = []
        
    setpoint_cols = [c for c in df.columns if '_sp_' in c]

    setpoint_key_map = {
        'cooling_sp_celsius': 'cooling_setpoint_temp_c',
        'heating_sp_celsius': 'heating_setpoint_temp_c',
        'hum_sp_pc': 'humidity_rel_pct',
        'dehum_sp_pc': 'dehumidification_rel_pct',
    }    
    setpoints = {}
    for col in setpoint_cols:
        series = df[col].dropna()
        if series.empty:
            continue
        key = setpoint_key_map.get(col, col)
        setpoints[key] = float(series.iloc[0])
    settings_key_map = {
        've_onoff': 'ventilation_onoff',
        've_supp_t_h_celsius': 'ventilation_supply_temp_heating_c',
        've_supp_t_c_celsius': 'ventilation_supply_temp_cooling_c',
        've_supp_sh_h_kgkg': 'ventilation_supply_humidity_heating_kgkg',
        've_supp_sh_c_kgkg': 'ventilation_supply_humidity_cooling_kgkg',
        'raditor_in_temp_celsius': 'radiator_inlet_temp_c',
        'radiator_out_temp_celsius': 'radiator_outlet_temp_c',
        'elec_pow_fac': 'electric_power_factor',
        'substation_voltage_a_volt': 'substation_voltage_volt',
        'susbtation_voltage_phase_angle_deg': 'substation_voltage_phase_angle_deg',
        'System Convective Fraction':'System Convective Fraction',
        'AHU humidity control':'AHU humidity control',
        'AHU sensible heat recovery':'AHU sensible heat recovery',
        'AHU latent heat recovery':'AHU latent heat recovery',
        'Outdoor Air Ratio':'Outdoor Air Ratio',
        'DomesticHotWater calculation':'DomesticHotWater calculation',
        'Appliances calculation':'Appliances calculation'        
    }

    settings_cols = [c for c in df.columns if c in settings_key_map]
    settings = {}
    for col in settings_cols:
        s = df[col].dropna()
        if not s.empty:
            key = settings_key_map[col]
            if not key in ["DomesticHotWater calculation","Appliances calculation", "AHU humidity control"]:
               settings[key] = float(s.iloc[0])
            if key in ["DomesticHotWater calculation","Appliances calculation"]:
                settings[key] = str(s.iloc[0])
            if key == "AHU humidity control":
                settings[key] = bool(s.iloc[0])
    return wd_matrix, we_matrix, holidays, setpoints, settings


def simulate_occupancy(
    wd_matrix: np.ndarray,
    we_matrix: np.ndarray,
    holidays: list,
    n_people: int,
    random_state: int,
    CONFIG = CONFIG
) -> pd.DataFrame:

    """
    Simulate occupancy with a Markov chain, using a warm-up period before start_date.

    Returns a DataFrame only for [start_date, end_date)
    """

    rng = np.random.default_rng(random_state)

    n_hours, n_states_from, n_states_to = wd_matrix.shape
    assert n_hours == 24, "Transition matrices must have 24 hourly slices."
    assert n_states_from == n_states_to, "Transition matrices must be square."
    n_states = n_states_from
    warmup_days = 5
    holidays_set = set(int(d) for d in holidays)
    start_date, end_date, ts_per_hour = CONFIG.start_date, CONFIG.final_date, CONFIG.ts_per_hour
    sim_start = start_date - dt.timedelta(days=warmup_days)
    total_hours = int((end_date - sim_start).total_seconds() // 3600)
    n_steps = total_hours * ts_per_hour
    step_delta = dt.timedelta(minutes=60 / ts_per_hour)
    current_states = np.zeros(n_people, dtype=int)

    results = []
    timestamps = []

    for step in range(n_steps):
        current_dt = sim_start + step * step_delta

        hour = current_dt.hour
        day_of_year = current_dt.timetuple().tm_yday
        weekday = current_dt.weekday()  # Monday=0, Sunday=6

        is_off_day = (weekday >= 5) or (day_of_year in holidays_set)
        T = we_matrix[hour] if is_off_day else wd_matrix[hour]

        new_states = np.empty_like(current_states)
        for i, s in enumerate(current_states):
            probs = T[s, :]
            if probs.sum() <= 0:
                probs = np.zeros_like(probs)
                probs[s] = 1.0
            else:
                probs = probs / probs.sum()
            new_states[i] = rng.choice(n_states, p=probs)

        current_states = new_states
        counts = np.bincount(current_states, minlength=n_states)

        results.append(counts)
        timestamps.append(current_dt)

    # full DF including warmup
    data = np.vstack(results)
    columns = [f"state_{i}" for i in range(n_states)]
    df_full = pd.DataFrame(data=data, index=pd.DatetimeIndex(timestamps), columns=columns)

    # slice to requested window [start_date, end_date)
    mask = (df_full.index >= start_date) & (df_full.index < end_date)
    df = df_full.loc[mask].copy()

    return df

def load_daily_sched_individual(name, ind_occupancy_schedule, holidays,  schedule_settings, schedule_setpoints):
    """

    """
    occ = ind_occupancy_schedule.copy()
    if not isinstance(occ.index, pd.DatetimeIndex):
            raise TypeError("ind_occupancy_schedule must have a DatetimeIndex")

    STATE_GAINS = {
        0: {"appliances": 50.0, "hot_water": 0.0,  "body_heat": 0.0},
        1: {"appliances": 50.0, "hot_water": 0.0,  "body_heat": 80.0},
        2: {"appliances": 100.0, "hot_water": 13,  "body_heat": 110.0},
    }

    state_cols = sorted([c for c in occ.columns if c.startswith("state_")])
    if not state_cols:
        raise ValueError("ind_occupancy_schedule must have columns like 'state_1', 'state_2', ...")

    n_steps = len(occ)
    appliances = np.zeros(n_steps)
    hot_water = np.zeros(n_steps)
    body_heat = np.zeros(n_steps)
    awake_people = np.zeros(n_steps)
    for col in state_cols:
        state_id = int(col.split("_")[1])
        if state_id not in STATE_GAINS:
            raise ValueError(f"No gains defined for state {state_id}")

        n_people_state = occ[col].to_numpy(dtype=float)
        gains = STATE_GAINS[state_id]

        appliances += n_people_state * gains["appliances"]
        hot_water += n_people_state * gains["hot_water"]
        body_heat += n_people_state * gains["body_heat"]
        if state_id == 2:
            awake_people += n_people_state 
    lighting = np.zeros(n_steps)
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
        'dhw',
    ]
    sched_df = pd.DataFrame(index=occ.index, columns=schedules_list, dtype=float)
    sched_df["appliances"] = appliances
    sched_df["lighting"]   = lighting
    sched_df["people"]     = body_heat
    sched_df["dhw"]        = awake_people
    sched_df["heatingTSP"] = schedule_setpoints.get("heating_temp_c", 20.0)
    sched_df["coolingTSP"] = schedule_setpoints.get("cooling_temp_c", 26.0)
    sched_df["HeatingRHSP"] = schedule_setpoints.get("humidity_rel_pct", 30)/100
    sched_df["CoolingRHSP"] = schedule_setpoints.get("dehumidification_rel_pct", 60)/100
    sched_df["ventFlowRate"] = schedule_settings.get("ventilation_onoff", 1.0)
    sched_df["infFlowRate"] = schedule_settings.get("infiltration_ach", 1.0)
    sched_df["ahu_supply_temperature_heating"] = schedule_settings.get(
        "ventilation_supply_temp_heating_c", 20.0
    )
    sched_df["ahu_supply_specific_humidity_heating"] = schedule_settings.get(
        "ventilation_supply_humidity_heating_kgkg", 0.007
    )
    sched_df["ahu_supply_temperature_cooling"] = schedule_settings.get(
        "ventilation_supply_temp_cooling_c", 16.0
    )
    sched_df["ahu_supply_specific_humidity_cooling"] = schedule_settings.get(
        "ventilation_supply_humidity_cooling_kgkg", 0.009
    )

    sched_df = EndUse.rescale_df(CONFIG.ts_per_hour, sched_df)

    scalar_data = {}
    try:
        scalar_data['conFrac'] = float(schedule_settings['System Convective Fraction'])
        scalar_data['AHUHUM'] = bool(schedule_settings['AHU humidity control'])
        scalar_data['sensRec'] = float(schedule_settings['AHU sensible heat recovery'])
        scalar_data['latRec'] = float(schedule_settings['AHU latent heat recovery'])
        scalar_data['outdoorAirRatio'] = float(schedule_settings['Outdoor Air Ratio'])
        scalar_data['DomesticHotWater calculation'] = str(schedule_settings['DomesticHotWater calculation'])
        scalar_data['Appliances calculation'] = str(schedule_settings['Appliances calculation'])
        
    except KeyError:
        raise KeyError(
            f"ERROR Loading end use {name}. GeneralData does not have the correct columns names: ConvFrac, AHUHum, SensRec, LatRec, OutAirRatio")
    except ValueError:
        raise ValueError(f"""ERROR 
                         Loading end use {name}. GeneralData
                         I'm not able to parse the General data. 
                             ConvFrac should be a float {schedule_settings['ConvFrac']}
                             AHUHum should be a boolean {schedule_settings['AHUHum']}
                             SensRec should be a float {schedule_settings['SensRec']}
                             LatRec  should be a float {schedule_settings['LatRec']}
                             OutAirRatio   should be a float {schedule_settings['OutAirRatio']}
                             DomesticHotWater calculation   should be a str {schedule_settings['DomesticHotWater calculation']}
                             Appliances calculation   should be a str {schedule_settings['Appliances calculation']}
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
    if scalar_data['DomesticHotWater calculation'] not in ["UNI-TS 11300-2", "Schedule", "DHW calc"]:
        raise ValueError(f"End use {name}, DomesticHotWater calculation not allowed: {scalar_data['DomesticHotWater calculation']}. Allowed values: UNI-TS 11300-2, Schedule")
    if scalar_data['Appliances calculation'] not in ["Italian Residential Building Stock", "Schedule"]:
        raise ValueError(f"End use {name}, Appliances calculation not allowed: {scalar_data['Appliances calculation']}. Allowed values: Italian Residential Building Stock, Schedule")

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
        unit='W',
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
        unit='W',
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
        unit='W',
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
        nominal_value=0.2,
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
        unit='Vol/h',
        nominal_value=0.5,
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
        "occupants",
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
        unit="px",
        schedule=domestic_hot_water_sched,
    )

    end_use_obj = EndUse(name)
    end_use_obj.scalar_data = scalar_data

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
        """
        Class method to create an EndUse object from an Excel worksheet.
        
        Parameters
        ----------
        name : str
        Name of the building archetype.
        daily_df_from_excel : pandas.DataFrame
        Hourly schedules for various systems and internal loads.
        scalar_df_from_excel : pandas.DataFrame
        Metadata such as convective fraction, recovery ratios, and calculation methods.
        holidays : list of int
        List of holidays (0–364) used to adjust schedules.
        
        Returns
        -------
        EndUse
        Configured EndUse object with schedule and system properties set.
        """

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
            scalar_data['Appliances calculation'] = str(scalar_df_from_excel.loc[name]['Appliances calculation'])
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
                                 Appliances calculation   should be a str {scalar_df_from_excel.loc[name]['Appliances calculation']}
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
        if scalar_data['DomesticHotWater calculation'] not in ["UNI-TS 11300-2", "Schedule", "DHW calc"]:
            raise ValueError(f"End use {name}, DomesticHotWater calculation not allowed: {scalar_data['DomesticHotWater calculation']}. Allowed values: UNI-TS 11300-2, Schedule")
        if scalar_data['Appliances calculation'] not in ["Italian Residential Building Stock", "Schedule"]:
            raise ValueError(f"End use {name}, Appliances calculation not allowed: {scalar_data['Appliances calculation']}. Allowed values: Italian Residential Building Stock, Schedule")

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
            unit='m3/(m2 s)',
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
        end_use_obj.scalar_data = scalar_data

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

