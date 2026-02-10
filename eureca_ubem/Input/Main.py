# Main.py - EUReCA Simulation main file 
import os 
import sys
import logging

#preload heavy dependencies
import numpy as np 
import pandas as pd 
import shapely 
import geopandas as gpd
import pyproj

# Import eureca libraries
from eureca_building.config import load_config
from eureca_ubem import node_calculator, assign_runs
from itertools import product
from time import time


_NO_CHANGE = "__NO_CHANGE__"

def create_scenarios(cityjson: pd.DataFrame, change_dict: dict) -> dict:
    """
    Generates scenarios including cases where some dimensions are not applied at all.
    Returns: dict[str, pd.DataFrame]
    """

    # Keep valid dimensions only
    valid_items = []
    for col, options in change_dict.items():
        if col not in cityjson.columns:
            print(f"column '{col}' is not in the city model, continuing without considering it")
            continue
        if not isinstance(options, (list, tuple)):
            options = [options]
        valid_items.append((col, list(options)))

    if not valid_items:
        return {"scenario_base": cityjson.copy(deep=True)}

    cols = [c for c, _ in valid_items]
    options_lists = [[_NO_CHANGE] + opts for _, opts in valid_items]  # <-- add no-change

    scenarios = {}

    for combo in product(*options_lists):
        # Build key: skip no-change parts
        key_parts = [str(v).replace(" ", "_") for v in combo if v != _NO_CHANGE]
        scenario_key = "scenario_" + ("base" if not key_parts else "_".join(key_parts))

        df = cityjson.copy(deep=True)
        for col, chosen in zip(cols, combo):
            if chosen == _NO_CHANGE:
                continue
            df[col] = chosen

        scenarios[scenario_key] = df

    return scenarios

def main():

    config_path = os.path.join(".","Example_District_Config.json")                 #Simulation Settings Given as JSON file    
    
    weather_file = os.path.join(".","ITA_Venezia-Tessera.161050_IGDG.epw")          #Path to weatherfile in epw energyplus format
    schedules_file = "Schedules_total.xlsx"                     #Path to the schedules for the end use
    materials_file = os.path.join(".","Materials.xlsx")                              #Path to the construction material information
    city_model_file = os.path.join(".","BelzoniDaNEST_ideal.geojson")                  #Path to the geoindexed file of the footprints of buildings


    systems_file = os.path.join(".","systems.xlsx")                                 #Path to the HVAC systems specifications
    output_folder = os.path.join(".","grasloken_check")                                #Path to the output folder
    log_file = os.path.join(output_folder , "run_log.txt")
    building_model = "2C"
    shading_calculation = True
    quasi_steady_state = False
    output_type = "csv"

    load_config(config_path) 
    from eureca_building.config import CONFIG
    from eureca_ubem.city import City

    os.makedirs(output_folder, exist_ok = True)

    city_reference = City(
        city_model=city_model_file,
        epw_weather_file=weather_file,
        end_uses_types_file=schedules_file,
        envelope_types_file=materials_file,
        systems_templates_file=systems_file,
        shading_calculation=shading_calculation,                                                   #Shading Calculation Requires Preprocessing Time
        building_model = building_model,                                                      #1C for 5R1C (ISO 13790), 2C for 7R2C (VDI6007)
        output_folder=output_folder                            
    )
    city_reference.simulate()
    city_reference.dhn_nodes = node_calculator.create_dhn_nodes_from_buildings(city_reference, city_reference)

    # scenario_changes = {"Envelope": ["A" , "B", "C" , "D"],
    #  "Heating system": ["District Heating", "Heat Pump"]}
    # scenarios = create_scenarios(city_reference.cityjson, scenario_changes)
    
    # extreme_scenarios ={}
    # extreme_scenarios["base"] = city_reference
    # for key, value in scenarios.items():
    #     city = City(
    #         city_model=value,
    #         epw_weather_file=weather_file,
    #         end_uses_types_file=city_reference,
    #         envelope_types_file=materials_file,
    #         systems_templates_file=systems_file,
    #         shading_calculation=shading_calculation,                              
    #         building_model = building_model,                                                      
    #         output_folder=os.path.join(output_folder,key)                             
    #     )    
        # extreme_scenarios[key] = city

    # for key, city in extreme_scenarios.items():
    #     city.simulate()
    #     print("simulated")
    #     city.dhn_nodes = node_calculator.create_dhn_nodes_from_buildings(city, city_reference)

    
    # return city_reference, extreme_scenarios
    return city_reference

if __name__ == "__main__":
    mycity = main()   
    a= time()
    
#     #%%
# runs = assign_runs.generate_heatnode_diffusion_runs(scenario_cities=sc,
#     scenario_changes={"Envelope": ["A" , "B", "C" , "D"]}
#      ,
#     scenario_diffusion={
#         "Envelope": [0.1, 0.2, 0.3, 0.3]             
#     },
#     n_runs=500,                
#     rng=np.random.default_rng(42))

# print(time()-a)
