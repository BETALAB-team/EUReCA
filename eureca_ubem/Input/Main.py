# Main.py - EUReCA Simulation main file 
import os 
# os.environ["PYTHONUTF8"] = "1"
# os.environ["PYTHONIOENCODING"] = "utf-8"

# set PYTHONUTF8=0
# set PYTHONIOENCODING=cp1257
import sys
import logging

#preload heavy dependencies
import numpy as np 
import geopandas as gpd

# Import eureca libraries
from eureca_building.config import load_config
from eureca_ubem import node_calculator, assign_runs
from itertools import product
from time import time
from eureca_pubem import scenario_process as sc

_NO_CHANGE = "__NO_CHANGE__"




def main():

    config_path = os.path.join(".","Example_District_Config.json")                 #Simulation Settings Given as JSON file    
    
    weather_file = os.path.join(".","SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw")          #Path to weatherfile in epw energyplus format
    schedules_file = "TransitionMatrixSchedule.csv"                     #Path to the schedules for the end use
    materials_file = os.path.join(".","tabula_sverige.xlsx")                              #Path to the construction material information
    city_model_file = os.path.join(".","soderman_limited_reproject.geojson")                  #Path to the geoindexed file of the footprints of buildings


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
    
    cities_envelopes = {}
    
    for i, x in enumerate(["none", "shallow", "medium", "deep"]):
        os.makedirs(output_folder, exist_ok = True)
        gdf = sc.load_geojson_with_envelope_prefix(city_model_file, i)
        scenario_city = City(
            city_model=gdf,
            epw_weather_file=weather_file,
            end_uses_types_file=schedules_file,
            envelope_types_file=materials_file,
            systems_templates_file=systems_file,
            shading_calculation=shading_calculation,                                                   #Shading Calculation Requires Preprocessing Time
            building_model = building_model,                                                      #1C for 5R1C (ISO 13790), 2C for 7R2C (VDI6007)
            output_folder=output_folder                            
        )
        scenario_city.simulate()
        gdf, value_store = sc.create_gdf_dictionary(scenario_city)
        cities_envelopes[x] = value_store


    return gdf, cities_envelopes

if __name__ == "__main__":
    mycity_gdf, mycity = main()   
    a= time()
    
#%%
from eureca_pubem import scenario_process as sc

import warnings
warnings.filterwarnings("ignore")
baseline = sc.create_baseline(input_gdf=mycity_gdf, 
                            input_city=mycity, 
                            distributions_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/distribution.geojson",
                            baseline_gdf="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson", 
                            weather_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw", 
                            street_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/roads.geojson")

#%%
from eureca_pubem import scenario_process as sc
from eureca_pubem import buildings_cost as bc
# config = {
#     1:  {"env": "deep",   "heat": "boiler", "dhw": "hp_he", "fuel": "gas", "pv_percentage": 60},
#     2:  {"env": "shallow","heat": "hp_me", "dhw": "hp_me", "fuel": None, "pv_percentage": 60},
#     3:  {"env": "shallow","heat": "dhn", "dhw": "hp_me", "fuel": None, "pv_percentage": 60},
#     4:  {"env": "deep",   "heat": "dhn", "dhw": "dhn", "fuel": None, "pv_percentage": 60},
#     5:  {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 70},
#     7:  {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     25: {"env": "deep",   "heat": "hp_he", "dhw": "dhn", "fuel": None, "pv_percentage": 60},
#     35: {"env": "shallow","heat": "dhn", "dhw": "hp_me", "fuel": None, "pv_percentage": 60},
#     36: {"env": "deep",   "heat": "dhn", "dhw": "hp_he", "fuel": None, "pv_percentage": 70},
#     38: {"env": "shallow","heat": "hp_me", "dhw": "dhn", "fuel": None, "pv_percentage": 60},
#     40: {"env": "medium", "heat": "dhn", "dhw": "hp_he", "fuel": None, "pv_percentage": 70},
#     41: {"env": "medium", "heat": "dhn", "dhw": "dhn", "fuel": None, "pv_percentage": 70},
#     45: {"env": "medium", "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     47: {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     48: {"env": "medium", "heat": "hp_me", "dhw": "hp_me", "fuel": None, "pv_percentage": 60},
#     52: {"env": "deep",   "heat": "dhn", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     65: {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     77: {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
# }
def load_buildings(input_data):

    if isinstance(input_data, gpd.GeoDataFrame):
        gdf = input_data.copy()
    elif isinstance(input_data, str):
        gdf = gpd.read_file(input_data)
    else:
        raise ValueError("invalid_input")

    required_columns = [
        "Name",
        "EEdepth",
        "SHSource",
        "DHWsource",
        "PVType",
        "PVpercentage"
    ]

    missing = [c for c in required_columns if c not in gdf.columns]

    if missing:
        raise ValueError(f"{missing}")

    return gdf

baseline_gdf = load_buildings("C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson")
config0 = bc.extract_config((baseline_gdf))

interv_dict, building_info = sc.make_dictionary(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                                               city=mycity,
                                               baseline_scenario=baseline,
                                               weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                                               intervention_dictionary = config0)
# retrofit, buildings, _ = sc.analyze_intervention(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
#                                                city=mycity,
#                                                baseline_scenario=baseline,
#                                                weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
#                                                intervention_dictionary = config)


Building_dict = bc.build_dict_gen(building_info,
                   interv_dict,
                   baseline_gdf_path = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                   configuration = config0, 
                   ee_measure_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\EE_measures_catalog.json",
                   pv_type_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\pv_config.json",
                   hp_catalog_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\hp_config.json",
                   grid_pricing_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\grid_pricing_formulas.json",
                   dhn_pricing_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\dhn_pricing_formulas.json",
                   fuels_path =r"C:\Works\EUReCA\EUReCA\eureca_pubem\fuels.json",
                   spot_price_path="C:\Works\EUReCA\EUReCA\eureca_pubem\spot_price_se3.csv"
                   )

Apparent_optimal_config, Apparent_optimal_set = bc.optimize_configuration_per_building(config0 = config0,
                                base_dictionary = Building_dict,
                               baseline_gdf_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                               ee_measure_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\EE_measures_catalog.json",
                               pv_type_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\pv_config.json",
                               hp_catalog_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\hp_config.json",
                               grid_pricing_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\grid_pricing_formulas.json",
                               dhn_pricing_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\dhn_pricing_formulas.json",
                               fuels_path =r"C:\Works\EUReCA\EUReCA\eureca_pubem\fuels.json",
                               spot_price_path="C:\Works\EUReCA\EUReCA\eureca_pubem\spot_price_se3.csv",
                               weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                               mycity=mycity,
                               baseline_scenario=baseline,
                               r = 0.04, 
                               T = 25
                               )
retrofit, buildings, _ = sc.analyze_intervention(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                                               city=mycity,
                                               baseline_scenario=baseline,
                                               weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                                               intervention_dictionary = Apparent_optimal_config)


# config = {
#     1:  {"env": "deep",   "heat": "boiler", "dhw": "hp_he", "fuel": "gas", "pv_percentage": 60},
#     2:  {"env": "shallow","heat": "hp_me", "dhw": "hp_me", "fuel": None, "pv_percentage": 60},
#     3:  {"env": "shallow","heat": "dhn", "dhw": "hp_me", "fuel": None, "pv_percentage": 60},
#     4:  {"env": "deep",   "heat": "dhn", "dhw": "dhn", "fuel": None, "pv_percentage": 60},
#     5:  {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 70},
#     7:  {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     25: {"env": "deep",   "heat": "hp_he", "dhw": "dhn", "fuel": None, "pv_percentage": 60},
#     35: {"env": "shallow","heat": "dhn", "dhw": "hp_me", "fuel": None, "pv_percentage": 60},
#     36: {"env": "deep",   "heat": "dhn", "dhw": "hp_he", "fuel": None, "pv_percentage": 70},
#     38: {"env": "shallow","heat": "hp_me", "dhw": "dhn", "fuel": None, "pv_percentage": 60},
#     40: {"env": "medium", "heat": "dhn", "dhw": "hp_he", "fuel": None, "pv_percentage": 70},
#     41: {"env": "medium", "heat": "dhn", "dhw": "dhn", "fuel": None, "pv_percentage": 70},
#     45: {"env": "medium", "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     47: {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     48: {"env": "medium", "heat": "hp_me", "dhw": "hp_me", "fuel": None, "pv_percentage": 60},
#     52: {"env": "deep",   "heat": "dhn", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     65: {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
#     77: {"env": "deep",   "heat": "hp_he", "dhw": "hp_he", "fuel": None, "pv_percentage": 60},
# }
# New = bc.compare_config_with_base(Building_dict,
#                                baseline_gdf_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
#                                configuration = config, 
#                                ee_measure_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\EE_measures_catalog.json",
#                                pv_type_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\pv_config.json",
#                                hp_catalog_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\hp_config.json",
#                                grid_pricing_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\grid_pricing_formulas.json",
#                                dhn_pricing_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\dhn_pricing_formulas.json",
#                                fuels_path =r"C:\Works\EUReCA\EUReCA\eureca_pubem\fuels.json",
#                                spot_price_path="C:\Works\EUReCA\EUReCA\eureca_pubem\spot_price_se3.csv",
#                                weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
#                                mycity=mycity,
#                                baseline_scenario=baseline
#                                )
# import json
# import pandas as pd
# with open(r"C:\Works\EUReCA\EUReCA\eureca_pubem\EE_measures_catalog.json", "r") as f:
#     ee_measure = json.load(f)
# with open(r"C:\Works\EUReCA\EUReCA\eureca_pubem\pv_config.json", "r") as f:
#     pv_installs = json.load(f)
# with open(r"C:\Works\EUReCA\EUReCA\eureca_pubem\hp_config.json", "r") as f:
#     hp_catalog = json.load(f)
# with open(r"C:\Works\EUReCA\EUReCA\eureca_pubem\grid_pricing_formulas.json", "r") as f:
#     pricing = json.load(f)
    
# with open(r"C:\Works\EUReCA\EUReCA\eureca_pubem\dhn_pricing_formulas.json", "r") as f:
#     dnn_pricing = json.load(f)
    
# with open(r"C:\Works\EUReCA\EUReCA\eureca_pubem\fuels.json", "r") as f:
#     fuels = json.load(f)
# baseline_gdf = load_buildings("C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson")

    
# prices = pd.read_csv(r"C:\Works\EUReCA\EUReCA\eureca_pubem\spot_price_se3.csv") 
# prices["price"] = prices["price"].astype(float)


# Buildings_dict = {}

# for idx, building in a.items():
#     EUR_to_SEK = 10.86
#     meta={}
#     interventions = {}
#     operation = {}
#     costs = {}
#     idx = int(idx)
#     if idx in interv_dict.keys():
#         interventions = interv_dict[idx]
#     meta["PV_type"] = baseline_gdf.loc[baseline_gdf["id"] == idx, "PVType"].values[0]
#     pv_production = building ["pv_production"]
#     electricity_need = building["hp_electricity"] + building["base"]["appliance_electricity"]
#     operation["electricity_bought"] = np.maximum(electricity_need - pv_production, 0)
#     operation["electricity_sold"] = np.maximum(pv_production - electricity_need , 0)
#     operation["dhn_bought"] = building["thermal"]["dhn_demand"]
#     fuel = config[idx]["fuel"]
#     if fuel != None:
#         operation["fuel_bought"] = building["thermal"][f"{fuel}_demand"]
#         meta["fuel"]=fuel
        
#     area = building["base"]["pv_available_area"]
#     floors = baseline_gdf.loc[baseline_gdf["id"] == idx, "Floors"]
#     used_area = area*floors

#     efficiency_measure_cost = 0
#     pv_install_cost = 0
#     envelope_type = baseline_gdf.loc[baseline_gdf["id"] == idx, "Envelope"].values[0]
#     building_type = envelope_type[:3]
#     costs["capital cost"]={}
#     costs["operational cost"] = {}
    
#     if 'envelope' in interv_dict[idx].keys():
#         efficiency_measure = interv_dict[idx]['envelope']
#         before_measure = efficiency_measure[0]
#         if before_measure in ["medium", "deep"]:
#             before_measure = before_measure + " " + building_type
#         after_measure = efficiency_measure[1]
#         if after_measure in ["medium", "deep"]:
#             after_measure = after_measure + " " + building_type
#         wall_area = building["base"]["opaque_exposed_area"] - building["base"]["pv_available_area"]
#         roof_area = building["base"]["pv_available_area"]
#         window_area = building ["base"]["glazing_area"]
#         wall_cost_fix=ee_measure[after_measure]["wall cost constant"] - ee_measure[before_measure]["wall cost constant"]
#         roof_cost_fix=ee_measure[after_measure]["roof cost constant"] - ee_measure[before_measure]["roof cost constant"]
#         window_cost_fix=ee_measure[after_measure]["window cost constant"] - ee_measure[before_measure]["window cost constant"]
#         wall_cost_persqm=ee_measure[after_measure]["roof cost per square meter"] - ee_measure[before_measure]["roof cost per square meter"]
#         roof_cost_persqm=ee_measure[after_measure]["wall cost per square meter"] - ee_measure[before_measure]["wall cost per square meter"]
#         window_cost_persqm=ee_measure[after_measure]["window cost per square meter"] - ee_measure[before_measure]["window cost per square meter"]
#         efficiency_measure_cost = wall_cost_fix + roof_cost_fix + window_cost_fix \
#                                 + wall_cost_persqm * wall_area \
#                                 + window_cost_persqm * window_area\
#                                 + roof_cost_persqm * roof_area
                                
#     if "PVpercentage" in interv_dict[idx].keys():
#         pv_install = interv_dict[idx]["PVpercentage"]
#         before_measure = pv_install[0]
#         after_measure = pv_install[1]
#         pv_type = meta["PV_type"]
#         pv_install_area = (after_measure - before_measure) * building["base"]["pv_available_area"]/100
#         pv_install_fix_cost = pv_installs[pv_type]["cost_fixed"]
#         pv_install_persqm_cost = pv_installs[pv_type]["cost_per_m2"]
#         pv_install_cost = pv_install_persqm_cost * pv_install_area + pv_install_fix_cost
        
#     if any(x in interv_dict[idx].keys() for x in ["dhw_source","sh_source"]):

#         hp_install_cost = hp_cost(building, interv_dict[idx], hp_catalog, factor=1.0) * EUR_to_SEK 
#         dhn_connection_cost =  dhn_cost(building, interv_dict[idx], connection_cost = 24000)
        

#     costs["capital cost"]["efficiency_measure_cost"]  =   efficiency_measure_cost     
#     costs["capital cost"]["pv_install_cost"]  =   pv_install_cost      
#     costs["capital cost"]["heating_systems"] = hp_install_cost + dhn_connection_cost 
#     costs["capital cost"]["total"] = sum(list(costs["capital cost"].values()))
#     costs["operational cost"]["electricity"] = np.pad(hourly_grid_cost(operation, prices, pricing["1"]), (8760 - 8759, 0), mode='edge')
#     costs["operational cost"]["electricity"] = np.pad(costs["operational cost"]["electricity"], (8760 - len(costs["operational cost"]["electricity"]), 0), mode='edge')
#     costs["operational cost"]["district_heating"] = hourly_dhn_cost(operation, dnn_pricing["1"]["buy"], used_area)
#     if fuel != None:
#         costs["operational cost"]["fuel"] = fuel_cost_array(meta["fuel"], operation["fuel_bought"], fuels)*EUR_to_SEK
#     costs["operational cost"]["hourly_total"] = np.sum(list(costs["operational cost"].values()), axis=0)
#     costs["operational cost"]["yearly_total"] = np.sum(costs["operational cost"]["hourly_total"])
        
#     Buildings_dict[idx] = {"meta": meta,
#                            "interventions" : interventions,
#                            "operation" : operation,
#                            "costs": costs}
        
        
        
   
    






#%%
# Building 
# interventions 
# operation 

# DH 
# intervention
# operation 
# selling pricing 

# Grid 
# intervention 
# operation 
# selling pricing 
# buying pricing 


