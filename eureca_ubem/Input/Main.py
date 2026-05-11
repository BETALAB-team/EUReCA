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
from eureca_pubem import dhn_costs as dc

# Import eureca libraries
from eureca_building.config import load_config
from eureca_ubem import node_calculator, assign_runs
from itertools import product
from time import time
from eureca_pubem import scenario_process as sc
import warnings
warnings.filterwarnings("ignore")

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

#%%
from eureca_pubem import scenario_process as sc
from eureca_pubem import buildings_cost as bc

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

#%%
from eureca_pubem import dhn_costs as dc
from eureca_pubem import grid_costs  as gc
from eureca_pubem import market
from eureca_pubem import scenario_process as sc
from eureca_pubem import buildings_cost as bc

assumptions = {
    "discount_rate": 0.04,
    "horizon_years": 20,

    "hp_scop": 3.0,
    "boiler_efficiency": 0.9,

    "pv_yield_kwh_per_m2_year": 160.0,
    "pv_self_consumption_ratio": 0.6,

    "dhn_connection_cost": 30_000.0,

    "pv_capex_per_m2": 2500.0,
    "hp_capex_sh": 90_000.0,
    "hp_capex_dhw": 40_000.0,
    "boiler_capex_sh": 35_000.0,
    "boiler_capex_dhw": 20_000.0,
    "generator_capex": 50_000.0,

    "boiler_fuel_cost_per_kwh": 1.2,
    "generator_cost_per_kwh": 2.5,
    "generator_kwh_per_year": 0.0,

    "electricity_spot_price_per_kwh": 1.0,
}
supplier_costs = {
    "grids": {
        0: {
            "electricity_purchase_cost_per_kwh": 1.0,
            "grid_fixed_cost_yearly": 0.0,
        },
    },
    "dhns": {
        "District Heating Supply 8": {
            "heat_supply_cost_per_mwh": 500.0,
            "fixed_cost_yearly": 0.0,
        },
    },
}
optimization_settings = {
    "max_rounds": 1000,
    "min_rounds": 2,
    "tolerance": 1e-4,

    "local_search": {
        "grid": {
            "active_levers": [
                "pricing.buy.grid local distribution cost monthly fix",
                "pricing.buy.grid local distribution cost monthly per kW peak",
                "pricing.buy.grid local distribution cost monthly per kWh usage",
            ],
            "step_sizes": {
                "pricing.buy.grid local distribution cost monthly fix": 25.0,
                "pricing.buy.grid local distribution cost monthly per kW peak": 2.5,
                "pricing.buy.grid local distribution cost monthly per kWh usage": 0.005,
            },
            "bounds": {
                "pricing.buy.grid local distribution cost monthly fix": (0.0, 3000.0),
                "pricing.buy.grid local distribution cost monthly per kW peak": (0.0, 500.0),
                "pricing.buy.grid local distribution cost monthly per kWh usage": (0.0, 2.0),
            },
        },

        "dhn": {
            "active_levers": [
                "pricing.area fee per m2",
                "pricing.fixed heat price per MWh",
                "pricing.variable heat price per MWh",
                "pricing.admin fee yearly",
                "pricing.subscription fixed yearly per unit",
                "pricing.subscription variable price per MWh",
                "connection_cost",
            ],
            "step_sizes": {
                "pricing.area fee per m2": 1.0,
                "pricing.fixed heat price per MWh": 5.0,
                "pricing.variable heat price per MWh": 10.0,
                "pricing.admin fee yearly": 50.0,
                "pricing.subscription fixed yearly per unit": 100.0,
                "pricing.subscription variable price per MWh": 5.0,
                "connection_cost": 1000.0,
            },
            "bounds": {
                "pricing.area fee per m2": (0.0, 200.0),
                "pricing.fixed heat price per MWh": (0.0, 1000.0),
                "pricing.variable heat price per MWh": (0.0, 5000.0),
                "pricing.admin fee yearly": (0.0, 10000.0),
                "pricing.subscription fixed yearly per unit": (0.0, 50000.0),
                "pricing.subscription variable price per MWh": (0.0, 1000.0),
                "connection_cost": (0.0, 50000.0),
            },
        },
    },
}

grid_regulation = {
    "allowed_revenue_yearly": 8_000_000.0,
    "background_similarity_factor": 0.30,

    "background_baseline_customer_count": 8_000.0,
    "background_baseline_electricity_bought_kwh_year": 32_000_000.0,
    "background_baseline_electricity_sold_kwh_year": 2_000_000.0,
    "background_baseline_peak_kw": 12_000.0,

    "revenue_gap_weight": 1.0,
    "tariff_change_weight": 1_000.0,
}
zero_study_area_baseline_grid_totals = {
    "customer_count": 0.0,
    "electricity_bought_kwh_year": 0.0,
    "electricity_sold_kwh_year": 0.0,
    "peak_kw": 0.0,
}





#%%
from eureca_pubem import scenario_process as sc

#Initializing

# 0.0. baseline generation 
baseline = sc.create_baseline(input_gdf=mycity_gdf, 
                            input_city=mycity, 
                            distributions_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/distribution.geojson",
                            baseline_gdf="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson", 
                            weather_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw", 
                            street_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/roads.geojson")
baseline_gdf = load_buildings("C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson")
config0 = bc.extract_config((baseline_gdf))
interv_dict, building_info = sc.make_dictionary(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                                               city=mycity,
                                               baseline_scenario=baseline,
                                               weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                                               intervention_dictionary = config0)
retrofit, buildings, _ = sc.analyze_intervention(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                                               city=mycity,
                                               baseline_scenario=baseline,
                                               weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                                               intervention_dictionary = interv_dict)
retro_1 = retrofit
#%%
# 0.1. system capital cost initial
# 0.1.1 dhn capital cost initial 
dhn_pipe_costs = dc.build_cost_table(pipe_data = "C:\Works\EUReCA\EUReCA\eureca_pubem\dhn\_pipe_diameters.json")
capital_costs_dhn = dc.compute_dhn_cost(
        dhn_pipe_changes = retrofit.dhn_pipe_changes,
        area_type = "urban",
        assumptions = {
        "replacement_factor_ground": 0.5,
        "removal_factor": 0.3
    },
        pipe_json = dhn_pipe_costs
        )
capital_costs_grid = gc.compute_grid_cost(
    grid_line_changes=retrofit.grid_line_changes,
    area_type="town",
    assumptions={
        "replacement_factor_ground": 1.0,
        "removal_factor": 0.2,
        "ground_share": 0.55,
        "rest_share": 0.45
    },
    cable_json="C:/Works/EUReCA/EUReCA/eureca_pubem/grid/swedish_cable.json",
    cable_key="name"
)

# 0.1.2 grid capital cost initial 

# 0.2. initial market 
buildings = market.reachables (buildings, retrofit.District_Heating_Systems, retrofit.Electrical_Network)

buildings = market.build_network_maps(buildings, 
                          retrofit.District_Heating_Systems,
                          retrofit.Electrical_Network)
dhn_levers, grid_levers = market.initialize_market_levers(dhns = retrofit.District_Heating_Systems, grids=retrofit.Electrical_Network)
building_option_cache = market.build_building_option_cache(
    buildings=buildings,
    grids=retrofit.Electrical_Network,
    dhns=retrofit.District_Heating_Systems,
    grid_levers=grid_levers,
    dhn_levers=dhn_levers,
    assumptions = assumptions
)
tech_data = market.prepare_global_assumptions(assumptions)
avg_n_occ = sum(d["meta"]["n_occ"] for d in buildings.values()) / len(buildings)
market_results = market.optimize_market_levers(
    building_option_cache=building_option_cache,
    initial_grid_levers=grid_levers,
    initial_dhn_levers=dhn_levers,
    tech_data=tech_data,
    supplier_costs=supplier_costs,
    optimization_settings=optimization_settings,
    grid_regulation=grid_regulation,
    zero_study_area_baseline_grid_totals=zero_study_area_baseline_grid_totals,
    avg_n_occ=avg_n_occ,
    study_area_grid_capex=capital_costs_grid,
)



#0.3. first building move 
grid_pricing = {"1":market_results["grid_levers"][0]["pricing"]}
dhn_pricing = {"1":{"buy": market_results["dhn_levers"][next(iter(market_results["dhn_levers"]))]["pricing"]}}

Building_dict = bc.build_dict_gen(building_info,
                   interv_dict,
                   baseline_gdf_path = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                   configuration = config0, 
                   ee_measure_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\EE_measures_catalog.json",
                   pv_type_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\pv_config.json",
                   hp_catalog_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\hp_config.json",
                   grid_pricing_path = grid_pricing,
                   dhn_pricing_path = dhn_pricing,
                   fuels_path =r"C:\Works\EUReCA\EUReCA\eureca_pubem\fuels.json",
                   spot_price_path="C:\Works\EUReCA\EUReCA\eureca_pubem\spot_price_se3.csv"
                   )


current_optimal_config, current_optimal_set = bc.optimize_configuration_per_building_one_step(config_current = config0,
                              current_dictionary=Building_dict,
                                baseline_dictionary = Building_dict,
                               baseline_gdf_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                               ee_measure_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\EE_measures_catalog.json",
                               pv_type_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\pv_config.json",
                               hp_catalog_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\hp_config.json",
                               grid_pricing_path = grid_pricing,
                               dhn_pricing_path = dhn_pricing,
                               fuels_path =r"C:\Works\EUReCA\EUReCA\eureca_pubem\fuels.json",
                               spot_price_path="C:\Works\EUReCA\EUReCA\eureca_pubem\spot_price_se3.csv",
                               weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                               mycity=mycity,
                               baseline_scenario=baseline,
                               r = 0.04, 
                               T = 25
                               )
#%%
buildings0 = buildings
import copy
previous_config = current_optimal_config
def diff_dicts(old, new, path=""):
    changes = []

    old_keys = set(old.keys())
    new_keys = set(new.keys())

    for key in old_keys - new_keys:
        changes.append((f"{path}{key}", old[key], "__MISSING__"))

    for key in new_keys - old_keys:
        changes.append((f"{path}{key}", "__MISSING__", new[key]))

    for key in old_keys & new_keys:
        old_val = old[key]
        new_val = new[key]
        new_path = f"{path}{key}"

        if isinstance(old_val, dict) and isinstance(new_val, dict):
            changes.extend(diff_dicts(old_val, new_val, path=new_path + "."))
        elif old_val != new_val:
            changes.append((new_path, old_val, new_val))

    return changes
for game_step in range(0,20):

    now_config = copy.deepcopy(current_optimal_config)
    # changes = diff_dicts(previous_config, now_config)

    # if not changes:
    #     print(f"Stopping at game_step={game_step}: no building changed configuration.")
    #     current_optimal_config = copy.deepcopy(now_config)
    #     break
    
    # print(f"\nChanges at game_step={game_step}:")
    # for key_path, old_value, new_value in changes:
    #     print(f"  {key_path}: {old_value} -> {new_value}")
    
    # current_optimal_config = copy.deepcopy(now_config)

    # if previous_config is not None and now_config == previous_config:
    #     print(f"Stopping at game_step={game_step}: configuration did not change.")
    #     break
    
    interv_dict, building_info = sc.make_dictionary(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                                                   city=mycity,
                                                   baseline_scenario=baseline,
                                                   weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                                                   intervention_dictionary = now_config)
    retrofit, buildings, _ = sc.analyze_intervention(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                                                   city=mycity,
                                                   baseline_scenario=baseline,
                                                   weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                                                   intervention_dictionary = interv_dict)
    
    retro_n = retrofit 
    # 0.1. system capital cost initial
    # 0.1.1 dhn capital cost initial 
    dhn_pipe_costs = dc.build_cost_table(pipe_data = "C:\Works\EUReCA\EUReCA\eureca_pubem\dhn\_pipe_diameters.json")
    capital_costs_dhn = dc.compute_dhn_cost(
            dhn_pipe_changes = retrofit.dhn_pipe_changes,
            area_type = "urban",
            assumptions = {
            "replacement_factor_ground": 0.5,
            "removal_factor": 0.3
        },
            pipe_json = dhn_pipe_costs
            )
    capital_costs_grid = gc.compute_grid_cost(
        grid_line_changes=retrofit.grid_line_changes,
        area_type="town",
        assumptions={
            "replacement_factor_ground": 1.0,
            "removal_factor": 0.2,
            "ground_share": 0.55,
            "rest_share": 0.45
        },
        cable_json="C:/Works/EUReCA/EUReCA/eureca_pubem/grid/swedish_cable.json",
        cable_key="name"
    )
    # 0.1.2 grid capital cost initial 
    
    # 0.2. initial market 
    buildings = market.reachables (buildings, retrofit.District_Heating_Systems, retrofit.Electrical_Network)
    
    buildings = market.build_network_maps(buildings, 
                              retrofit.District_Heating_Systems,
                              retrofit.Electrical_Network)
    dhn_levers, grid_levers = market.initialize_market_levers(dhns = retrofit.District_Heating_Systems, grids=retrofit.Electrical_Network)
    building_option_cache = market.build_building_option_cache(
        buildings=buildings,
        grids=retrofit.Electrical_Network,
        dhns=retrofit.District_Heating_Systems,
        grid_levers=grid_levers,
        dhn_levers=dhn_levers,
        assumptions = assumptions
    )
    tech_data = market.prepare_global_assumptions(assumptions)
    avg_n_occ = sum(d["meta"]["n_occ"] for d in buildings.values()) / len(buildings)
    market_results = market.optimize_market_levers(
        building_option_cache=building_option_cache,
        initial_grid_levers=grid_levers,
        initial_dhn_levers=dhn_levers,
        tech_data=tech_data,
        supplier_costs=supplier_costs,
        optimization_settings=optimization_settings,
        grid_regulation=grid_regulation,
        zero_study_area_baseline_grid_totals=zero_study_area_baseline_grid_totals,
        avg_n_occ=avg_n_occ,
        study_area_grid_capex=capital_costs_grid,
    )
    
    #0.3. first building move 
    grid_pricing = {"1":market_results["grid_levers"][0]["pricing"]}
    dhn_pricing = {"1":{"buy": market_results["dhn_levers"][next(iter(market_results["dhn_levers"]))]["pricing"]}}
    print(grid_pricing)
    print(dhn_pricing)
    
    current_dict = bc.build_dict_gen(building_info,
                       interv_dict,
                       baseline_gdf_path = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                       configuration = now_config, 
                       ee_measure_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\EE_measures_catalog.json",
                       pv_type_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\pv_config.json",
                       hp_catalog_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\hp_config.json",
                       grid_pricing_path = grid_pricing,
                       dhn_pricing_path = dhn_pricing,
                       fuels_path =r"C:\Works\EUReCA\EUReCA\eureca_pubem\fuels.json",
                       spot_price_path="C:\Works\EUReCA\EUReCA\eureca_pubem\spot_price_se3.csv"
                       )
    
    previous_config = copy.deepcopy(current_optimal_config)
    current_optimal_config, current_optimal_set = bc.optimize_configuration_per_building_one_step(config_current = now_config,
                                  current_dictionary=current_dict,
                                    baseline_dictionary = Building_dict,
                                   baseline_gdf_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
                                   ee_measure_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\EE_measures_catalog.json",
                                   pv_type_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\pv_config.json",
                                   hp_catalog_path = "C:\Works\EUReCA\EUReCA\eureca_pubem\hp_config.json",
                                   grid_pricing_path = grid_pricing,
                                   dhn_pricing_path = dhn_pricing,
                                   fuels_path =r"C:\Works\EUReCA\EUReCA\eureca_pubem\fuels.json",
                                   spot_price_path="C:\Works\EUReCA\EUReCA\eureca_pubem\spot_price_se3.csv",
                                   weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                                   mycity=mycity,
                                   baseline_scenario=baseline,
                                   r = 0.05, 
                                   T = 50
                                   )
    now_config = copy.deepcopy(current_optimal_config)
    changes = diff_dicts(previous_config, now_config)

    if not changes:
        print(f"Stopping at game_step={game_step}: no building changed configuration.")
        current_optimal_config = copy.deepcopy(now_config)
        break
    
    print(f"\nChanges at game_step={game_step}:")
    for key_path, old_value, new_value in changes:
        print(f"  {key_path}: {old_value} -> {new_value}")
    
    current_optimal_config = copy.deepcopy(now_config)

    if previous_config is not None and now_config == previous_config:
        print(f"Stopping at game_step={game_step}: configuration did not change.")
        break
    
    
# #%%
# from eureca_pubem import scenario_process as sc

# now_config = current_optimal_config 
# interv_dict, building_info = sc.make_dictionary(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
#                                                city=mycity,
#                                                baseline_scenario=baseline,
#                                                weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
#                                                intervention_dictionary = now_config)
# retrofit, buildings, _ = sc.analyze_intervention(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson",
#                                                city=mycity,
#                                                baseline_scenario=baseline,
#                                                weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
#                                                intervention_dictionary = interv_dict)
