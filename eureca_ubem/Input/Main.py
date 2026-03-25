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
    city_model_file = os.path.join(".","soderman_limited_reproject_2.geojson")                  #Path to the geoindexed file of the footprints of buildings


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
                            baseline_gdf="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject_2.geojson", 
                            weather_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw", 
                            street_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/roads.geojson")

#%%

retrofit, buildings, interv_dict = sc.analyze_intervention(baseline_geojson = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject_2.geojson",
                                               city=mycity,
                                               baseline_scenario=baseline,
                                               weatherfile_path="C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
                                               mode="app",
                                               intervention_dictionary = None)
# q, inter = app.building_editor("C:/Works/EUReCA/EUReCA/eureca_ubem/Input/soderman_limited_reproject.geojson")
# qa = sc.demand_aggregation(q["gdf"], mycity)
# qb = sc.dataframes(qa) 

#%%
import numpy as np
import matplotlib.pyplot as plt

A = retrofit.Electrical_Network[0]
# collect the 8760 arrays from the dictionary objects
arrays = [obj.V for obj in A.nodes]   # list of arrays

# stack -> shape (n_objects, 8760)
data = np.vstack(arrays)

# hourly max and min across all objects
max_vals = data.max(axis=0)
min_vals = data.min(axis=0)

# plotting
hours = np.arange(8760)

plt.figure()
plt.plot(hours, max_vals, label="Max")
plt.plot(hours, min_vals, label="Min")
plt.legend()
plt.xlabel("Hour")
plt.ylabel("Value")
plt.title("Hourly Min/Max Across Objects")
plt.show()
