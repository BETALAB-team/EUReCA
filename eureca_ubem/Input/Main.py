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


def main():

    #User Settings
    config_path = os.path.join(".","Example_District_Config.json")                 #Simulation Settings Given as JSON file 
    weather_file = os.path.join(".","ITA_Venezia-Tessera.161050_IGDG.epw")          #Path to weatherfile in epw energyplus format
    schedules_file = os.path.join(".","Schedules_total.xlsx")                       #Path to the schedules for the end use
    materials_file = os.path.join(".","Materials.xlsx")                              #Path to the construction material information
    city_model_file = os.path.join(".","Example_District.geojson")                  #Path to the geoindexed file of the footprints of buildings
    systems_file = os.path.join(".","systems.xlsx")                                 #Path to the HVAC systems specifications
    output_folder = os.path.join(".","output_folder")                                #Path to the output folder
    log_file = os.path.join(output_folder , "run_log.txt")
    
    building_model = "1C"
    shading_calculation = True
    quasi_steady_state = False
    output_type = "csv"

    load_config(config_path) 

    from eureca_ubem.city import City

    os.makedirs(output_folder, exist_ok = True)

    city = City(
        city_model=city_model_file,
        epw_weather_file=weather_file,
        end_uses_types_file=schedules_file,
        envelope_types_file=materials_file,
        systems_templates_file=systems_file,
        shading_calculation=shading_calculation,                                                   #Shading Calculation Requires Preprocessing Time
        building_model = building_model,                                                      #1C for 5R1C (ISO 13790), 2C for 7R2C (VDI6007)
        output_folder=output_folder                            
    )
    if quasi_steady_state:
        city.simulate_quasi_steady_state()
    else:
        city.simulate(output_type=output_type)

    
if __name__ == "__main__":
    main()   


