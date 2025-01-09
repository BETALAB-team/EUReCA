''' IMPORTING MODULES '''

import os
import glob
import time as tm
import logging

import pandas as pd
import numpy as np

# import matplotlib
# matplotlib.use('TkAgg')
# matplotlib.interactive(True)

from eureca_building.config import load_config
load_config(os.path.join("eureca_ubem","InputbelzoniPerCoolingDown","config.json"))

from eureca_ubem.city import City

files = glob.glob(os.path.join("eureca_ubem","InputbelzoniPerCoolingDown","*.epw"))

for scen in ["_CS", "_ideal"]:
    for file in files:
        if "REV03b" in file:
            sim_name=file.split(os.sep)[-1][:-4] + scen
            weather_file = file
            schedules_file = os.path.join("eureca_ubem","InputbelzoniPerCoolingDown","Schedules.xlsx")
            materials_file = os.path.join("eureca_ubem","InputbelzoniPerCoolingDown","total envelope types.xlsx")
            city_model_file = os.path.join("eureca_ubem","InputbelzoniPerCoolingDown","BelzoniDaNEST" + scen + ".geojson")
            
            city_geojson = City(
                city_model=city_model_file,
                epw_weather_file=weather_file,
                end_uses_types_file=schedules_file,
                envelope_types_file=materials_file,
                shading_calculation=False,
                building_model = "2C",
                output_folder=os.path.join("eureca_ubem","InputbelzoniPerCoolingDown",sim_name)
            )
            city_geojson.loads_calculation(region="Veneto")
            city_geojson.simulate(print_single_building_results=False, output_type="csv")
            
            t_ext = city_geojson.weather_file.hourly_data["out_air_db_temperature"] 
            
            np.savetxt(os.path.join("eureca_ubem","InputbelzoniPerCoolingDown",sim_name,'t_ext.csv'), t_ext, delimiter=",")

