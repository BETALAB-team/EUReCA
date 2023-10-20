''' IMPORTING MODULES '''

import os
import time as tm
import logging

import pandas as pd

# import matplotlib
# matplotlib.use('TkAgg')
# matplotlib.interactive(True)

from eureca_building.config import load_config
load_config("config.json")

from eureca_ubem.city import City

weather_file = os.path.join(".","ITA_Venezia-Tessera.161050_IGDG.epw")
schedules_file = os.path.join(".","Schedules.xlsx")
materials_file = os.path.join(".","materials_and_construction_test.xlsx")
city_model_file = os.path.join(".","PiovegoRestricted_with_holes.geojson")
#
city_geojson = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
    shading_calculation=True,
    building_model = "1C",
    output_folder=os.path.join(".","geojson")
)
city_geojson.loads_calculation(region="Veneto")
city_geojson.simulate(print_single_building_results=True)
#
# materials_file = os.path.join(".","total envelope types.xlsx")
# city_model_file = os.path.join(".","Belzoni_2023_July_Update.json")
#
# start = tm.time()
# # Creation of the City object exit
# belzoni = City(
#     city_model=city_model_file,
#     epw_weather_file=weather_file,
#     end_uses_types_file=schedules_file,
#     envelope_types_file=materials_file,
#     output_folder=os.path.join(".","belzoni_new_hvac_DHW_calc")
# )
# print(f"Belzoni creation : {(tm.time() - start)/60:.2f} min")
# start = tm.time()
# belzoni.loads_calculation(region="Veneto")
# print(f"Belzoni loads calc : {(tm.time() - start)/60:0.2f} min")
# start = tm.time()
# belzoni.simulate(print_single_building_results=True)
# print(f"Belzoni simulation : {(tm.time() - start)/60:0.2f} min")

# city_model_file = os.path.join(".","PaduaRestricted.json")
#
# # Creation of the City object exit
# city_json = City(
#     city_model=city_model_file,
#     epw_weather_file=weather_file,
#     end_uses_types_file=schedules_file,
#     envelope_types_file=materials_file,
#     output_folder=os.path.join(".","cityjson")
# )
# city_json.loads_calculation()


# import json
# import numpy as np
#
# with open(".\\Belzoni.json", "r") as file:
#     cityjson = json.load(file)
#
# for bd in cityjson["CityObjects"].values():
#     if bd["type"] == "Building":
#         heating_residential = np.random.choice(np.arange(1, 5), p=[0.4, 0.5, 0.05, 0.05])
#         heating_residential = {
#             1:"Traditional Gas Boiler, Centralized, High Temp Radiator",
#             2:"Condensing Gas Boiler, Centralized, Low Temp Radiator",
#             3:"Oil Boiler, Centralized, High Temp Radiator",
#             4:"A-W Heat Pump, Single, Fan coil",
#         }[heating_residential]
#         heating_not_residential = np.random.choice(np.arange(1, 3), p=[0.9, 0.1])
#         heating_not_residential = {
#             1:"Condensing Gas Boiler, Centralized, Fan coil",
#             2:"A-W Heat Pump, Centralized, Radiant surface",
#         }[heating_not_residential]
#
#         if bd["attributes"]["End Use"] == "Residential":
#             bd["attributes"]["Heating System"] = heating_residential
#         else:
#             bd["attributes"]["Heating System"] = heating_not_residential
#
#         cooling_residential = np.random.choice(np.arange(1, 3), p=[0.95,0.05])
#         cooling_residential = {
#             1:"A-A split",
#             2:"A-W chiller, Centralized, Fan coil",
#         }[cooling_residential]
#         cooling_not_residential = np.random.choice(np.arange(1, 3), p=[0.1, 0.9])
#         cooling_not_residential = {
#             1:"A-A split",
#             2:"A-W chiller, Centralized, Fan coil",
#         }[cooling_not_residential]
#
#         if bd["attributes"]["End Use"] == "Residential":
#             bd["attributes"]["Cooling System"] = cooling_residential
#         else:
#             bd["attributes"]["Cooling System"] = cooling_not_residential
#
#
# with open(".\\Belzoni_new_hvac.json", "w") as file:
#     json.dump(cityjson,file)

