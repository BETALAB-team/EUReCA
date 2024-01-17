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
materials_file = os.path.join(".","total envelope types.xlsx")
city_model_file = os.path.join(".","InputPD_19Gen2023.geojson")

city_geojson_1 = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
    shading_calculation=True,
    building_model = "1C",
    output_folder=os.path.join(".","results")
)
city_geojson_1.loads_calculation(region="Veneto")
city_geojson_1.simulate(print_single_building_results=True, output_type="csv")
# city_model_file = os.path.join(".","InputPD_.geojson")
#
# city_geojson_2 = City(
#     city_model=city_model_file,
#     epw_weather_file=weather_file,
#     end_uses_types_file=schedules_file,
#     envelope_types_file=materials_file,
#     shading_calculation=True,
#     building_model = "1C",
#     output_folder=os.path.join(".","results_2")
# )
# city_geojson_2.loads_calculation(region="Veneto")
# city_geojson_2.simulate(print_single_building_results=True, output_type="csv")
#
# city_model_file = os.path.join(".","InputPD_small__wo.geojson")
#
# city_geojson_3 = City(
#     city_model=city_model_file,
#     epw_weather_file=weather_file,
#     end_uses_types_file=schedules_file,
#     envelope_types_file=materials_file,
#     shading_calculation=True,
#     building_model = "1C",
#     output_folder=os.path.join(".","results_3")
# )
# city_geojson_3.loads_calculation(region="Veneto")
# city_geojson_3.simulate(print_single_building_results=True, output_type="csv")
