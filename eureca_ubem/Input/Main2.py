''' IMPORTING MODULES '''

import os
import time as tm
import pandas as pd

# import matplotlib
# matplotlib.use('TkAgg')
# matplotlib.interactive(True)

from eureca_building.config import load_config
load_config("config.json")

from eureca_ubem.city import City

weather_file = os.path.join(".","ITA_Venezia-Tessera.161050_IGDG.epw")

schedules_file = os.path.join(".","Schedules_total.xlsx")
materials_file = os.path.join(".","eureca_input_final_blank.xlsx")
city_model_file = os.path.join(".","EURECA_Input_V4_BIS.geojson")

city_geojson = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
)
city_geojson.loads_calculation()
city_geojson.simulate()
