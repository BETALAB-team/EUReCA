''' IMPORTING MODULES '''

import os
import time as tm

import matplotlib
matplotlib.use('TkAgg')
matplotlib.interactive(True)

from eureca_building.config import load_config
load_config("config.json")

from eureca_ubem.city import City

weather_file = os.path.join(".","ITA_Venezia-Tessera.161050_IGDG.epw")
schedules_file = os.path.join(".","Schedules.xlsx")
materials_file = os.path.join(".","materials_and_construction_test.xlsx")
city_model_file = os.path.join(".","PiovegoRestricted_with_holes.geojson")

# Creation of the City object exit
city_geojson = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
)
city_geojson.loads_calculation()


city_model_file = os.path.join(".","PaduaRestricted.json")

# Creation of the City object exit
city_json = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
)
city_json.loads_calculation()

materials_file = os.path.join(".","total envelope types.xlsx")
city_model_file = os.path.join(".","Belzoni.json")

start = tm.time()
# Creation of the City object exit
belzoni = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
)
belzoni.loads_calculation()
belzoni.simulate()

print(f"Belzoni loading and calcs: {(tm.time() - start)/60} min")
