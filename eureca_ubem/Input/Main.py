''' IMPORTING MODULES '''

import os
import time as tm
import pandas as pd

# import matplotlib
# matplotlib.use('TkAgg')
# matplotlib.interactive(True)

from eureca_building.config import load_config
load_config(os.path.join(".","config2.json"))

from eureca_ubem.city import City

weather_file = os.path.join(".","ITA_Venezia-Tessera.161050_IGDG.epw")

schedules_file = os.path.join(".","Schedules_total.xlsx")
materials_file = os.path.join(".","materials_and_construction_test.xlsx")
city_model_file = os.path.join(".","PiovegoRestricted_with_holes_corr_coef.geojson")
city_model_file = os.path.join(".","PiovegoRestricted_with_holes_corr_coef_sysmod.geojson")
systems_file = os.path.join(".","systems.xlsx")

city_geojson = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
    systems_templates_file=systems_file,
    shading_calculation=True,
    building_model = "2C",
    output_folder=os.path.join(".","geojson_corr_sysmod")
)
city_geojson.loads_calculation(region="Veneto")
city_geojson.simulate( output_type="csv")
# city_geojson.simulate_quasi_steady_state()

