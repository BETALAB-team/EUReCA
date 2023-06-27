import os
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('TkAgg')
matplotlib.interactive(True)

from eureca_building.config import load_config

config_path = os.path.join('.', 'config.json')
load_config(config_path)
from eureca_building.config import CONFIG

from eureca_building.schedule import Schedule
from eureca_building.domestic_hot_water import DomesticHotWater
from eureca_building.weather import WeatherFile

epw_path = os.path.join('..', 'example_scripts', 'ITA_Venezia-Tessera.161050_IGDG.epw')
weather_file = WeatherFile(epw_path,
                           time_steps=CONFIG.ts_per_hour,
                           azimuth_subdivisions=CONFIG.azimuth_subdivisions,
                           height_subdivisions=CONFIG.height_subdivisions, )

dhw_flow_rate = Schedule(
    "dhw_flow_rate",
    "mass_flow_rate",
    np.array(([.9] * 8 * 2 + [.5] * 2 * 2 + [.9] * 4 * 2 + [.5] * 10 * 2) * 365*6)[:-11],
)

dhw_1 = DomesticHotWater(
    "dhw_1",
    calculation_method="Schedule",
    unit = "L/s",
    schedule=dhw_flow_rate,

)

dhw_2 = DomesticHotWater(
    "dhw_2",
    calculation_method="UNI-TS 11300-2",
)


dhw_3 = DomesticHotWater(
    "dhw_3",
    calculation_method="DHW calc",
)

area = 89.
n_units = 1

volumes = []
demands = []

for i in [dhw_1,dhw_2, dhw_3]:#
    volume, demand = i.get_dhw_yearly_mass_flow_rate(area, n_units, weather_file)
    volumes.append(volume)
    demands.append(demand)

volumes = np.array(volumes).T
demands = np.array(demands).T

fig, ax1 = plt.subplots(nrows = 1)

ax1.plot(volumes)



























