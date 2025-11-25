# -*- coding: utf-8 -*-
"""

---
Created by: Mohamad
Betalab - DII, University of Padua
---

"""

import pandas as pd
import numpy as np
from scipy.ndimage import shift
import pvlib
from eureca_building.weather import WeatherFile
# from eureca_building.config import CONFIG
#class
area=10000
W_annual=area*3600*area**(-0.18)
W_dailyy=W_annual/(365)
normal_annual_consumption = 230 #kWh/m2/a
Load_ratio = pd.DataFrame({
    'hour': np.arange(25),
    'medium_temperature': 
    [0.988982241,0.926722733,0.856251752,0.819348731,0.864454333,
     0.841634661,0.773614634,0.849485855,1.014482565,1.037806412,
     0.956058409,0.987516480,1.100522635,1.104290722,1.034807861,
     1.058795741,1.139486332,1.121850772,1.020394755,0.984865811,
     1.068386570,1.165873083,1.153045516,1.049329578,0.988982241],
    'low_temperature':
    [1.021834447,1.130729321,0.941327968,0.863925212,0.889266827,
     0.841423316,0.897480049,0.867075202,0.917750535,1.331083444,
     1.088187006,1.034058157,1.207069106,0.925160330,0.981194843,
     0.891645807,0.877585329,1.106916961,1.036566668,0.916558151,
     0.997308728,1.193684741,0.994156594,1.014320512,1.021834447]
})
dataframe=pd.DataFrame({
    "hour":[1,13.25,18.5],
    "temp":[2,14,18.2],
    "sched":[0,1,1]
    })   
FR=0.7
Med_LR = np.interp(dataframe["hour"], Load_ratio["hour"], Load_ratio["medium_temperature"])
Low_LR = np.interp(dataframe["hour"], Load_ratio["hour"], Load_ratio["low_temperature"])
a=np.array(dataframe["sched"]) * 0.75 + np.array(1 - dataframe["sched"]) * 0.7
b=np.array(dataframe["sched"]) * 0.18 + np.array(1 - dataframe["sched"]) * 0.05

Work_modifier = (Low_LR * FR + Med_LR)/(FR + 1)
Schedule_modifier = (a + b * Work_modifier)/(1 + Work_modifier)
Load_multiplier = 1 + Schedule_modifier * (np.array(dataframe["temp"]) - 0) / (30 - 0)
