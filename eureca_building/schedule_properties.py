__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


import numpy as np


"""
This dictionary defines for each schedule type 
the properties using a JSON/Dictionary structure
"""
schedule_types = {
    "unit_type": [
        "dimensionless",
        "percent",
        "temperature",
        "capacity",
        "power",
        "availability",
        "mass_flow_rate"
    ]
}

internal_loads_prop = {
    "people": {
        "unit": ["W", "W/m2", "px", "px/m2", ],
        "fraction_latent": [0., 1.],
        "fraction_radiant": [0., 1.],
        "fraction_convective": [0., 1.],
        "MET": ["W/px"],
        "tags": [],
    },
    "light": {
        "unit": ["W", "W/m2", ],
        "fraction_to_zone": [0., 1.],
        "fraction_radiant": [0., 1.],
        "fraction_convective": [0., 1.],
        "tags": [],
    },
    "electric": {
        "unit": ["W", "W/m2", "W/px"],
        "fraction_to_zone": [0., 1.],
        "fraction_radiant": [0., 1.],
        "fraction_convective": [0., 1.],
        "tags": [],
    },
    "vapour": {
        "unit": ["W", "W/m2", "g", "g/m2", ],
        "fraction_to_zone": [0., 1.],
        "tags": [],
    },
}

setpoint_prop = {
    "temperature": {
        "unit": "°C",
        "limit": [-50., 60.],  # Throws just a warning
        "tags": [],
    },
    "relative_humidity": {
        "unit": "-",
        "limit": [0., 1.],  # Throws just a warning
        "tags": [],
    },
    # Eventually add the specific humidity
}

ventilation_prop = {
    "mechanical": {
        "unit": ["Vol/h", "kg/s", "kg/(m2 s)", "m3/s", "m3/(m2 s)"],
        "temperature limit": [-40., 70.],
        "specific humidity limit": [0.0005, 0.04],
    },
    "infiltration": {
        "unit": ["Vol/h", "kg/s", "kg/(m2 s)", "m3/s", "m3/(m2 s)"],
    },
    "natural": {
        "unit": ["-", "%"],
    },
}

domestic_hot_water_prop = {
    "calculation_method": ["Schedule", "UNI-TS 11300-2", "DHW calc"],
    "unit": ["L/s", "m3/s", "L/(m2 h)"],
    "tags": [],
    "total_days": 365,  # number of days for the simulation
    "nuses" : 4,  # number of uses
    "target temperature [°C]": 45.,

    # Data from DHW calc
    "vol_aver_drawoff" : np.array([15 / 60, 360 / 60, 480 / 60, 480 / 60], dtype=float),  # Liters/min mu
    "vol_desv_drawoff" : np.array([60 / 60, 120 / 60, 120 / 60, 120 / 60], dtype=float),  # Liters/min sigma
    "time_aver_drawoff" : np.array([1, 10, 10, 5], dtype=float),  # minutes
    "proportion_event" : np.array([0.14, 0.36, 0.10, 0.40],
                                dtype=float),  # probability in percentile of 1 for each event to happen from the total daily volume
    "dist" : np.array([[0.01, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.01],
                     [0.01, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.01],
                     [0, 0.0000001, 0.0000002, 0.0000003, 0.0000004, 0.0000005, 0.0103092783505155, 0.0206185567010309,
                      0.0257731958762887, 0.0309278350515464, 0.0463917525773196, 0.0515463917525773, 0.0515463917525773,
                      0.0515463917525773, 0.0515463917525773, 0.0515463917525773, 0.0515463917525773, 0.134020618556701,
                      0.22680412371134, 0.134020618556701, 0.0309278350515464, 0.0206185567010309, 0.0103092783505155, 0],
                     [0, 0.00001, 0.00002, 0.00003, 0.0380952380952381, 0.142857142857143, 0.238095238095238,
                      0.142857142857143, 0.0380952380952381, 0.019047619047619, 0.019047619047619, 0.019047619047619,
                      0.019047619047619, 0.019047619047619, 0.019047619047619, 0.019047619047619, 0.019047619047619,
                      0.0285714285714286, 0.0761904761904762, 0.0761904761904762, 0.0285714285714286, 0.019047619047619,
                      0.019047619047619, 0]], dtype=float).T,
}


