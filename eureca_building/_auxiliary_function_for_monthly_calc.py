"""This module includes a function to return monthly average or sums
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging

import pvlib
import numpy as np
import pandas as pd

from eureca_building.config import CONFIG

def get_monthly_value_from_annual_vector(
        array,
        method='sum',
        # start = CONFIG.start_date,
        # end = CONFIG.final_date,
    ):
    """Returns the monthly sum or average of a simulation array

    Parameters
    ----------
    array : numpy.array
        The numpy array to be resample.
    method : str, default 'sum'
        The method can be sum or mean

    Returns
    -------
    numpy.array
        an array with the resampled array
    """
    array_year = pd.Series(array, index = pd.date_range(start="01/01/2023 00:00",end="31/12/2023 23:00",freq=f"{CONFIG.time_step}S"))
    monthly_array = {
        'sum': array_year.resample('1M').sum(),
        'integral': array_year.resample('1M').sum()*CONFIG.time_step, # [W]*[s] = [J]
        'mean': array_year.resample('1M').mean()
    }[method]
    return monthly_array.values
