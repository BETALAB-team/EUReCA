"""This module includes a container class for schedule end-uses
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

'''IMPORTING MODULES'''

import os
import logging

import pandas as pd
import numpy as np

from eureca_building.config import CONFIG
from eureca_building.systems import hvac_heating_systems_classes, hvac_cooling_systems_classes, HeatingFromParams, CoolingFromParams


# %% ---------------------------------------------------------------------------------------------------
# %% Useful functions to create the schedule EndUses

def load_system_templates(path):
    """
    Loads heating and cooling system templates from an Excel file.
    
    Parameters
    ----------
    path : str
        Path to the spreadsheet containing heating and cooling system definitions.
        
    Returns
    -------
    dict
        Dictionary with two keys: 'heating_systems_templates' and 'cooling_systems_templates',
        each mapping template names to parameter dictionaries.
        
    Notes
    -----
    Populates `hvac_heating_systems_classes` and `hvac_cooling_systems_classes`
    with template keys mapped to parametric constructors.
    """
    # Check input data type

    # read systems templates from excel
    dfs = pd.read_excel(path, sheet_name=None, header=[0], index_col=[0], skiprows=0)

    info_templates = {}

    df = dfs["HeatingSystems"]
    df["name"] = df.index
    info_templates["heating_systems_templates"] = df.to_dict(orient = "index")
    df = dfs["CoolingSystems"]
    df["name"] = df.index
    info_templates["cooling_systems_templates"] = df.to_dict(orient = "index")

    # This adds the templates to the dict of heating/cooling systems classes
    for hs_k in info_templates["heating_systems_templates"].keys():
        hvac_heating_systems_classes[hs_k] = HeatingFromParams
    for cs_k in info_templates["cooling_systems_templates"].keys():
        hvac_cooling_systems_classes[cs_k] = CoolingFromParams

    return info_templates

if __name__ == "__main__":
    path = os.path.join(".","Input","systems.xlsx")
    final_info = load_system_templates(path)
