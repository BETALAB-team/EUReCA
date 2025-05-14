# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 17:45:31 2025

@author: khajmoh18975
"""

import json

# Specify the path to your JSON file
file_path = 'C:/Works/PNRR/8.4.7/Eureca/EUReCA/PostProcess/Supermarkets/3.0_15_MPEurospinIperrossettoViaSarpi.json'

# Open the JSON file and load its content
with open(file_path, 'r') as file:
    data = json.load(file)

# Now, 'data' is a Python object (dict or list, depending on the JSON structure)
print(data)
