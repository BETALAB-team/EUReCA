''' IMPORTING MODULES '''

import os
import numpy as np
from RC_classes.Simulation import Sim


# Creation of the Sim objectexit
city = Sim()
city_2 = Sim()

# Loading the input data
city.set_input_from_text_file(os.path.join('.','Input','SimInput'))
city_2.set_input_from_excel_file(os.path.join('.','Input','SimInput.xlsx'))

# # Loading weather data, envelopes and schedules
# city.preprocessing()

# # Creation of the district (geometrical processing)
# city.city_creation()

# # Evaluating Urban shadings between buildings
# city.surfaces_and_shading()

# # Calculation buildings parameters
# city.buildings_params_and_loads()

# # Evaluating Urban shadings between buildings
# city.urban_canopy()

# # Design power of buildings and plants creation
# city.plants_design_and_creation()

# # Annual simulation
# city.simulation()

# # Output processing
# city.output()

# %%
