''' IMPORTING MODULES '''

import os
import numpy as np
from RC_classes.Simulation import Sim


# Creation of the Sim object
city = Sim()

# Loading the input data
city.set_input_from_text_file(os.path.join('.','Input','SimInput'))

# Loading weather data, envelopes and schedules
city.preprocessing()

# Creation of the district (geometrical processing)
city.city_creation()

# Evaluating Urban shadings between buildings
city.urban_shading()

# Calculation buildings parameters
city.buildings_params_and_loads()

# Design power of buildings and plants creation
city.plants_design_and_creation()

# Annual simulation
city.simulation()

# Output processing
city.output()

