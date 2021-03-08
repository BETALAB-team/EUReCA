''' IMPORTING MODULES '''

import os
import numpy as np
from RC_classes.Simulation import Sim

exec(open(os.path.join(".","Input","SimInput")).read())

city = Sim()
city.set_input_from_dictionary(Input_files,Sim_input,UWG_data)
city.preprocessing()
city.city_creation()
city.urban_shading()
city.buildings_params_and_loads()
city.plants_design_and_creation()
city.simulation()
city.output()

