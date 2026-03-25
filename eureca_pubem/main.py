"""
Created on Wed Dec 17 14:23:04 2025 
@author: Mohamad Khajedehi
"""

from scenario.orchestrator import  create_scenario, apply_dhn_changes
from time import time 
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import psutil, os

from eureca_building.config import load_config
from eureca_ubem.city import City 

# S = City ()

# generate scenarios based on the archetype files 
# buildings = {}
# for each scenario 
    
#     S = City.simulate() 
#     for building in s.buildings :
#         building = Building()
#         ==> building heating demand [ scenario]
#         building cooling demand [scenario]
#         building dhw demand 
#         building electricty appliances 
#         building surface area 
#         building shading coefficients 
#         building x and y 
#         buildings[ id ] = building 

# scenario = baselien -> 
# for b in buildings :
#     building.heating_demand = heating demand [baseline]
#     building.cooling_demand = cooling demand [baseline]
#     building.dhn_demand = dhw? + sh?
#     building.elec_consumption = appl + sh? + cooling? + dhw ? 
#     building.elec_generation = surface_area * percentage * sharing * something related to the PV stuff 
#     building.gas_consumption = dhw? + sh?
#     building.stuff_consumption = dhw? + sh? 
#     generate building dfs needed 
    
    
# scenario = new 
# for b in buildings :
#     building.heating_demand = heating demand [baseline]
#     building.cooling_demand = cooling demand [baseline]
#     building.dhn_demand = dhw? + sh?
#     building.elec_consumption = appl + sh? + cooling? + dhw ? 
#     building.elec_generation = surface_area * percentage * sharing * something related to the PV stuff 
#     building.gas_consumption = dhw? + sh?
#     building.stuff_consumption = dhw? + sh? 
#     generate building dfs needed 

    


#%%
df = pd.read_csv ("inputs/dhn_demands.csv", sep = ";")
a = time()
DESIGN = create_scenario (  dhn_demand_path = df,
                            dhn_nodes_path = "inputs/dhn_nodes_2.geojson",
                            elec_consumption_path = "inputs/elec_consumption.csv",
                            elec_production_path = "inputs/elec_production.csv",
                            elec_nodes_path = "inputs/elec_nodes.geojson",
                            climate_file_path = "inputs/SWE_UP_Uppsala.Univ.024620_2015.epw", 
                            streets_path = "inputs/elec_streets.geojson")


#%%
a=time()

BASELINE = create_scenario ( 
        dhn_demand_path = "inputs/dhn_demands.csv",
        climate_file_path = "inputs/SWE_UP_Uppsala.Univ.024620_2015.epw", 
        mode = "BASELINE",
        base = DESIGN
        )


#%%
a=time()

NEW = create_scenario ( 
        dhn_demand_path = "inputs/dhn_demands2.csv",
        elec_consumption_path ="inputs/elec_consumption1.csv",
        elec_production_path = "inputs/elec_production1.csv",
        climate_file_path = "inputs/SWE_UP_Uppsala.Univ.024620_2015.epw", 
        mode = "RETROFIT",
        base = DESIGN
        )





