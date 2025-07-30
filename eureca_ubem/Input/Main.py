##% Import Necessary Modules of Main.py 
import os
from eureca_building.config import load_config
from eureca_ubem.city import City


#%%  Load Inputs 
load_config(os.path.join(".","Example_District_Config.json"))                   #Simulation Settings Given as JSON file 
weather_file = os.path.join(".","ITA_Venezia-Tessera.161050_IGDG.epw")          #Path to weatherfile in epw energyplus format
schedules_file = os.path.join(".","Schedules_total.xlsx")                       #Path to the schedules for the end use
materials_file = os.path.join(".","Materials.xlsx")       #Path to the construction material information
city_model_file = os.path.join(".","Example_District.geojson")                  #Path to the geoindexed file of the footprints of buildings
systems_file = os.path.join(".","systems.xlsx")                                 #Path to the HVAC systems specifications


#%% Generation of the City Object 
city_geojson = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
    systems_templates_file=systems_file,
    shading_calculation=True,                                                   #Shading Calculation Requires Preprocessing Time
    building_model = "2C",                                                      #1C for 5R1C (ISO 13790), 2C for 7R2C (VDI6007)
    output_folder=os.path.join(".","Output_folder")                             
)

#%% Simulations
city_geojson.simulate(output_type="csv")                                        #Comment/Uncomment for Dynamic Simulation   (1C, 2C) 
#city_geojson.simulate_quasi_steady_state()                                      #Comment/Uncomment for Quasi-Steady-State Simulation

