"""
This module includes a method to calculate the cost of scenarios applied to a base scenario of the buildings. 
"""

__author__ = "Mohamad H. Khajedehi"
__credits__ = ["Mohamad H. Khajedehi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Mohamad H. Khajedehi"

import pandas as pd
# from eureca_ubem.city import City

def calculate_intervention_cost(City):

    if not hasattr( City.parent_city , 'base_scenarios' ) or not isinstance ( City.parent_city.base_scenarios, dict ):
        raise AttributeError ("City object has no valid base scenarios, create the scenarios before trying to do the cost anlysis.")
    # if not hasattr(City, 'output_json') or not isinstance(City.output_json, pd.DataFrame):
    #     raise AttributeError ("Scenario object has no valid simulated output, simulate the scenarios before trying to do the cost anlysis.")
               
    
    
    main_building_dict = City.parent_city.buildings_info
    base_buildings = City.buildings_info
    output_df = City.output_geojson.copy()

    
    for b_id, main_b_data in main_building_dict.items():
        
        if b_id not in base_buildings:
            continue
        
        base_b_data = base_buildings [b_id]
        if main_b_data != base_b_data :
            building_object = City.buildings_objects[b_id]

            cost = compute_cost (building_object, City.interventions_applied)
            building_object.intervention_cost = cost
            building_unique_id = main_b_data.get("id")
            
            
            if building_unique_id in output_df ["id"].values:
                output_df.loc[ output_df["id"] == building_unique_id, "intervention_cost [euro]"] = cost
                output_df.loc[ output_df["id"] == building_unique_id, "scenario"] = " ".join(City.interventions_applied)
        
    City.output_geojson = output_df 
    
    
def compute_cost (Building, interventions):
    
    cost = 0 
    
    # area = sum(zone._net_floor_area for zone in Building._Building__thermal_zones_list) #m2
    # volume = sum(zone._volume for zone in Building._Building__thermal_zones_list) #m3
    power_photovoltaic = 0
    if hasattr (Building, "pv_system"):
        power_photovoltaic = Building.pv_system.total_power_installed / 1000 #kW

    area_solarthermal = 0 
    if hasattr (Building, "_heating_system"):
        if hasattr (Building._heating_system, "solar_thermal_system"):
            area_solarthermal=Building._heating_system.solar_thermal_system.installed_area #m2
    design_heating_load = sum(zone.design_heating_system_power for zone in Building._Building__thermal_zones_list)
    design_cooling_load = abs(sum(zone.design_sensible_cooling_system_power for zone in Building._Building__thermal_zones_list))
    opaque_area = sum(zone.ext_wall_opaque_area + zone.ext_roof_area for zone in Building._Building__thermal_zones_list)
    window_area = sum(zone.ext_wall_glazed_area for zone in Building._Building__thermal_zones_list)

    
    if "Envelope" in interventions:
        cost = cost + 180 * opaque_area + 660 * window_area
    if "deep retrofit" in interventions:
        cost = cost + 180 * opaque_area + 660 * window_area + 0.72 * max (design_heating_load,design_cooling_load)
    if "PV" in interventions:
        cost = cost + 1920 * power_photovoltaic
    if "HVAC" in interventions:
        cost = cost + 0.72 * max (design_heating_load,design_cooling_load)
    if "TS" in interventions:
        cost = cost + 1200 * area_solarthermal
    if "PV_TS" in interventions:
        cost = cost + 1200 * area_solarthermal + 2320 * power_photovoltaic

    return cost

    
    
            
        
        
