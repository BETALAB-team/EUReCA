from eureca_ubem.city import City 
from eureca_building.config import CONFIG
import numpy as np
import re

def _thermal_zones_aggregator(building, variable):
    total = None 
    ref_shape = None 
    zones_dict = getattr(building, "_thermal_zones_list")
    for zone in zones_dict:
        arr = getattr(zone, variable)
        if total is None: 
            total = np.zeros_like(arr, dtype=np.result_type(arr, np.float64))
            ref_shape = arr.shape 
        
        total += arr
    return total

def _thermal_zones_picker(building, variable):
    total = None 
    ref_shape = None 
    zones_dict = getattr(building, "_thermal_zones_list")
    for zone in zones_dict:
        arr = getattr(zone, variable)
        if total is None: 
            total = np.zeros_like(arr, dtype=np.result_type(arr, np.float64))
            total += arr 
            return total
def _radiator_temp(Ta, Tmin=23.0, T_switch=12.0):
    def cubic(T):
        return 0.0026*T**3 + 0.019*T**2 - 1.6*T + 56
    y0 = cubic(T_switch)
    m0 = 3*0.0026*T_switch**2 + 2*0.019*T_switch - 1.6  
    def linear(T):
        return y0 + m0*(T - T_switch)
    y = np.where(Ta <= T_switch, cubic(Ta), linear(Ta))
    return np.maximum(Tmin, y)     

    
def reconstruct(A, B, C):
    A = np.asarray(A).ravel()
    B = np.asarray(B).ravel()
    C = np.asarray(C).ravel()

    dB = np.diff(B, prepend=B[0])
    dC = np.diff(C, prepend=C[0])

    def auto_bw(v):
        s = np.median(np.abs(v - np.median(v))) / 0.6745
        return 1.06 * max(s, 1e-12) * (len(v) ** (-1/5))

    bw_B = auto_bw(B)
    bw_dB = auto_bw(dB)

    A_prime = np.empty_like(C)

    for i, (x, dx) in enumerate(zip(C, dC)):
        u1 = (B - x) / bw_B
        u2 = (dB - dx) / bw_dB
        w = np.exp(-0.5 * (u1*u1 + u2*u2))

        if w.sum() < 1e-12:
            A_prime[i] = np.nan
        else:
            A_prime[i] = (w @ A) / w.sum()

    return A_prime

def radiator_temperature_calculation (bd_id, building, space_heating_demand, external_temperature, is_reference, city_ref):
    if (is_reference):
        radiator_temp = _radiator_temp(external_temperature)
        building.ref_radiator_supply_temperature = radiator_temp 
        building.ref_space_heating = space_heating_demand
        
    else: 
        
        reference_sup_temp = city_ref.buildings_objects[bd_id].ref_radiator_supply_temperature
        reference_sh_demand = city_ref.buildings_objects[bd_id].ref_space_heating
        radiator_temp = reconstruct (reference_sup_temp, reference_sh_demand, space_heating_demand)
    return radiator_temp

def create_dhn_nodes_from_buildings(city: City,
                                    reference_city: City):
    
    if city == reference_city:
        is_reference = True
    else:
        is_reference = False
    acceptable_heating_systems = ["District Heating", "IdealLoad"]
    
    print(CONFIG.heating_season_end_time_step)
    dhn_nodes = {}
 
    for bd_id, bd_obj in city.buildings_objects.items():
        building_heating_system = city.buildings_info[bd_id]["Heating System"]
        match = re.match(r"District Heating\s*(\d+)?$", building_heating_system) 
        if match and match.group(1) is not None:
            network_id = int(match.group(1))
            multiple_network_systems = True
        else:
            network_id = 9999999
            multiple_network_systems = False
        if (building_heating_system in acceptable_heating_systems) or multiple_network_systems:
            
            space_heating_load = _thermal_zones_aggregator(bd_obj, "sensible_load_hourly")
            domestic_hot_water_load = _thermal_zones_aggregator(bd_obj, "domestic_hot_water_demand")
            external_temperature = _thermal_zones_picker(bd_obj, "external_temperature")
            radiator_temperature = radiator_temperature_calculation (bd_id, bd_obj, space_heating_load, external_temperature, is_reference, reference_city)
            domestic_hot_water_temperature = 30 
            point = city.buildings_info[bd_id]["geometry"].centroid
            
            dhn_nodes[bd_id] = DHN_node (space_heating_load, domestic_hot_water_load, radiator_temperature,domestic_hot_water_temperature, point, network_id, node_type = "consumer")
            
        
    
    return dhn_nodes
        
        
class DHN_node():
    
    def __init__(self,space_heating_load, domestic_hot_water_load, radiator_temperature, domestic_hot_water_temperature, point, network_id, node_type):
        self.space_heating_load = space_heating_load
        self.domestic_hot_water_load = domestic_hot_water_load 
        self.radiator_supply_temperature = radiator_temperature 
        self.radiator_return_temperature = 0.64 * radiator_temperature + 9.3
        self.domestic_hot_water_supply_temperature = radiator_temperature 
        self.point = point 
        self.network_id = network_id 
        self.node_type = node_type