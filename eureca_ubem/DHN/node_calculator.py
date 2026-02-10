from eureca_ubem.city import City 






def create_dhn_nodes_from_buildings(city: City,
                                    is_reference: bool):
    s=1
    for bd_id, bd_obj in city.buildings_objects:
        x = bd_id 
        b = bd_obj.thermal_zones_list