'''IMPORTING MODULES'''
import copy
import math
import os
import json
import concurrent.futures
# import psutil

import pickle as pickle
import pandas as pd
import shapely
import geopandas as gpd
import numpy as np
from cjio import cityjson
from scipy.spatial import cKDTree
from shapely.validation import make_valid
from shapely.geometry import MultiPolygon
from eureca_building.config import CONFIG
from eureca_building.pv_system import PV_system
from eureca_building.weather import WeatherFile
from eureca_building.thermal_zone import ThermalZone
from eureca_building.building import Building
from eureca_building.surface import Surface, SurfaceInternalMass
from eureca_building._geometry_auxiliary_functions import normal_versor_2
from eureca_building.air_handling_unit import AirHandlingUnit
from eureca_ubem.end_uses import load_schedules
from eureca_ubem.envelope_types import load_envelopes
from eureca_ubem.systems_templates import load_system_templates
from eureca_ubem.electric_load_italian_distribution import get_italian_random_el_loads

#%% ---------------------------------------------------------------------------------------------------
#%% City class

class City():
    '''This class manages the city simulation via json input file: geojson or cityjson
    '''
    
    # class variables
    T_w_0 = 15.          # Starting average temperature of the walls [°C]
    T_out_inf = 15.      # Starting average temperature outgoing the buildings [°C]
    T_out_AHU = 15.      # Starting average temperature outgoing the Air Handling units [°C]
    V_0_inf = 0.         # Starting average volumetric flow rate outgoing the buildings due to the inflitrations[m3/s]
    V_0_vent = 0.        # Starting average volumetric flow rate outgoing the Air Handling units [m3/s]
    H_waste_0 = 0.       # Starting waste heating rejected by external condensers [kW]

    def __init__(self,
                 city_model:str,
                 envelope_types_file:str,
                 end_uses_types_file:str,
                 epw_weather_file:str,
                 output_folder: str,
                 building_model = "2C",
                 shading_calculation = False,
                 systems_templates_file = None,
                 ):
        """Creates the city from all the input files

        Parameters
        ----------
        city_model : str
            path to the json/cityjson file
        envelope_types_file : str
            path to the envelope_types file
        end_uses_types_file : str
            path to the end_uses file
        epw_weather_file : str
            path to the epw weather file
        output_folder : str
            folder to save results
        building_model : str, default "2C"
            1C or 2C string
        shading_calculation : bool, default False
            whether to do or not the mutual shading calculation
        """

        self.__city_surfaces = []  # List of all the external surfaces of a city
        # Loading weather file
        self.weather_file = WeatherFile(
            epw_weather_file,
            year = CONFIG.simulation_reference_year,
            time_steps = CONFIG.ts_per_hour,
            irradiances_calculation = CONFIG.do_solar_radiation_calculation,
            azimuth_subdivisions = CONFIG.azimuth_subdivisions,
            height_subdivisions = CONFIG.height_subdivisions,
            urban_shading_tol = CONFIG.urban_shading_tolerances
        )

        # Loading Envelope and Schedule Data
        self.envelopes_dict = load_envelopes(envelope_types_file)  # Envelope file loading
        self.end_uses_dict = load_schedules(end_uses_types_file)
        self.systems_templates = load_system_templates(systems_templates_file) if systems_templates_file is not None else None

        self.building_model = building_model
        self.shading_calculation = shading_calculation

        if city_model.endswith('.json'):
            self.buildings_creation_from_cityjson(city_model)
        elif city_model.endswith('.geojson'):
            self.buildings_creation_from_geojson(city_model)
        else:
            raise TypeError(f"City object creation: city model file must be a cityjson or a geojson. City_model: {city_model}")

        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
        self.output_folder = output_folder

    @property
    def building_model(self) -> str:
        return self._building_model

    @building_model.setter
    def building_model(self, value: str):
        if not isinstance(value, str):
            raise TypeError(
                f"City class, building_model is not a string: {value}"
            )
        if value not in ["1C","2C"]:
            # Check if unreasonable values provided
            raise ValueError(
                f"City class, building_model is not a allowed: {value}. please provide an allowed model."
            )
        self._building_model = value

    def buildings_creation_from_cityjson(self,json_path):
        '''This method creates the city from the json file (CityJSON 3D)
        
        Parameters
        ----------
        json_path : str
            path of the cityJSON file

        '''

        # Case of cityJSON file availability:
        self.n_Floors = 0
        # Function to fix invalid geometries
        def fix_geometry(geom):
            if geom.is_valid:
                return geom
            else:
                try:
                    # Try to fix invalid geometry
                    fixed_geom = make_valid(geom)
                    return fixed_geom
                except:
                    return None

        try:
            with open(json_path, 'r') as f:
                self.cityjson = cityjson.CityJSON(file=f)
                #fix geometries
                self.cityjson['geometry'] = self.cityjson['geometry'].apply(fix_geometry)
                self.cityjson = self.cityjson[~self.cityjson['geometry'].isna()]
        except FileNotFoundError:
            raise FileNotFoundError(f'ERROR Cityjson object not found: {json_path}. Give a proper path')

        if not(isinstance(self.cityjson.j['vertices'], list)):
            raise ValueError('json file vertices are not a list')

        self.output_geojson = {}
        self.output_geojson["type"] = "FeatureCollection"
        self.output_geojson["name"] = "Output_geojson"
        self.output_geojson["crs"] = {
            "type": "name",
            "properties":{
                "name": self.cityjson.j["metadata"]["referenceSystem"]
            }
        }

        self.output_geojson["features"] = []


        self.json_buildings= {}
        [(self.json_buildings.update({i:self.cityjson.j['CityObjects'][i]})) for i in self.cityjson.j['CityObjects'] if self.cityjson.j['CityObjects'][i]['type']=='Building']
        self.buildings_objects = {}
        self.buildings_info = {}

        for bd_key, bd_data in self.json_buildings.items():
            # Setting the attributes of the building
            vertices_list = self.cityjson.j['vertices']

            name = bd_key  # Heating plant of the building
            envelope = self.envelopes_dict[bd_data['attributes']['Envelope']]  # Age-class of the building

            surf_counter = 0
            max_height = 0.
            min_height = 1000.
            footprint_area = 0.
            surfaces_list = []
            for geo in bd_data['geometry']:
                if geo['type'] == 'MultiSurface':
                    boundaries = geo['boundaries']
                    for surface in boundaries:
                        for subsurface in surface:
                            surf = []
                            for vert_id in subsurface:
                                surf.append(tuple(vertices_list[vert_id]))

                            surface = Surface(
                                name = f"Bd {name}: surface {surf_counter}",
                                vertices = surf,
                            )

                            if surface.surface_type != "GroundFloor":
                                self.__city_surfaces.append(surface)

                            # TODO: Update wwr calculation

                            if surface.surface_type == "ExtWall":
                                surface._wwr = 0.125

                            surface.construction = {
                                "ExtWall": envelope.external_wall,
                                "Roof": envelope.roof,
                                "GroundFloor": envelope.ground_floor,
                            }[surface.surface_type]

                            surface.window = envelope.window

                            max_height = max(max_height, surface.max_height())
                            min_height = min(min_height, surface.min_height())

                            surfaces_list.append(surface)
                            surf_counter += 1

                            if surface.surface_type == "GroundFloor":
                                footprint_area += surface._area

            # Add internal walls and ceilings 3.3 m height
            floor_height = 3.3
            n_floors = int((max_height - min_height) // floor_height)
            surf_counter = 0
            for i in range(1, n_floors):
                surfaces_list.append(SurfaceInternalMass(
                    name = f"Bd {name}: internal surface {surf_counter}",
                    area = footprint_area,
                    surface_type = "IntCeiling",
                    construction = envelope.interior_ceiling
                ))

                surfaces_list.append(SurfaceInternalMass(
                    name=f"Bd {name}: internal surface {surf_counter + 1}",
                    area=footprint_area,
                    surface_type="IntFloor",
                    construction=envelope.interior_floor
                ))

                surf_counter += 2

            surfaces_list.append(
                SurfaceInternalMass(
                    name=f"Bd {name}: internal surface {surf_counter}",
                    area=footprint_area * n_floors  * 2.5,
                    surface_type="IntWall",
                    construction=envelope.interior_wall
                )
            )

            # Creation of thermal zone and building
            n_units = int(np.around(footprint_area * n_floors / 77.))
            if n_units == 0: n_units = 1

            thermal_zone = ThermalZone(
                name=f"Bd {name} thermal zone",
                surface_list=surfaces_list,
                net_floor_area=footprint_area * n_floors,
                volume=footprint_area * n_floors * floor_height,
                number_of_units=n_units, # 77 average flor area of an appartment according to ISTAT
            )




            self.buildings_info[bd_key] = bd_data['attributes']
            self.buildings_info[bd_key]['Name'] = bd_key
            self.buildings_objects[bd_key] = Building(name=f"Bd {name}", thermal_zones_list=[thermal_zone], model=self.building_model)
            geojson_feature = self.buildings_objects[bd_key].get_geojson_feature_parser()
            geojson_feature["properties"]["Name"] = name
            geojson_feature["properties"]["id"] = name
            geojson_feature["properties"]["new_id"] = name
            geojson_feature["properties"]["End Use"] = self.buildings_info[bd_key]["End Use"]
            geojson_feature["properties"]["Envelope"] = self.buildings_info[bd_key]["Envelope"]
            geojson_feature["properties"]["Heating System"] = self.buildings_info[bd_key]["Heating System"]
            geojson_feature["properties"]["Cooling System"] = self.buildings_info[bd_key]["Cooling System"]
            geojson_feature["properties"]["Height"] = max_height - min_height
            geojson_feature["properties"]["Floors"] = n_floors
            geojson_feature["properties"]["ExtWallCoeff"] = 1.
            geojson_feature["properties"]["VolCoeff"] = 1.

            self.output_geojson["features"].append(geojson_feature)


        self.geometric_preprocessing()

        for bd_k, bd in self.buildings_objects.items():
            thermal_zone = bd._thermal_zones_list[0]
            {
                "1C": thermal_zone._ISO13790_params,
                "2C": thermal_zone._VDI6007_params,
            }[self.building_model]()

        # with open(os.path.join("output_geojson.geojson"), 'w') as outfile:
        #     json.dump(self.output_geojson, outfile)
        self.output_geojson = gpd.read_file(json.dumps(self.output_geojson)).explode(index_parts=True)

    def buildings_creation_from_geojson(self, json_path):
        '''Function to create buildings from geojson file (2D file).

        Parameters
        ----------
        json_path : str
            Path to geojson file.

        '''
        # Function to fix invalid geometries
        def fix_geometry(geom):
            if geom.is_valid:
                return geom
            else:
                try:
                    # Try to fix invalid geometry
                    fixed_geom = make_valid(geom)
                    return fixed_geom
                except:
                    return None
        # Case of GeoJSON file availability:
        
        self.cityjson = gpd.read_file(json_path)# .explode(index_parts=True)
        # for p in self.cityjson.index:
        #     if shapely.geometry.mapping(self.cityjson.loc[p].geometry)["type"] == "Polygon":
        #         self.cityjson["geometry"].loc[p] = MultiPolygon([self.cityjson.loc[p].geometry])
        #fix geometries
        self.cityjson['geometry'] = self.cityjson['geometry'].apply(fix_geometry)
        self.cityjson = self.cityjson[~self.cityjson['geometry'].isna()]
        if "Simulate" not in self.cityjson.columns:
            self.cityjson["Simulate"] = True
        if "ExtWallCoeff" not in self.cityjson.columns:
            self.cityjson["ExtWallCoeff"] = 1.
        if "VolCoeff" not in self.cityjson.columns:
            self.cityjson["VolCoeff"] = 1.
        if "Lower End Use" not in self.cityjson.columns:
            self.cityjson["Lower End Use"] = ''
        if "Upper End Use" not in self.cityjson.columns:
            self.cityjson["Upper End Use"] = ''
        if "Solar technologies" not in self.cityjson.columns:
            self.cityjson["Solar technologies"] = ''
            
            
        if "Windows S Area" not in self.cityjson.columns:
            self.cityjson["Windows S Area"] = np.nan
        if "Windows N Area" not in self.cityjson.columns:
            self.cityjson["Windows N Area"] = np.nan
        if "Windows W Area" not in self.cityjson.columns:
            self.cityjson["Windows W Area"] = np.nan
        if "Windows E Area" not in self.cityjson.columns:
            self.cityjson["Windows E Area"] = np.nan
        
        for direct in ["N","S","E","W"]:
            if f"WWR {direct}" not in self.cityjson.columns:
                self.cityjson[f"WWR {direct}"] = np.nan
            if f"Total wall area {direct}" not in self.cityjson.columns:
                self.cityjson[f"Total wall area {direct}"] = np.nan
        
        self.cityjson["Solar technologies"] = self.cityjson["Solar technologies"].fillna('')
        self.output_geojson = self.cityjson
        self.json_buildings= {}
        self.buildings_objects = {}
        self.buildings_info = {}

        # Extrusion from the footprint operation
        for i in self.cityjson.index:
            id = str(self.cityjson.loc[i]['id']) #  + "_" + str(i[1])
            self.cityjson.loc[i,"new_id"] = id
            self.json_buildings[id] = self.cityjson.loc[i].to_dict()
            bd_data = self.json_buildings[id]
            try:
                n_floors = int(self.cityjson.loc[i]['Floors'])
            except ValueError:
                n_floors = 1
            floor_height = self.cityjson.loc[i]['Height'] / n_floors
            # https://gis.stackexchange.com/questions/287306/list-all-polygon-vertices-coordinates-using-geopandas
            name = str(bd_data["Name"])#  + "_" + str(i[1])
            try:
                building_parts = list(self.cityjson.loc[i].geometry.geoms)
            except AttributeError:
                building_parts = [self.cityjson.loc[i].geometry] # Case for polygon geometry
            surf_counter = 0
            footprint_area = 0.

            building_parts_data = {
                "single End Use": {
                    "surfaces objs": [],
                    "end use": self.cityjson["End Use"],
                    "surfaces counter": 0.,
                    "subpart counter": 0.,
                    "footprint area": 0.,
                },
                
                "upper End Use": {
                    "surfaces objs": [],
                    "end use": self.cityjson["Upper End Use"],
                    "surfaces counter": 0.,
                    "subpart counter": 0.,
                    "footprint area": 0.,
                },
                
                "lower End Use": {
                    "surfaces objs": [],
                    "end use": self.cityjson["Lower End Use"],
                    "surfaces counter": 0.,
                    "subpart counter": 0.,
                    "footprint area": 0.,
                },
            }

            total_areas = {
                "S": 0.,
                "W": 0.,
                "N": 0.,
                "E": 0.,
            }

            for g in building_parts:



                x,y = g.exterior.coords.xy
                coords = np.dstack((x,y)).tolist()
                coords = coords[0]
                coords.pop()
                build_surf = []
                pavimento = []
                soffitto = []
                lower_build_surf=[]
                lower_pavimento = []
                upper_build_surf=[]
                upper_soffitto = []
                z_pav = 0
                z_soff = self.cityjson.loc[i]['Height']
                first_floor_z = floor_height - 0.01 if n_floors == 1 else floor_height
                normal = normal_versor_2(tuple([tuple(coords[n] + [z_pav]) for n in range(len(coords))]))
                if normal[2] > 0.:
                    # Just to adjust in case of anticlockwise perimeter
                    coords.reverse()
                for n in range(len(coords)):
                    pavimento.append(tuple(coords[n]+[z_pav]))
                    soffitto.append(tuple(coords[-n]+[z_soff]))
                    lower_pavimento.append(tuple(coords[n]+[z_pav]))
                    upper_soffitto.append(tuple(coords[-n]+[z_soff]))
                for n in range(len(coords)):
                    build_surf.append(tuple([tuple(coords[n-1]+[z_soff]),\
                                        tuple(coords[n]+[z_soff]),\
                                        tuple(coords[n]+[z_pav]),\
                                        tuple(coords[n-1]+[z_pav])]))
                    lower_build_surf.append(tuple([tuple(coords[n-1]+[first_floor_z]),\
                                        tuple(coords[n]+[first_floor_z]),\
                                        tuple(coords[n]+[z_pav]),\
                                        tuple(coords[n-1]+[z_pav])]))
                    upper_build_surf.append(tuple([tuple(coords[n-1]+[z_soff]),\
                                        tuple(coords[n]+[z_soff]),\
                                        tuple(coords[n]+[z_pav + first_floor_z]),\
                                        tuple(coords[n-1]+[z_pav + first_floor_z])]))

                list_of_int_rings = []
                area_of_int_rings = []
                for int_rings in g.interiors:
                    x,y = int_rings.coords.xy
                    coords_int = np.dstack((x,y)).tolist()[0]
                    coords_int.pop()
                    normal = normal_versor_2(tuple([tuple(coords_int[n] + [z_pav]) for n in range(len(coords_int))]))
                    if normal[2] > 0.:
                        # Just to adjust in case of anticlockwise perimeter
                        coords_int.reverse()
                    list_of_int_rings.append(tuple(coords_int))
                    area_of_int_rings.append(shapely.geometry.Polygon(int_rings).area)

                    # aggiunta delle superfici dei cortili interni nell'edificio (muri verticali)
                    for n in range(len(coords_int)):
                        build_surf.append(tuple([tuple(coords_int[n-1]+[z_soff]),\
                                            tuple(coords_int[n-1]+[z_pav]),\
                                            tuple(coords_int[n]+[z_pav]),\
                                            tuple(coords_int[n]+[z_soff]),]))
                        lower_build_surf.append(tuple([tuple(coords_int[n-1]+[first_floor_z]),\
                                            tuple(coords_int[n-1]+[z_pav]),\
                                            tuple(coords_int[n]+[z_pav]),\
                                            tuple(coords_int[n]+[first_floor_z]),]))
                        upper_build_surf.append(tuple([tuple(coords_int[n-1]+[z_soff]),\
                                            tuple(coords_int[n-1]+[z_pav + first_floor_z]),\
                                            tuple(coords_int[n]+[z_pav + first_floor_z]),\
                                            tuple(coords_int[n]+[z_soff]),]))

                build_surf.append(tuple(pavimento))
                build_surf.append(tuple(soffitto))

                lower_build_surf.append(tuple(lower_pavimento))
                upper_build_surf.append(tuple(upper_soffitto))

                area_of_int_rings = np.array(area_of_int_rings).sum()

                building_parts_data["single End Use"]["surfaces coord"] = build_surf
                building_parts_data["lower End Use"]["surfaces coord"] = lower_build_surf
                building_parts_data["upper End Use"]["surfaces coord"] = upper_build_surf

                # Creation of surfaces
                if bd_data["Simulate"]:
                    envelope = self.envelopes_dict[bd_data['Envelope']]  # Age-class of the building

                maximum_footprint_area = 0.
                for part_k, part_v in building_parts_data.items():
                    for vertices in part_v["surfaces coord"]:
                        surface = Surface(
                                name = f"Bd {name}: surface {surf_counter}",
                                vertices = vertices,
                            )

                        if surface.surface_type != "GroundFloor":
                            self.__city_surfaces.append(surface)

                        if surface.surface_type in ["GroundFloor","Roof"]:
                            surface.reduce_area(area_of_int_rings)

                        # TODO: Update wwr calculation

                        if surface.surface_type == "ExtWall":
                            if surface._azimuth_round == 0:
                                try:
                                    surface._wwr = self.cityjson.loc[i]["WWR S"]
                                except TypeError:
                                    surface._wwr = 0.125
                            elif surface._azimuth_round == 90:
                                try:
                                    surface._wwr = self.cityjson.loc[i]["WWR W"]
                                except TypeError:
                                    surface._wwr = 0.125
                            elif surface._azimuth_round == -90:
                                try:
                                    surface._wwr = self.cityjson.loc[i]["WWR E"]
                                except TypeError:
                                    surface._wwr = 0.125
                            elif surface._azimuth_round == -180:
                                try:
                                    surface._wwr = self.cityjson.loc[i]["WWR N"]
                                except TypeError:
                                    surface._wwr = 0.125

                        if bd_data["Simulate"]:
                            if surface.surface_type in ["GroundFloor","ExtWall", "Roof"]:
                                surface.apply_ext_surf_coeff(bd_data["ExtWallCoeff"])
                            surface.construction = {
                                "ExtWall": envelope.external_wall,
                                "Roof": envelope.roof,
                                "GroundFloor": envelope.ground_floor,
                            }[surface.surface_type]

                            surface.window = envelope.window

                        part_v["surfaces objs"].append(surface)
                        
                        part_v["surfaces counter"] += 1

                        if surface.surface_type == "GroundFloor":
                            part_v["footprint area"] += surface._area
                            
                        part_v["subpart counter"] += 1
                    maximum_footprint_area = np.max([maximum_footprint_area, part_v["footprint area"]])

            for s in building_parts_data["single End Use"]["surfaces objs"]:
                if s.surface_type == "ExtWall":
                    if s._azimuth_round == 0:
                        total_areas["S"] += s._area
                    elif s._azimuth_round == 90:
                        total_areas["W"] += s._area
                    elif s._azimuth_round == -90:
                        total_areas["E"] += s._area
                    elif s._azimuth_round == -180:
                        total_areas["N"] += s._area
                        
            self.output_geojson["Total wall area E"].loc[i] = total_areas["E"]
            self.output_geojson["Total wall area W"].loc[i] = total_areas["W"]
            self.output_geojson["Total wall area N"].loc[i] = total_areas["N"]
            self.output_geojson["Total wall area S"].loc[i] = total_areas["S"]

            # total_wwr = {
            #     "S": 0.,
            #     "N": 0.,
            #     "E": 0.,
            #     "W": 0.,
            # }
            
            # for k, v in total_wwr.items():
            #     try:
            #         total_wwr[k] = bd_data[f'Windows {k} Area'] / total_areas[k]
            #     except ZeroDivisionError:
            #         total_wwr[k] = 0.
            #     if np.isnan(total_wwr[k]):
            #         total_wwr[k] = 0.125


            # for s in building_parts_data["single End Use"]["surfaces objs"]:
            #     if s.surface_type == "ExtWall":
            #         if s._azimuth_round == 0:
            #             s._wwr = total_wwr["S"] if total_wwr["S"] < 0.9 else 0.9
            #             self.output_geojson["WWR S"].loc[i] = s._wwr
            #         elif s._azimuth_round == 90:
            #             s._wwr = total_wwr["W"] if total_wwr["W"] < 0.9 else 0.9
            #             self.output_geojson["WWR W"].loc[i] = s._wwr
            #         elif s._azimuth_round == -90:
            #             s._wwr = total_wwr["E"] if total_wwr["E"] < 0.9 else 0.9
            #             self.output_geojson["WWR E"].loc[i] = s._wwr
            #         elif s._azimuth_round == -180:
            #             s._wwr = total_wwr["N"] if total_wwr["N"] < 0.9 else 0.9
            #             self.output_geojson["WWR N"].loc[i] = s._wwr


            if bd_data["Simulate"]:
                # Add internal walls and ceilings 3.3 m height
                floor_height = self.cityjson.loc[i]['Height'] / n_floors
                surf_counter = 0

                building_parts_data["single End Use"]["number of floors"] = n_floors
                building_parts_data["lower End Use"]["number of floors"] = 1
                building_parts_data["upper End Use"]["number of floors"] = n_floors - 1

                for part_k, part_v in building_parts_data.items():
                    for _ in range(1, part_v["number of floors"]):
                        part_v["surfaces objs"].append(SurfaceInternalMass(
                            name=f"Bd {name}: internal surface {part_v['surfaces counter']}",
                            area=part_v["footprint area"],
                            surface_type="IntCeiling",
                            construction=envelope.interior_ceiling
                        ))

                        part_v["surfaces objs"].append(SurfaceInternalMass(
                            name=f"Bd {name}: internal surface {part_v['surfaces counter'] + 1}",
                            area=part_v["footprint area"],
                            surface_type="IntFloor",
                            construction=envelope.interior_floor
                        ))

                        part_v["surfaces counter"] += 2

                    part_v["surfaces objs"].append(
                        SurfaceInternalMass(
                            name=f"Bd {name}: internal surface {part_v['surfaces counter']}",
                            area=part_v["footprint area"] * part_v["number of floors"] * 2.5,
                            surface_type="IntWall",
                            construction=envelope.interior_wall
                        )
                    )

                    n_units = int(np.around(maximum_footprint_area * part_v["number of floors"] / 77.))
                    if n_units == 0: n_units = 1

                    part_v["tz"] = ThermalZone(
                        name=f"Bd {name} thermal zone {part_k}",
                        surface_list= part_v["surfaces objs"],
                        net_floor_area=maximum_footprint_area * part_v["number of floors"],
                        volume=maximum_footprint_area * part_v["number of floors"] * floor_height * bd_data["VolCoeff"],
                        number_of_units=n_units, # 77 average flor area of an appartment according to ISTAT
                    )

                
                if self.cityjson.loc[i]['Lower End Use'] != '' and (self.cityjson.loc[i]['Lower End Use'] != self.cityjson.loc[i]['Upper End Use']):
                    tz_list = [building_parts_data["upper End Use"]["tz"], building_parts_data["lower End Use"]["tz"]]
                else:
                    tz_list = [building_parts_data["single End Use"]["tz"]]


                
                self.buildings_info[id] = bd_data
                self.buildings_objects[id] = Building(name=f"Bd {name}", thermal_zones_list=tz_list,
                                                              model=self.building_model)

        # Geometric preprocessing
        self.geometric_preprocessing()

        for bd_k, bd in self.buildings_objects.items():
            for thermal_zone in bd._thermal_zones_list:
                {
                    "1C": thermal_zone._ISO13790_params,
                    "2C": thermal_zone._VDI6007_params,
                }[self.building_model]()

    def loads_calculation(self, region = None,  ext_wall_coef = None, t_set = None):
        '''This method does the internal heat gains and solar calculation, as well as it sets the setpoints, ventilation and systems to each building
        '''
        
        if isinstance(region, str):
            italian_el_loads = get_italian_random_el_loads(len(self.buildings_info.values()),region)
            italian_el_loads["Index"] = list(self.buildings_info.keys())
            italian_el_loads.set_index("Index", drop=True, inplace = True)
        
        # iiii=0
        # jjjj=3608
        for bd_id, building_info in self.buildings_info.items():
            # iiii=iiii+1
            building_obj = self.buildings_objects[bd_id]

            if len(building_obj._thermal_zones_list) > 1:
                zones_objs = [
                    [building_obj._thermal_zones_list[0], self.end_uses_dict[building_info["Upper End Use"]]],
                    [building_obj._thermal_zones_list[1], self.end_uses_dict[building_info["Lower End Use"]]],
                ]
            else:
                zones_objs = [
                    [building_obj._thermal_zones_list[0], self.end_uses_dict[building_info["End Use"]]],
                ]

            for tz, use in zones_objs:

                # TODO: copy.deepcopy
                if use.scalar_data["Appliances calculation"] == "Italian Residential Building Stock":
                    try:
                        italian_el_loads
                    except NameError:
                        raise ValueError("""If Italian Residential Building Stock is used as Appliances calculation method, then a region must be passed in this function:
For example: city.loads_calculation(region='Piemonte').
Available regions:
ValleDAosta, Piemonte, Liguria, Lombardia, Veneto, TrentinoAltoAdige,
FriuliVeneziaGiulia, EmiliaRomagna, Umbria, Toscana, Marche, Abruzzo,
Lazio, Campania, Basilicata, Molise, Puglia, Calabria, Sicilia, Sardegna
    """)
                    # TODO: Update with real values.
                    app_nv = italian_el_loads["Tot"].loc[bd_id] * 1000 # W/kW
                    app = copy.deepcopy(use.heat_gains['appliances'])
                    app.unit = "W"
                    app.nominal_value = app_nv / (app.schedule.schedule.sum() / CONFIG.ts_per_hour) * tz.number_of_units

                    tz.add_internal_load(
                        app,
                        use.heat_gains['people'],
                    )
                else:
                    tz.add_internal_load(
                        use.heat_gains['appliances'],
                        use.heat_gains['people'],
                        use.heat_gains['lighting'],
                    )
                
                tz.extract_convective_radiative_latent_electric_load()
                {
                    "1C": tz.calculate_zone_loads_ISO13790,
                    "2C": tz.calculate_zone_loads_VDI6007
                }[self.building_model](self.weather_file)

                tz.add_temperature_setpoint(use.zone_system['temperature_setpoint'])
                tz.add_humidity_setpoint(use.zone_system['humidity_setpoint'])

                tz.add_infiltration(use.infiltration['infiltration'])
                tz.calc_infiltration(self.weather_file)

                ahu = AirHandlingUnit(
                    name = f"ahu Bd {building_info['Name']}",
                    mechanical_vent = use.air_handling_unit_system['ventilation_flow_rate'],
                    supply_temperature = use.air_handling_unit_system['ahu_supply_temperature'],
                    supply_specific_humidity = use.air_handling_unit_system['ahu_supply_humidity'],
                    ahu_operation = use.air_handling_unit_system['ahu_availability'],
                    humidity_control = use.air_handling_unit_system['ahu_humidity_control'],
                    sensible_heat_recovery_eff = use.air_handling_unit_system['ahu_sensible_heat_recovery'],
                    latent_heat_recovery_eff = use.air_handling_unit_system['ahu_latent_heat_recovery'],
                    outdoor_air_ratio = use.air_handling_unit_system['outdoor_air_ratio'],
                    weather = self.weather_file,
                    thermal_zone = tz,
                )
                

                tz.design_sensible_cooling_load(self.weather_file, model=self.building_model)
                
                tz.design_heating_load(-5.)
                
                tz.add_domestic_hot_water(self.weather_file, use.domestic_hot_water['domestic_hot_water'])


            hs_params = None
            cs_params = None

            # Get HVAC paramas for manually set havc systems
            if self.systems_templates is not None:
                if building_info["Heating System"] in self.systems_templates["heating_systems_templates"].keys():
                    hs_params = self.systems_templates["heating_systems_templates"][building_info["Heating System"]]
                if building_info["Cooling System"] in self.systems_templates["cooling_systems_templates"].keys():
                    cs_params = self.systems_templates["cooling_systems_templates"][building_info["Cooling System"]]


            building_obj.set_hvac_system(
                building_info["Heating System"],
                building_info["Cooling System"],
                heating_system_params = hs_params,
                cooling_system_params = cs_params)
            
            building_obj.set_hvac_system_capacity(self.weather_file)
            # print(f"{int(100*iiii/jjjj)} %: load calc")

            if t_set is not None and ext_wall_coef is not None:
                building_obj.update_calibration_params(self.weather_file, ext_wall_coef = ext_wall_coef[bd_id],  h_t_set=t_set[bd_id])
            
            if "PV" in building_info["Solar technologies"]:
                building_obj.add_pv_system(weather_obj=self.weather_file)
                # TODO: add battery/non battery config in solar techologies column ["Only PV, PV and battery"]
            if "ST" in building_info["Solar technologies"]:
                building_obj.add_solar_thermal(weather_obj=self.weather_file)

    def simulate(self, print_single_building_results = False, output_type = "parquet"):
        """Simulation of the whole city, and memorization and stamp of results.

        Parameters
        ----------
        print_single_building_results : bool, default False
            If True, the prints a file with time step results for each building.
            USE CAREFULLY: It might fill a lot of disk space
        output_type : str, default 'parquet'
            It can be either csv (more suitable for excel but more disk space) or parquet (more efficient but not readable from excel)
        """
       
        import time
        start = time.time()
        # parallel simulation commented
        # bd_parallel_list = [[bd, self.weather_file, self.output_folder] for bd in self.buildings_objects.values()]
        # def bd_parallel_solve(x):
        #     bd, weather, out_fold = x
        #     bd.simulate(weather, output_folder=out_fold)
        # with concurrent.futures.ThreadPoolExecutor() as executor:
        #     bd_executor = executor.map(bd_parallel_solve, bd_parallel_list)
        #
        # print(f"Parallel simulation : {(time.time() - start)/60:0.2f} min")
        # start = time.time()

        final_results = {}

        index = pd.date_range(start = CONFIG.start_date,periods = CONFIG.number_of_time_steps, freq = f"{CONFIG.time_step}s")
        district_hourly_results = pd.DataFrame(0., index = index, columns = [
            "Gas consumption [Nm3]",
            "Electric consumption [Wh]",
            "Oil consumption [L]",
            "Wood consumption [kg]",
            "TZ sensible load [W]",
            "TZ latent load [W]",
            "TZ AHU pre heater load [W]",
            "TZ AHU post heater load [W]",
            "TZ DHW demand [W]",
            "Net PV production to grid [Wh]"
        ])
        n_buildings = len(self.buildings_objects)
        counter = 0
        for bd_id, building_info in self.buildings_info.items():
            if counter%10 == 0:
                print(f"{counter} buildings simulated out of {n_buildings}")
            counter += 1

            info = {}
            info["TZ info"] = {}
            for tz in self.buildings_objects[bd_id]._thermal_zones_list:
                info["TZ info"][f"TZ {tz.name}"] = tz.get_zone_info()
            info["TZ info"] = pd.DataFrame(info["TZ info"])
            info["TZ info"]["Total"] = info["TZ info"].sum(axis = 1)
            for i in info["TZ info"].index:
                info[i] = info["TZ info"].loc[i]["Total"]
            info["TZ info"] = info["TZ info"].to_dict()
            info["Name"] = self.buildings_objects[bd_id].name
            if print_single_building_results:
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=self.output_folder, output_type=output_type)
            else:
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=None)
            results.index = index

            demand = results[["TZ sensible load [W]",
                              "TZ latent load [W]",
                              "TZ AHU pre heater load [W]",
                              "TZ AHU post heater load [W]"]]

            dhw_demand = results[["TZ DHW demand [W]"]]

            results["Heating Demand [Wh]"] = demand[demand >= 0].sum(axis = 1) / CONFIG.ts_per_hour
            results["Cooling Demand [Wh]"] = demand[demand < 0].sum(axis = 1) / CONFIG.ts_per_hour
            results["DHW Demand [Wh]"] = dhw_demand[dhw_demand >= 0].sum(axis = 1) / CONFIG.ts_per_hour
            monthly = results.resample("M").sum()

            heat_demand = monthly["Heating Demand [Wh]"]
            cooling_demand = monthly["Cooling Demand [Wh]"]
            dhw_demand = monthly["DHW Demand [Wh]"]
            taken_from_grid = monthly["Taken from the Gird [Wh]"].sum(axis =1)
            given_to_grid = monthly["Given to Grid [Wh]"].sum(axis =1)
            gas_consumption = monthly[[col for col in monthly.columns if "gas consumption" in col[0]]].sum(axis=1)
            el_consumption = monthly[[col for col in monthly.columns if "electric consumption" in col[0]]].sum(axis=1)
            oil_consumption = monthly[[col for col in monthly.columns if "oil consumption" in col[0]]].sum(axis=1)
            wood_consumption = monthly[[col for col in monthly.columns if "wood consumption" in col[0]]].sum(axis=1)
            heat_demand["Total"] = heat_demand.sum()
            cooling_demand["Total"] = cooling_demand.sum()
            dhw_demand["Total"] = dhw_demand.sum()
            gas_consumption["Total"] = gas_consumption.sum()
            el_consumption["Total"] = el_consumption.sum()
            oil_consumption["Total"] = oil_consumption.sum()
            wood_consumption["Total"] = wood_consumption.sum()
            taken_from_grid["Total"] = taken_from_grid.sum()
            given_to_grid["Total"] = given_to_grid.sum()
            for i in gas_consumption.index:
                if i == "Total":
                    info[f"{i} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i} wood consumption [kg]"] = wood_consumption.loc[i]
                    info[f"{i} Heating Demand [Wh]"] = heat_demand.loc[i]
                    info[f"{i} DHW Demand [Wh]"] = dhw_demand.loc[i]
                    info[f"{i} Cooling Demand [Wh]"] = cooling_demand.loc[i]
                    info[f"{i} Taken from grid [Wh]"] = taken_from_grid.loc[i]
                    info[f"{i} Given to grid [Wh]"] = given_to_grid.loc[i]
                else:
                    info[f"{i.month_name()} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i.month_name()} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i.month_name()} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i.month_name()} wood consumption [kg]"] = wood_consumption.loc[i]
                    info[f"{i.month_name()} Heating Demand [Wh]"] = heat_demand.loc[i]
                    info[f"{i.month_name()} DHW Demand [Wh]"] = dhw_demand.loc[i]
                    info[f"{i.month_name()} Cooling Demand [Wh]"] = cooling_demand.loc[i]
                    info[f"{i.month_name()} Taken from grid [Wh]"] = taken_from_grid.loc[i]
                    info[f"{i.month_name()} Given to grid [Wh]"] = given_to_grid.loc[i]
            final_results[bd_id] = info
            district_hourly_results["Gas consumption [Nm3]"] += results["Heating system gas consumption [Nm3]"].iloc[:,0]
            district_hourly_results["Oil consumption [L]"] += results["Heating system oil consumption [L]"].iloc[:,0]
            district_hourly_results["Wood consumption [kg]"] += results["Heating system wood consumption [kg]"].iloc[:,0]
            district_hourly_results["Electric consumption [Wh]"] += results["Heating system electric consumption [Wh]"].iloc[:,0] \
                                                                    + results["Cooling system electric consumption [Wh]"].iloc[:,0] \
                                                                    + results["Appliances electric consumption [Wh]"].iloc[:,0] \
                                                                    + results['AHU electric consumption [Wh]'].iloc[:,0]
            district_hourly_results["TZ sensible load [W]"] += results["TZ sensible load [W]"].iloc[:,0:1].sum(axis = 1)
            district_hourly_results["TZ latent load [W]"] += results["TZ latent load [W]"].iloc[:,0:1].sum(axis = 1)
            district_hourly_results["TZ AHU pre heater load [W]"] += results["TZ AHU pre heater load [W]"].iloc[:,0:1].sum(axis = 1)
            district_hourly_results["TZ AHU post heater load [W]"] += results["TZ AHU post heater load [W]"].iloc[:,0:1].sum(axis = 1)
            district_hourly_results["TZ DHW demand [W]"] += results["TZ DHW demand [W]"].iloc[:,0:1].sum(axis = 1)
            district_hourly_results["Net PV production to grid [Wh]"] += results["Given to Grid [Wh]"].iloc[:,0]


        district_hourly_results.to_csv(os.path.join(self.output_folder,"District_hourly_summary.csv"), sep =";")
        bd_summary = pd.DataFrame.from_dict(final_results,orient="index")
        # bd_summary.to_csv(os.path.join(self.output_folder,"Buildings_summary.csv"), sep =";")
        bd_summary.drop(["Name"], axis = 1, inplace = True)
        self.output_geojson.set_index("new_id", drop=True, inplace = True)
        new_geojson = pd.concat([self.output_geojson.loc[bd_summary.index],bd_summary],axis=1)
        new_geojson.to_file(os.path.join(self.output_folder,"Buildings_summary.geojson"), driver = "GeoJSON")
        new_geojson.drop("geometry",axis = 1).to_csv(os.path.join(self.output_folder, "Buildings_summary.csv"), sep =";")

        print(f"Standard simulation : {(time.time() - start)/60:0.2f} min")

        # def bd_parallel_solve(x):
        #     # Local function to be used in ThreadPoolExecutor
        #     t, bd, weather = x
        #     bd.solve_timestep(t, weather)
        # for t in range(-preprocessing_timesteps, 31*24 - 1):
            # bd_parallel_list = [[t, bd, self.weather_file] for bd in self.buildings_objects.values()]
            # with concurrent.futures.ThreadPoolExecutor() as executor:
            #     bd_executor = executor.map(bd_parallel_solve, bd_parallel_list)

    def geometric_preprocessing(self):
        '''This method firstly reduces the area of coincidence surfaces in the city. This first part must be done to get consistent results
        Moreover, it takes into account the shading effect between buildings surfaces, if shading_calculation is set to True at the city creation
        '''
        
        toll_az = CONFIG.urban_shading_tolerances[0]
        toll_dist = CONFIG.urban_shading_tolerances[1]
        toll_theta = CONFIG.urban_shading_tolerances[2]

        # Each surface is compared with all the others
        # qqq=len(self.__city_surfaces)
        # qqqq=qqq*(qqq-1)/2
        # iiii=0
        # jjjj=0
        list_of_centroids = [obj._centroid for obj in self.__city_surfaces]
        plg_kdtree=cKDTree(list_of_centroids)
        max_number_of_neighborhoods = len(self.__city_surfaces) if len(self.__city_surfaces) < 500 else 500
        for x in range(len(self.__city_surfaces)):
            # jjjj=jjjj+1
            # print(f"{int(10000*jjjj/len(self.__city_surfaces))/100} percent: geometry")
           
            # Function to filter using scipy the nearest points (centroids)
            _, filtered_indices = plg_kdtree.query(list_of_centroids[x], k=max_number_of_neighborhoods)
            
            # iiiii=0
            try:
                self.__city_surfaces[x].shading_coupled_surfaces
            except AttributeError:
                self.__city_surfaces[x].shading_coupled_surfaces = []
            for y in filtered_indices:
                # iiii=iiii+1
                # iiiii=iiiii+1
                # print(f" {iiii}:{jjjj}//{len(self.__city_surfaces)},{iiiii}:{len(check_indices)}")

                try:
                    self.__city_surfaces[y].shading_coupled_surfaces
                except AttributeError:
                    self.__city_surfaces[y].shading_coupled_surfaces = []

                # Calculation of the distance between the centroids of the two surfaces under examination
                dist = math.sqrt(
                    (self.__city_surfaces[x]._centroid[0]-self.__city_surfaces[y]._centroid[0])**2 +
                    (self.__city_surfaces[x]._centroid[1]-self.__city_surfaces[y]._centroid[1])**2
                )

                # Reducing the area of the coincidence surfaces
                if dist < 15.:
                    if self.__city_surfaces[x].check_surface_coincidence(self.__city_surfaces[y]):
                        intersectionArea = self.__city_surfaces[x].calculate_intersection_area(self.__city_surfaces[y])
                        self.__city_surfaces[y].reduce_area(intersectionArea)
                        self.__city_surfaces[x].reduce_area(intersectionArea)

                # print(f" {int(iiii/jjjj)*100} percent: geometry")
                try:
                    self.__city_surfaces[y].shading_coupled_surfaces
                except AttributeError:
                    self.__city_surfaces[y].shading_coupled_surfaces = []

                if self.shading_calculation:
                    if dist == 0.0:
                        theta=0
                    else:

                        # Calculation of the vector direction between the centroids of the two surfaces under examination
                        theta_xy = np.degrees(np.arccos((self.__city_surfaces[y]._centroid[0]-self.__city_surfaces[x]._centroid[0])/dist))
                        theta = -(theta_xy + 90)
                        if self.__city_surfaces[y]._centroid[1] < self.__city_surfaces[x]._centroid[1]:
                            theta = theta + 2*theta_xy
                        if theta < -180:
                            theta = theta + 360
                        if theta > 180:
                            theta = theta - 360

                    # Conditions:
                    #    1. the distance between surfaces must be less than toll_dist
                    #    2. the theta angle between surfaces must be within the range
                    #    3. the azimuth angle of the second surface must be within the range

                    if dist < toll_dist:
                        if self.__city_surfaces[x]._azimuth < 0:
                            azimuth_opp = self.__city_surfaces[x]._azimuth + 180
                            azimuth_opp_max = azimuth_opp + toll_az
                            azimuth_opp_min = azimuth_opp - toll_az
                            theta_max = self.__city_surfaces[x]._azimuth + toll_theta
                            theta_min = self.__city_surfaces[x]._azimuth - toll_theta
                            if theta_min < -180:
                                theta_min = theta_min + 360
                                if theta_min < theta < 180 or -180 < theta < theta_max:
                                    if azimuth_opp_max > 180:
                                        azimuth_opp_max = azimuth_opp_max - 360
                                        if azimuth_opp_min < self.__city_surfaces[y]._azimuth <= 180 or -180 <= self.__city_surfaces[y]._azimuth < azimuth_opp_max:
                                            self.__city_surfaces[x].shading_coupled_surfaces.append([dist,y])
                                            self.__city_surfaces[y].shading_coupled_surfaces.append([dist,x])
                                    else:
                                        if azimuth_opp_min < self.__city_surfaces[y]._azimuth < azimuth_opp_max:
                                            self.__city_surfaces[x].shading_coupled_surfaces.append([dist,y])
                                            self.__city_surfaces[y].shading_coupled_surfaces.append([dist,x])
                            else:
                                if theta_min < theta < theta_max:
                                    if azimuth_opp_max > 180:
                                        azimuth_opp_max = azimuth_opp_max - 360
                                        if azimuth_opp_min < self.__city_surfaces[y]._azimuth <= 180 or -180 <= self.__city_surfaces[y]._azimuth < azimuth_opp_max:
                                            self.__city_surfaces[x].shading_coupled_surfaces.append([dist,y])
                                            self.__city_surfaces[y].shading_coupled_surfaces.append([dist,x])
                                    else:
                                        if azimuth_opp_min < self.__city_surfaces[y]._azimuth < azimuth_opp_max:
                                            self.__city_surfaces[x].shading_coupled_surfaces.append([dist,y])
                                            self.__city_surfaces[y].shading_coupled_surfaces.append([dist,x])
                        else:
                            azimuth_opp = self.__city_surfaces[x]._azimuth - 180
                            azimuth_opp_max = azimuth_opp + toll_az
                            azimuth_opp_min = azimuth_opp - toll_az
                            theta_max = self.__city_surfaces[x]._azimuth + toll_theta
                            theta_min = self.__city_surfaces[x]._azimuth - toll_theta
                            if theta_max > 180:
                                theta_max = theta_max - 360
                                if theta_min < theta < 180 or -180 <= theta < theta_max:
                                    if azimuth_opp_min < -180:
                                        azimuth_opp_min = azimuth_opp_min + 360
                                        if -180 <= self.__city_surfaces[y]._azimuth < azimuth_opp_max or azimuth_opp_min < self.__city_surfaces[y]._azimuth <= 180:
                                            self.__city_surfaces[x].shading_coupled_surfaces.append([dist,y])
                                            self.__city_surfaces[y].shading_coupled_surfaces.append([dist,x])
                                    else:
                                        if azimuth_opp_min < self.__city_surfaces[y]._azimuth < azimuth_opp_max:
                                            self.__city_surfaces[x].shading_coupled_surfaces.append([dist,y])
                                            self.__city_surfaces[y].shading_coupled_surfaces.append([dist,x])
                            else:
                                if theta_min < theta < theta_max:
                                    if azimuth_opp_min < -180:
                                        azimuth_opp_min = azimuth_opp_min + 360
                                        if -180 <= self.__city_surfaces[y]._azimuth < azimuth_opp_max or azimuth_opp_min < self.__city_surfaces[y]._azimuth <= 180:
                                            self.__city_surfaces[x].shading_coupled_surfaces.append([dist,y])
                                            self.__city_surfaces[y].shading_coupled_surfaces.append([dist,x])
                                    else:
                                        if azimuth_opp_min < self.__city_surfaces[y]._azimuth < azimuth_opp_max:
                                            self.__city_surfaces[x].shading_coupled_surfaces.append([dist,y])
                                            self.__city_surfaces[y].shading_coupled_surfaces.append([dist,x])
        
        if self.shading_calculation:
            # SECTION 2: Calculation of the shading effect
            for x in range(len(self.__city_surfaces)):
                self.__city_surfaces[x].shading_coefficient = 1.
                # iiiii=0
                if self.__city_surfaces[x].shading_coupled_surfaces != []:
                    self.__city_surfaces[x].shading_calculation = 'On'
                    shading = [0]*len(self.__city_surfaces[x].shading_coupled_surfaces)
                    for y in range(len(self.__city_surfaces[x].shading_coupled_surfaces)):
                        # iiii=iiii+1
                        # iiiii=iiiii+1
                        # print(f" {iiii}:{jjjj}//{len(self.__city_surfaces)},{iiiii}:{len(self.__city_surfaces[x].shading_coupled_surfaces)}")

                        distance = self.__city_surfaces[x].shading_coupled_surfaces[y][0]
                        surface_opposite_index = self.__city_surfaces[x].shading_coupled_surfaces[y][1]

                        # Calculation of the solar height limit
                        if distance == 0:
                            # Case of distance = 0
                            sol_h_lim = 90.
                        else:
                            sol_h_lim = np.degrees(
                                np.arctan(
                                    (self.__city_surfaces[surface_opposite_index].max_height() -
                                     self.__city_surfaces[x]._centroid[2])/distance
                                )
                            )
                        self.__city_surfaces[x].shading_coupled_surfaces[y].append(sol_h_lim)
                        
                        # Calculation of the solar azimuth limits
                        sol_az_lim1_xy = np.degrees(
                            np.arccos(
                                (self.__city_surfaces[surface_opposite_index]._vertices[0][0] -
                                 self.__city_surfaces[x]._centroid[0])/
                                math.sqrt(
                                    (self.__city_surfaces[x]._centroid[0] -
                                     self.__city_surfaces[surface_opposite_index]._vertices[0][0])**2
                                    + (self.__city_surfaces[x]._centroid[1] -
                                       self.__city_surfaces[surface_opposite_index]._vertices[0][1])**2)
                            )
                        )
                        sol_az_lim2_xy = np.degrees(
                            np.arccos(
                                (self.__city_surfaces[surface_opposite_index]._vertices[2][0]
                                 - self.__city_surfaces[x]._centroid[0])
                                /math.sqrt(
                                    (self.__city_surfaces[x]._centroid[0] -
                                     self.__city_surfaces[surface_opposite_index]._vertices[2][0])**2 +
                                    (self.__city_surfaces[x]._centroid[1]
                                     - self.__city_surfaces[surface_opposite_index]._vertices[2][1])**2
                                )
                            )
                        )
                        sol_az_lim1 = -(sol_az_lim1_xy + 90)
                        sol_az_lim2 = -(sol_az_lim2_xy + 90)
                        if self.__city_surfaces[surface_opposite_index]._vertices[0][1] < self.__city_surfaces[x]._centroid[1]:
                            sol_az_lim1 = sol_az_lim1 + 2*sol_az_lim1_xy
                        if sol_az_lim1 < -180:
                            sol_az_lim1 = sol_az_lim1 + 360
                        if sol_az_lim1 > 180:
                            sol_az_lim1 = sol_az_lim1 - 360
                        if self.__city_surfaces[surface_opposite_index]._vertices[2][1] < self.__city_surfaces[x]._centroid[1]:
                            sol_az_lim2 = sol_az_lim2 + 2*sol_az_lim2_xy
                        if sol_az_lim2 < -180:
                            sol_az_lim2 = sol_az_lim2 + 360
                        if sol_az_lim2 > 180:
                            sol_az_lim2 = sol_az_lim2 - 360
                        
                        # Necessary conditions:
                        #    1. solar height less than the solar height limit
                        #    2. solar azimuth between the solar azimuth limits
                        
                        shading_sol_h = np.less(
                            self.weather_file.hourly_data["solar_position_elevation"],
                            sol_h_lim
                        )
                        sol_az_lim_inf = min(sol_az_lim1,sol_az_lim2)
                        sol_az_lim_sup = max(sol_az_lim1,sol_az_lim2)
                        if abs(sol_az_lim_inf - sol_az_lim_sup) < 180:
                            self.__city_surfaces[x].shading_coupled_surfaces[y].append([sol_az_lim_inf,sol_az_lim_sup])
                            shading_sol_az = [
                                np.less(
                                self.weather_file.hourly_data["solar_position_azimuth"],sol_az_lim_sup
                            ),np.greater(
                                self.weather_file.hourly_data["solar_position_azimuth"],sol_az_lim_inf
                            )
                            ]
                            shading_tot = shading_sol_h & shading_sol_az[0] & shading_sol_az[1]
                        else:
                            sol_az_lim_inf = max(sol_az_lim1,sol_az_lim2)
                            sol_az_lim_sup = min(sol_az_lim1,sol_az_lim2)
                            self.__city_surfaces[x].shading_coupled_surfaces[y].append([sol_az_lim_inf,sol_az_lim_sup])
                            shading_sol_az1 = [np.less_equal(self.weather_file.hourly_data["solar_position_azimuth"],180),
                                               np.greater(self.weather_file.hourly_data["solar_position_azimuth"],sol_az_lim_inf)]
                            shading_sol_az2 = [np.less(self.weather_file.hourly_data["solar_position_azimuth"],sol_az_lim_sup),
                                               np.greater_equal(self.weather_file.hourly_data["solar_position_azimuth"],-180)]
                            shading_az1 = shading_sol_az1[0] & shading_sol_az1[1]
                            shading_az2 = shading_sol_az2[0] & shading_sol_az2[1]
                            shading_az = shading_az1 | shading_az2
                            shading_tot = shading_sol_h & shading_az
                        shading[y] = shading_tot
                    for y in range(len(self.__city_surfaces[x].shading_coupled_surfaces)):
                        if y == 0:
                            shading_eff = shading[y]
                        else:
                            shading_eff = shading_eff | shading[y]
                    shading_eff_01 = (1 - np.where(shading_eff==True,1,shading_eff))
                    shading_eff_01 = np.where(shading_eff_01==0, 0. ,shading_eff_01)
                    self.__city_surfaces[x].shading_coefficient = shading_eff_01
                    
  
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
  
       
                    

    def calibrate_electricity_winter(self, print_single_building_results=False, output_type="parquet",
                  limits_ele_load = [0.2,1.8],
                  ):
        """Simulation of the whole city, and memorization and stamp of results.
    
        Parameters
        ----------
        print_single_building_results : bool, default False
            If True, the prints a file with time step results for each building.
            USE CAREFULLY: It might fill a lot of disk space
        output_type : str, default 'parquet'
            It can be either csv (more suitable for excel but more disk space) or parquet (more efficient but not readable from excel)
        """
    
        import time
        start = time.time()
        # parallel simulation commented
        # bd_parallel_list = [[bd, self.weather_file, self.output_folder] for bd in self.buildings_objects.values()]
        # def bd_parallel_solve(x):
        #     bd, weather, out_fold = x
        #     bd.simulate(weather, output_folder=out_fold)
        # with concurrent.futures.ThreadPoolExecutor() as executor:
        #     bd_executor = executor.map(bd_parallel_solve, bd_parallel_list)
        #
        # print(f"Parallel simulation : {(time.time() - start)/60:0.2f} min")
        # start = time.time()
    
        final_results = {}
    
        index = pd.date_range(start=CONFIG.start_date, periods=CONFIG.number_of_time_steps, freq=f"{CONFIG.time_step}s")
        district_hourly_results = pd.DataFrame(0., index=index, columns=[
            "Gas consumption [Nm3]",
            "Electric consumption [Wh]",
            "Oil consumption [L]",
            "Wood consumption [kg]",
            "TZ sensible load [W]",
            "TZ latent load [W]",
            "TZ AHU pre heater load [W]",
            "TZ AHU post heater load [W]",
            "TZ DHW demand [W]",
        ])
        n_buildings = len(self.buildings_objects)
        counter = 0
        for bd_id, building_info in self.buildings_info.items():
            if counter % 10 == 0:
                print(f"{counter} buildings simulated out of {n_buildings}")
            counter += 1
            
            self.buildings_objects[bd_id].info_geojson = building_info
    
            info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
            info["Name"] = self.buildings_objects[bd_id].name
            if print_single_building_results:
                building_info
                
                if not np.isnan(building_info["Ele Non Cooling [kWh]"]):
                    cal_res = self.buildings_objects[bd_id].calibrate_el_load_wint(building_info["Ele Non Cooling [kWh]"],
                                                                      self.weather_file,
                                                                      limits_ele_loads=limits_ele_load,
                                                                      )
                else:
                    cal_res = [[1.], None, None, None, None, None]
                    
                self.buildings_objects[bd_id].update_calibration_params(self.weather_file,
                                                                        elec_loads_m = cal_res[0][0]) 
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=self.output_folder,
                                                                  output_type=output_type)
                info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
                info["Name"] = self.buildings_objects[bd_id].name
                info["ElLoad_cal"] = cal_res[0][0]
                info["Calibration output"] = cal_res[4]
                info["Calibration iterations"] = cal_res[3]
                info["Calibration status"] = cal_res[5]
            else:
                if not np.isnan(building_info["Ele Non Cooling [kWh]"]):
                    cal_res = self.buildings_objects[bd_id].calibrate_el_load_wint(building_info["Ele Non Cooling [kWh]"], self.weather_file,
                                                                  limits_ele_loads=limits_ele_load,
                                                                )
                else:
                    cal_res = [[1.], None, None, None, None, None]
                self.buildings_objects[bd_id].update_calibration_params(self.weather_file,
                                                                       elec_loads_m = cal_res[0][0]) 
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=None)
                info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
                info["Name"] = self.buildings_objects[bd_id].name
                info["ElLoad_cal"] = cal_res[0][0]
                info["Calibration output"] = cal_res[4]
                info["Calibration iterations"] = cal_res[3]
                info["Calibration status"] = cal_res[5]
            results.index = index
            monthly = results.resample("M").sum()
            gas_consumption = monthly[[col for col in monthly.columns if "gas consumption" in col[0]]].sum(axis=1)
            el_consumption = monthly[[col for col in monthly.columns if "electric consumption" in col[0]]].sum(axis=1)
            oil_consumption = monthly[[col for col in monthly.columns if "oil consumption" in col[0]]].sum(axis=1)
            wood_consumption = monthly[[col for col in monthly.columns if "wood consumption" in col[0]]].sum(axis=1)
            gas_consumption["Total"] = gas_consumption.sum()
            el_consumption["Total"] = el_consumption.sum()
            oil_consumption["Total"] = oil_consumption.sum()
            wood_consumption["Total"] = wood_consumption.sum()
            for i in gas_consumption.index:
                if i == "Total":
                    info[f"{i} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i} wood consumption [kg]"] = wood_consumption.loc[i]
                else:
                    info[f"{i.month_name()} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i.month_name()} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i.month_name()} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i.month_name()} wood consumption [kg]"] = wood_consumption.loc[i]
            final_results[bd_id] = info
            district_hourly_results["Gas consumption [Nm3]"] += results["Heating system gas consumption [Nm3]"].iloc[:,
                                                                0]
            district_hourly_results["Oil consumption [L]"] += results["Heating system oil consumption [L]"].iloc[:, 0]
            district_hourly_results["Wood consumption [kg]"] += results["Heating system wood consumption [kg]"].iloc[:,
                                                                0]
            district_hourly_results["Electric consumption [Wh]"] += results[
                                                                        "Heating system electric consumption [Wh]"].iloc[
                                                                    :, 0] \
                                                                    + results[
                                                                          "Cooling system electric consumption [Wh]"].iloc[
                                                                      :, 0] \
                                                                    + results[
                                                                          "Appliances electric consumption [Wh]"].iloc[
                                                                      :, 0]
            district_hourly_results["TZ sensible load [W]"] += results["TZ sensible load [W]"].iloc[:, 0]
            district_hourly_results["TZ latent load [W]"] += results["TZ latent load [W]"].iloc[:, 0]
            district_hourly_results["TZ AHU pre heater load [W]"] += results["TZ AHU pre heater load [W]"].iloc[:, 0]
            district_hourly_results["TZ AHU post heater load [W]"] += results["TZ AHU post heater load [W]"].iloc[:, 0]
            district_hourly_results["TZ DHW demand [W]"] += results["TZ DHW demand [W]"].iloc[:, 0]
    
        district_hourly_results.to_csv(os.path.join(self.output_folder, "District_hourly_summary.csv"), sep=";")
        bd_summary = pd.DataFrame.from_dict(final_results, orient="index")
        # bd_summary.to_csv(os.path.join(self.output_folder,"Buildings_summary.csv"), sep =";")
        bd_summary.drop(["Name"], axis=1, inplace=True)
        self.output_geojson.set_index("new_id", drop=True, inplace=True)
        new_geojson = pd.concat([self.output_geojson.loc[bd_summary.index],bd_summary],axis=1)
        new_geojson.to_file(os.path.join(self.output_folder, "Buildings_summary.geojson"), driver="GeoJSON")
        new_geojson.drop("geometry", axis=1).to_csv(os.path.join(self.output_folder, "Buildings_summary.csv"), sep=";")
    
        print(f"Standard simulation : {(time.time() - start) / 60:0.2f} min")

    

    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
  
    # def simulate_electricity_summer(self, print_single_building_results=False, output_type="parquet",
    #               c_setpoint_column = "TSetC_cal",
    #               mfr_m_column = "MFR_M_cal"
    #               ):
    #     """Simulation of the whole city, and memorization and stamp of results.
    
    #     Parameters
    #     ----------
    #     print_single_building_results : bool, default False
    #         If True, the prints a file with time step results for each building.
    #         USE CAREFULLY: It might fill a lot of disk space
    #     output_type : str, default 'parquet'
    #         It can be either csv (more suitable for excel but more disk space) or parquet (more efficient but not readable from excel)
    #     """
    
    #     import time
    #     start = time.time()
    #     # parallel simulation commented
    #     # bd_parallel_list = [[bd, self.weather_file, self.output_folder] for bd in self.buildings_objects.values()]
    #     # def bd_parallel_solve(x):
    #     #     bd, weather, out_fold = x
    #     #     bd.simulate(weather, output_folder=out_fold)
    #     # with concurrent.futures.ThreadPoolExecutor() as executor:
    #     #     bd_executor = executor.map(bd_parallel_solve, bd_parallel_list)
    #     #
    #     # print(f"Parallel simulation : {(time.time() - start)/60:0.2f} min")
    #     # start = time.time()
    
    #     final_results = {}
    
    #     index = pd.date_range(start=CONFIG.start_date, periods=CONFIG.number_of_time_steps, freq=f"{CONFIG.time_step}s")
    #     district_hourly_results = pd.DataFrame(0., index=index, columns=[
    #         "Gas consumption [Nm3]",
    #         "Electric consumption [Wh]",
    #         "Oil consumption [L]",
    #         "Wood consumption [kg]",
    #         "TZ sensible load [W]",
    #         "TZ latent load [W]",
    #         "TZ AHU pre heater load [W]",
    #         "TZ AHU post heater load [W]",
    #         "TZ DHW demand [W]",
    #     ])
    #     n_buildings = len(self.buildings_objects)
    #     counter = 0
    #     for bd_id, building_info in self.buildings_info.items():
    #         if counter % 10 == 0:
    #             print(f"{counter} buildings simulated out of {n_buildings}")
    #         counter += 1
    
    #         info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
    #         info["Name"] = self.buildings_objects[bd_id].name
    #         if print_single_building_results:
    #             self.buildings_objects[bd_id].update_calibration_params(self.weather_file,
    #                                                                     mfr_m = building_info["MFR_M_cal"], 
    #                                                                     c_t_set = building_info["TSetC_cal"]) 
    #             results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=self.output_folder,
    #                                                               output_type=output_type)
    #             info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
    #             info["Name"] = self.buildings_objects[bd_id].name

    #         else:
    #             self.buildings_objects[bd_id].update_calibration_params(self.weather_file,
    #                                                                     mfr_m = building_info["MFR_M_cal"], 
    #                                                                     c_t_set = building_info["TSetC_cal"]) 
    #             results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=None)
    #             info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
    #             info["Name"] = self.buildings_objects[bd_id].name
    #         results.index = index
    #         monthly = results.resample("M").sum()
    #         gas_consumption = monthly[[col for col in monthly.columns if "gas consumption" in col[0]]].sum(axis=1)
    #         el_consumption = monthly[[col for col in monthly.columns if "electric consumption" in col[0]]].sum(axis=1)
    #         oil_consumption = monthly[[col for col in monthly.columns if "oil consumption" in col[0]]].sum(axis=1)
    #         wood_consumption = monthly[[col for col in monthly.columns if "wood consumption" in col[0]]].sum(axis=1)
    #         gas_consumption["Total"] = gas_consumption.sum()
    #         el_consumption["Total"] = el_consumption.sum()
    #         oil_consumption["Total"] = oil_consumption.sum()
    #         wood_consumption["Total"] = wood_consumption.sum()
    #         for i in gas_consumption.index:
    #             if i == "Total":
    #                 info[f"{i} gas consumption [Nm3]"] = gas_consumption.loc[i]
    #                 info[f"{i} electric consumption [Wh]"] = el_consumption.loc[i]
    #                 info[f"{i} oil consumption [L]"] = oil_consumption.loc[i]
    #                 info[f"{i} wood consumption [kg]"] = wood_consumption.loc[i]
    #             else:
    #                 info[f"{i.month_name()} gas consumption [Nm3]"] = gas_consumption.loc[i]
    #                 info[f"{i.month_name()} electric consumption [Wh]"] = el_consumption.loc[i]
    #                 info[f"{i.month_name()} oil consumption [L]"] = oil_consumption.loc[i]
    #                 info[f"{i.month_name()} wood consumption [kg]"] = wood_consumption.loc[i]
    #         final_results[bd_id] = info
    #         district_hourly_results["Gas consumption [Nm3]"] += results["Heating system gas consumption [Nm3]"].iloc[:,
    #                                                             0]
    #         district_hourly_results["Oil consumption [L]"] += results["Heating system oil consumption [L]"].iloc[:, 0]
    #         district_hourly_results["Wood consumption [kg]"] += results["Heating system wood consumption [kg]"].iloc[:,
    #                                                             0]
    #         district_hourly_results["Electric consumption [Wh]"] += results[
    #                                                                     "Heating system electric consumption [Wh]"].iloc[
    #                                                                 :, 0] \
    #                                                                 + results[
    #                                                                       "Cooling system electric consumption [Wh]"].iloc[
    #                                                                   :, 0] \
    #                                                                 + results[
    #                                                                       "Appliances electric consumption [Wh]"].iloc[
    #                                                                   :, 0]
    #         district_hourly_results["TZ sensible load [W]"] += results["TZ sensible load [W]"].iloc[:, 0]
    #         district_hourly_results["TZ latent load [W]"] += results["TZ latent load [W]"].iloc[:, 0]
    #         district_hourly_results["TZ AHU pre heater load [W]"] += results["TZ AHU pre heater load [W]"].iloc[:, 0]
    #         district_hourly_results["TZ AHU post heater load [W]"] += results["TZ AHU post heater load [W]"].iloc[:, 0]
    #         district_hourly_results["TZ DHW demand [W]"] += results["TZ DHW demand [W]"].iloc[:, 0]
    
    #     district_hourly_results.to_csv(os.path.join(self.output_folder, "District_hourly_summary.csv"), sep=";")
    #     bd_summary = pd.DataFrame.from_dict(final_results, orient="index")
    #     # bd_summary.to_csv(os.path.join(self.output_folder,"Buildings_summary.csv"), sep =";")
    #     bd_summary.drop(["Name"], axis=1, inplace=True)
    #     self.output_geojson.set_index("new_id", drop=True, inplace=True)
    #     new_geojson = pd.concat([self.output_geojson.loc[bd_summary.index],bd_summary],axis=1)
    #     new_geojson.to_file(os.path.join(self.output_folder, "Buildings_summary.geojson"), driver="GeoJSON")
    #     new_geojson.drop("geometry", axis=1).to_csv(os.path.join(self.output_folder, "Buildings_summary.csv"), sep=";")
    
    #     print(f"Standard simulation : {(time.time() - start) / 60:0.2f} min")             
                    

    def calibrate_electricity_summer(self, print_single_building_results=False, output_type="parquet",
                  limits_setpoint = [24,28],
                  limits_mfr_m = [0.3,1.]
                  ):
        """Simulation of the whole city, and memorization and stamp of results.
    
        Parameters
        ----------
        print_single_building_results : bool, default False
            If True, the prints a file with time step results for each building.
            USE CAREFULLY: It might fill a lot of disk space
        output_type : str, default 'parquet'
            It can be either csv (more suitable for excel but more disk space) or parquet (more efficient but not readable from excel)
        """
    
        import time
        start = time.time()
        # parallel simulation commented
        # bd_parallel_list = [[bd, self.weather_file, self.output_folder] for bd in self.buildings_objects.values()]
        # def bd_parallel_solve(x):
        #     bd, weather, out_fold = x
        #     bd.simulate(weather, output_folder=out_fold)
        # with concurrent.futures.ThreadPoolExecutor() as executor:
        #     bd_executor = executor.map(bd_parallel_solve, bd_parallel_list)
        #
        # print(f"Parallel simulation : {(time.time() - start)/60:0.2f} min")
        # start = time.time()
    
        final_results = {}
    
        index = pd.date_range(start=CONFIG.start_date, periods=CONFIG.number_of_time_steps, freq=f"{CONFIG.time_step}s")
        district_hourly_results = pd.DataFrame(0., index=index, columns=[
            "Gas consumption [Nm3]",
            "Electric consumption [Wh]",
            "Oil consumption [L]",
            "Wood consumption [kg]",
            "TZ sensible load [W]",
            "TZ latent load [W]",
            "TZ AHU pre heater load [W]",
            "TZ AHU post heater load [W]",
            "TZ DHW demand [W]",
        ])
        n_buildings = len(self.buildings_objects)
        counter = 0
        for bd_id, building_info in self.buildings_info.items():
            if counter % 10 == 0:
                print(f"{counter} buildings simulated out of {n_buildings}")
            counter += 1
            
            self.buildings_objects[bd_id].info_geojson = building_info
            if building_info["ElLoad_cal"] != 1.:
                el_load_m = building_info["ElLoad_cal"]
            else:
                el_load_m = 1.
    
            info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
            info["Name"] = self.buildings_objects[bd_id].name
            if print_single_building_results:
                
                
                if not np.isnan(building_info["Ele Cooling [kWh]"]):
                    cal_res = self.buildings_objects[bd_id].calibrate_summer_mfr_c_setp(building_info["Ele Cooling [kWh]"],
                                                                  self.weather_file,
                                                                  limits_setpoint=limits_setpoint,
                                                                  limits_mfr_m=limits_mfr_m
                                                                  )
                else: 
                    cal_res=[[1.,26.],None,None,None,None,None]
                
                    
                self.buildings_objects[bd_id].update_calibration_params(self.weather_file,
                                                                        mfr_m = cal_res[0][0], 
                                                                        c_t_set=cal_res[0][1], 
                                                                        elec_loads_m=el_load_m) 
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=self.output_folder,
                                                                  output_type=output_type)
                info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
                info["Name"] = self.buildings_objects[bd_id].name
                info["MFR_M_cal"] = cal_res[0][0]
                info["TSetC_cal"] = cal_res[0][1]
                info["Calibration output"] = cal_res[4]
                info["Calibration iterations"] = cal_res[3]
                info["Calibration status"] = cal_res[5]
            else:
                
                if not np.isnan(building_info["Ele Cooling [kWh]"]):
                    cal_res = self.buildings_objects[bd_id].calibrate_summer_mfr_c_setp(building_info["Ele Cooling [kWh]"], self.weather_file,
                                                                  limits_setpoint=limits_setpoint,
                                                                  limits_mfr_m=limits_mfr_m
                                                                )
                 
                else: 
                    cal_res=[[1.,26.],None,None,None,None,None]
                    
                self.buildings_objects[bd_id].update_calibration_params(self.weather_file,
                                                                        mfr_m = cal_res[0][0], 
                                                                        c_t_set=cal_res[0][1], 
                                                                        elec_loads_m=el_load_m) 
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=None)
                info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
                info["Name"] = self.buildings_objects[bd_id].name
                info["MFR_M_cal"] = cal_res[0][0]
                info["TSetC_cal"] = cal_res[0][1]
                info["Calibration output"] = cal_res[4]
                info["Calibration iterations"] = cal_res[3]
                info["Calibration status"] = cal_res[5]
            results.index = index
            monthly = results.resample("M").sum()
            gas_consumption = monthly[[col for col in monthly.columns if "gas consumption" in col[0]]].sum(axis=1)
            el_consumption = monthly[[col for col in monthly.columns if "electric consumption" in col[0]]].sum(axis=1)
            oil_consumption = monthly[[col for col in monthly.columns if "oil consumption" in col[0]]].sum(axis=1)
            wood_consumption = monthly[[col for col in monthly.columns if "wood consumption" in col[0]]].sum(axis=1)
            gas_consumption["Total"] = gas_consumption.sum()
            el_consumption["Total"] = el_consumption.sum()
            oil_consumption["Total"] = oil_consumption.sum()
            wood_consumption["Total"] = wood_consumption.sum()
            for i in gas_consumption.index:
                if i == "Total":
                    info[f"{i} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i} wood consumption [kg]"] = wood_consumption.loc[i]
                else:
                    info[f"{i.month_name()} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i.month_name()} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i.month_name()} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i.month_name()} wood consumption [kg]"] = wood_consumption.loc[i]
            final_results[bd_id] = info
            district_hourly_results["Gas consumption [Nm3]"] += results["Heating system gas consumption [Nm3]"].iloc[:,
                                                                0]
            district_hourly_results["Oil consumption [L]"] += results["Heating system oil consumption [L]"].iloc[:, 0]
            district_hourly_results["Wood consumption [kg]"] += results["Heating system wood consumption [kg]"].iloc[:,
                                                                0]
            district_hourly_results["Electric consumption [Wh]"] += results[
                                                                        "Heating system electric consumption [Wh]"].iloc[
                                                                    :, 0] \
                                                                    + results[
                                                                          "Cooling system electric consumption [Wh]"].iloc[
                                                                      :, 0] \
                                                                    + results[
                                                                          "Appliances electric consumption [Wh]"].iloc[
                                                                      :, 0]
            district_hourly_results["TZ sensible load [W]"] += results["TZ sensible load [W]"].iloc[:, 0]
            district_hourly_results["TZ latent load [W]"] += results["TZ latent load [W]"].iloc[:, 0]
            district_hourly_results["TZ AHU pre heater load [W]"] += results["TZ AHU pre heater load [W]"].iloc[:, 0]
            district_hourly_results["TZ AHU post heater load [W]"] += results["TZ AHU post heater load [W]"].iloc[:, 0]
            district_hourly_results["TZ DHW demand [W]"] += results["TZ DHW demand [W]"].iloc[:, 0]
    
        district_hourly_results.to_csv(os.path.join(self.output_folder, "District_hourly_summary.csv"), sep=";")
        bd_summary = pd.DataFrame.from_dict(final_results, orient="index")
        # bd_summary.to_csv(os.path.join(self.output_folder,"Buildings_summary.csv"), sep =";")
        bd_summary.drop(["Name"], axis=1, inplace=True)
        self.output_geojson.set_index("new_id", drop=True, inplace=True)
        new_geojson = pd.concat([self.output_geojson.loc[bd_summary.index],bd_summary],axis=1)
        new_geojson.to_file(os.path.join(self.output_folder, "Buildings_summary.geojson"), driver="GeoJSON")
        new_geojson.drop("geometry", axis=1).to_csv(os.path.join(self.output_folder, "Buildings_summary.csv"), sep=";")
    
        print(f"Standard simulation : {(time.time() - start) / 60:0.2f} min")

    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################

    def calibrate_gas(self, print_single_building_results=False, output_type="parquet",
                  limits_setpoint = [18,22],
                  limits_fwalls = [0.75,1.25]
                  ):
        """Simulation of the whole city, and memorization and stamp of results.
    
        Parameters
        ----------
        print_single_building_results : bool, default False
            If True, the prints a file with time step results for each building.
            USE CAREFULLY: It might fill a lot of disk space
        output_type : str, default 'parquet'
            It can be either csv (more suitable for excel but more disk space) or parquet (more efficient but not readable from excel)
        """
    
        import time
        start = time.time()
        # parallel simulation commented
        # bd_parallel_list = [[bd, self.weather_file, self.output_folder] for bd in self.buildings_objects.values()]
        # def bd_parallel_solve(x):
        #     bd, weather, out_fold = x
        #     bd.simulate(weather, output_folder=out_fold)
        # with concurrent.futures.ThreadPoolExecutor() as executor:
        #     bd_executor = executor.map(bd_parallel_solve, bd_parallel_list)
        #
        # print(f"Parallel simulation : {(time.time() - start)/60:0.2f} min")
        # start = time.time()
    
        final_results = {}
    
        index = pd.date_range(start=CONFIG.start_date, periods=CONFIG.number_of_time_steps, freq=f"{CONFIG.time_step}s")
        district_hourly_results = pd.DataFrame(0., index=index, columns=[
            "Gas consumption [Nm3]",
            "Electric consumption [Wh]",
            "Oil consumption [L]",
            "Wood consumption [kg]",
            "TZ sensible load [W]",
            "TZ latent load [W]",
            "TZ AHU pre heater load [W]",
            "TZ AHU post heater load [W]",
            "TZ DHW demand [W]",
        ])
        n_buildings = len(self.buildings_objects)
        counter = 0
        for bd_id, building_info in self.buildings_info.items():
            if counter % 10 == 0:
                print(f"{counter} buildings simulated out of {n_buildings}")
            counter += 1
            
            self.buildings_objects[bd_id].info_geojson = building_info
            if building_info["ElLoad_cal"] != 1.:
                el_load_m = building_info["ElLoad_cal"]
            else:
                el_load_m = 1.
            if building_info["MFR_M_cal"] != 1.:
                mfr = building_info["MFR_M_cal"]
            else:
                mfr = 1.
            if building_info["TSetC_cal"] != 26.:
                c_t_set = building_info["TSetC_cal"]
            else:
                c_t_set = 26.
                
    
            info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
            info["Name"] = self.buildings_objects[bd_id].name
            
            gas_cons = building_info["Gas Heating [Nm3]"] + building_info["Gas Non Heating [Nm3]"]
            
            if print_single_building_results:
                
                
                if not np.isnan(gas_cons):
                    cal_res = self.buildings_objects[bd_id].calibrate_winter(gas_cons,
                                                                  self.weather_file,
                                                                  limits_setpoint=limits_setpoint,
                                                                  limits_fwalls=limits_fwalls
                                                                  )
                else: 
                    cal_res=[[1.,20.],None,None,None,None,None]
                
                    
                self.buildings_objects[bd_id].update_calibration_params(self.weather_file,
                                                                        ext_wall_coef = cal_res[0][0], 
                                                                        h_t_set = cal_res[0][1],
                                                                        mfr_m = mfr, 
                                                                        c_t_set=c_t_set, 
                                                                        elec_loads_m=el_load_m) 
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=self.output_folder,
                                                                  output_type=output_type)
                info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
                info["Name"] = self.buildings_objects[bd_id].name
                info["FWalls_M_cal"] = cal_res[0][0]
                info["TSetH_cal"] = cal_res[0][1]
                info["Calibration output"] = cal_res[4]
                info["Calibration iterations"] = cal_res[3]
                info["Calibration status"] = cal_res[5]
            else:
                
                if not np.isnan(gas_cons):
                    cal_res = self.buildings_objects[bd_id].calibrate_winter(gas_cons,
                                                                  self.weather_file,
                                                                  limits_setpoint=limits_setpoint,
                                                                  limits_fwalls=limits_fwalls
                                                                )
                 
                else: 
                    cal_res=[[1.,20.],None,None,None,None,None]
                    
                self.buildings_objects[bd_id].update_calibration_params(self.weather_file,
                                                                        ext_wall_coef = cal_res[0][0], 
                                                                        h_t_set = cal_res[0][1],
                                                                        mfr_m = mfr, 
                                                                        c_t_set=c_t_set, 
                                                                        elec_loads_m=el_load_m) 
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=None)
                info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
                info["Name"] = self.buildings_objects[bd_id].name
                info["FWalls_M_cal"] = cal_res[0][0]
                info["TSetH_cal"] = cal_res[0][1]
                info["Calibration output"] = cal_res[4]
                info["Calibration iterations"] = cal_res[3]
                info["Calibration status"] = cal_res[5]
            results.index = index
            monthly = results.resample("M").sum()
            gas_consumption = monthly[[col for col in monthly.columns if "gas consumption" in col[0]]].sum(axis=1)
            el_consumption = monthly[[col for col in monthly.columns if "electric consumption" in col[0]]].sum(axis=1)
            oil_consumption = monthly[[col for col in monthly.columns if "oil consumption" in col[0]]].sum(axis=1)
            wood_consumption = monthly[[col for col in monthly.columns if "wood consumption" in col[0]]].sum(axis=1)
            gas_consumption["Total"] = gas_consumption.sum()
            el_consumption["Total"] = el_consumption.sum()
            oil_consumption["Total"] = oil_consumption.sum()
            wood_consumption["Total"] = wood_consumption.sum()
            for i in gas_consumption.index:
                if i == "Total":
                    info[f"{i} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i} wood consumption [kg]"] = wood_consumption.loc[i]
                else:
                    info[f"{i.month_name()} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i.month_name()} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i.month_name()} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i.month_name()} wood consumption [kg]"] = wood_consumption.loc[i]
            final_results[bd_id] = info
            district_hourly_results["Gas consumption [Nm3]"] += results["Heating system gas consumption [Nm3]"].iloc[:,
                                                                0]
            district_hourly_results["Oil consumption [L]"] += results["Heating system oil consumption [L]"].iloc[:, 0]
            district_hourly_results["Wood consumption [kg]"] += results["Heating system wood consumption [kg]"].iloc[:,
                                                                0]
            district_hourly_results["Electric consumption [Wh]"] += results[
                                                                        "Heating system electric consumption [Wh]"].iloc[
                                                                    :, 0] \
                                                                    + results[
                                                                          "Cooling system electric consumption [Wh]"].iloc[
                                                                      :, 0] \
                                                                    + results[
                                                                          "Appliances electric consumption [Wh]"].iloc[
                                                                      :, 0]
            district_hourly_results["TZ sensible load [W]"] += results["TZ sensible load [W]"].iloc[:, 0]
            district_hourly_results["TZ latent load [W]"] += results["TZ latent load [W]"].iloc[:, 0]
            district_hourly_results["TZ AHU pre heater load [W]"] += results["TZ AHU pre heater load [W]"].iloc[:, 0]
            district_hourly_results["TZ AHU post heater load [W]"] += results["TZ AHU post heater load [W]"].iloc[:, 0]
            district_hourly_results["TZ DHW demand [W]"] += results["TZ DHW demand [W]"].iloc[:, 0]
    
        district_hourly_results.to_csv(os.path.join(self.output_folder, "District_hourly_summary.csv"), sep=";")
        bd_summary = pd.DataFrame.from_dict(final_results, orient="index")
        # bd_summary.to_csv(os.path.join(self.output_folder,"Buildings_summary.csv"), sep =";")
        bd_summary.drop(["Name"], axis=1, inplace=True)
        self.output_geojson.set_index("new_id", drop=True, inplace=True)
        new_geojson = pd.concat([self.output_geojson.loc[bd_summary.index],bd_summary],axis=1)
        if 'level_0' in new_geojson.columns:
            new_geojson.drop('level_0', axis = 1, inplace = True)
        new_geojson.to_file(os.path.join(self.output_folder, "Buildings_summary.geojson"), driver="GeoJSON")
        new_geojson.drop("geometry", axis=1).to_csv(os.path.join(self.output_folder, "Buildings_summary.csv"), sep=";")
    
        print(f"Standard simulation : {(time.time() - start) / 60:0.2f} min")

    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
    ##########################################################################
    
    
    def calibrate(self, print_single_building_results=False, output_type="parquet",
                  limits_setpoint = [18,22],
                  limits_fwalls = [0.75,1.25]):
        """Simulation of the whole city, and memorization and stamp of results.

        Parameters
        ----------
        print_single_building_results : bool, default False
            If True, the prints a file with time step results for each building.
            USE CAREFULLY: It might fill a lot of disk space
        output_type : str, default 'parquet'
            It can be either csv (more suitable for excel but more disk space) or parquet (more efficient but not readable from excel)
        """

        import time
        start = time.time()
        # parallel simulation commented
        # bd_parallel_list = [[bd, self.weather_file, self.output_folder] for bd in self.buildings_objects.values()]
        # def bd_parallel_solve(x):
        #     bd, weather, out_fold = x
        #     bd.simulate(weather, output_folder=out_fold)
        # with concurrent.futures.ThreadPoolExecutor() as executor:
        #     bd_executor = executor.map(bd_parallel_solve, bd_parallel_list)
        #
        # print(f"Parallel simulation : {(time.time() - start)/60:0.2f} min")
        # start = time.time()

        final_results = {}

        index = pd.date_range(start=CONFIG.start_date, periods=CONFIG.number_of_time_steps, freq=f"{CONFIG.time_step}s")
        district_hourly_results = pd.DataFrame(0., index=index, columns=[
            "Gas consumption [Nm3]",
            "Electric consumption [Wh]",
            "Oil consumption [L]",
            "Wood consumption [kg]",
            "TZ sensible load [W]",
            "TZ latent load [W]",
            "TZ AHU pre heater load [W]",
            "TZ AHU post heater load [W]",
            "TZ DHW demand [W]",
        ])
        n_buildings = len(self.buildings_objects)
        counter = 0
        for bd_id, building_info in self.buildings_info.items():
            if counter % 10 == 0:
                print(f"{counter} buildings simulated out of {n_buildings}")
            counter += 1

            info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
            info["Name"] = self.buildings_objects[bd_id].name
            if print_single_building_results:
                cal_res = self.buildings_objects[bd_id].calibrate(self.buildings_info[bd_id]["STDMC_2021"],
                                                                  self.weather_file,
                                                                  limits_setpoint=limits_setpoint,
                                                                  limits_fwalls=limits_fwalls
                                                                  )
                self.buildings_objects[bd_id].update_calibration_params(self.weather_file, ext_wall_coef = cal_res[0][0], h_t_set = cal_res[0][1])
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=self.output_folder,
                                                                  output_type=output_type)
                info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
                info["Name"] = self.buildings_objects[bd_id].name
                info["ExtWallCoef_cal"] = cal_res[0][0]
                info["TSet_cal"] = cal_res[0][1]
                info["Calibration output"] = cal_res[4]
                info["Calibration iterations"] = cal_res[3]
                info["Calibration status"] = cal_res[5]
            else:
                cal_res = self.buildings_objects[bd_id].calibrate(self.buildings_info[bd_id]["STDMC_2021"], self.weather_file,
                                                                limits_setpoint = limits_setpoint,
                                                              limits_fwalls = limits_fwalls
                                                                )
                self.buildings_objects[bd_id].update_calibration_params(self.weather_file, ext_wall_coef = cal_res[0][0], h_t_set = cal_res[0][1])
                results = self.buildings_objects[bd_id].simulate(self.weather_file, output_folder=None)
                info = self.buildings_objects[bd_id]._thermal_zones_list[0].get_zone_info()
                info["Name"] = self.buildings_objects[bd_id].name
                info["ExtWallCoef_cal"] = cal_res[0][0]
                info["TSet_cal"] = cal_res[0][1]
                info["Calibration output"] = cal_res[4]
                info["Calibration iterations"] = cal_res[3]
                info["Calibration status"] = cal_res[5]
            results.index = index
            monthly = results.resample("M").sum()
            gas_consumption = monthly[[col for col in monthly.columns if "gas consumption" in col[0]]].sum(axis=1)
            el_consumption = monthly[[col for col in monthly.columns if "electric consumption" in col[0]]].sum(axis=1)
            oil_consumption = monthly[[col for col in monthly.columns if "oil consumption" in col[0]]].sum(axis=1)
            wood_consumption = monthly[[col for col in monthly.columns if "wood consumption" in col[0]]].sum(axis=1)
            gas_consumption["Total"] = gas_consumption.sum()
            el_consumption["Total"] = el_consumption.sum()
            oil_consumption["Total"] = oil_consumption.sum()
            wood_consumption["Total"] = wood_consumption.sum()
            for i in gas_consumption.index:
                if i == "Total":
                    info[f"{i} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i} wood consumption [kg]"] = wood_consumption.loc[i]
                else:
                    info[f"{i.month_name()} gas consumption [Nm3]"] = gas_consumption.loc[i]
                    info[f"{i.month_name()} electric consumption [Wh]"] = el_consumption.loc[i]
                    info[f"{i.month_name()} oil consumption [L]"] = oil_consumption.loc[i]
                    info[f"{i.month_name()} wood consumption [kg]"] = wood_consumption.loc[i]
            final_results[bd_id] = info
            district_hourly_results["Gas consumption [Nm3]"] += results["Heating system gas consumption [Nm3]"].iloc[:,
                                                                0]
            district_hourly_results["Oil consumption [L]"] += results["Heating system oil consumption [L]"].iloc[:, 0]
            district_hourly_results["Wood consumption [kg]"] += results["Heating system wood consumption [kg]"].iloc[:,
                                                                0]
            district_hourly_results["Electric consumption [Wh]"] += results[
                                                                        "Heating system electric consumption [Wh]"].iloc[
                                                                    :, 0] \
                                                                    + results[
                                                                          "Cooling system electric consumption [Wh]"].iloc[
                                                                      :, 0] \
                                                                    + results[
                                                                          "Appliances electric consumption [Wh]"].iloc[
                                                                      :, 0]
            district_hourly_results["TZ sensible load [W]"] += results["TZ sensible load [W]"].iloc[:, 0]
            district_hourly_results["TZ latent load [W]"] += results["TZ latent load [W]"].iloc[:, 0]
            district_hourly_results["TZ AHU pre heater load [W]"] += results["TZ AHU pre heater load [W]"].iloc[:, 0]
            district_hourly_results["TZ AHU post heater load [W]"] += results["TZ AHU post heater load [W]"].iloc[:, 0]
            district_hourly_results["TZ DHW demand [W]"] += results["TZ DHW demand [W]"].iloc[:, 0]

        district_hourly_results.to_csv(os.path.join(self.output_folder, "District_hourly_summary.csv"), sep=";")
        bd_summary = pd.DataFrame.from_dict(final_results, orient="index")
        # bd_summary.to_csv(os.path.join(self.output_folder,"Buildings_summary.csv"), sep =";")
        bd_summary.drop(["Name"], axis=1, inplace=True)
        self.output_geojson.set_index("new_id", drop=True, inplace=True)
        new_geojson = pd.concat([self.output_geojson.loc[bd_summary.index],bd_summary],axis=1)
        if 'level_0' in new_geojson.columns:
            new_geojson.drop('level_0', axis = 1, inplace = True)
        new_geojson.to_file(os.path.join(self.output_folder, "Buildings_summary.geojson"), driver="GeoJSON")
        new_geojson.drop("geometry", axis=1).to_csv(os.path.join(self.output_folder, "Buildings_summary.csv"), sep=";")

        print(f"Standard simulation : {(time.time() - start) / 60:0.2f} min")