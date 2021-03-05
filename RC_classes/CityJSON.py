'''IMPORTING MODULES'''

import sys
import math
import os
import geopandas as gpd
import numpy as np
import time as tm
import pandas as pd
from warnings import warn as wrn
from cjio import cityjson
from RC_classes.WeatherData import SolarPosition, Weather
from RC_classes.thermalZone import Building, Complex
from RC_classes.Envelope import loadEnvelopes
from RC_classes.Geometry import Surface
from RC_classes.Climate import UrbanCanyon

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions
 
def cityobj(p):
    
    '''
    This function permits to create the city object when mode == cityjson
    
    Parameters
    ----------
    p : str
        path of the cityJSON file
    
    Returns
    -------
    c : cjio.cityjson.CityJSON
        city object
        
    '''
    
    # Check input data type
    
    if not isinstance(p, str):
        raise TypeError(f'Ops... cityobj function, p is not a str: p {p}')

    # CityJSON object creation
    try:
        with open(p, 'r') as f:
            c = cityjson.CityJSON(file=f)
            return c
    except FileNotFoundError:
        raise FileNotFoundError(f'Ops...Cityjson object not found: {p}. Give a proper path')
        


def createBuilding(name,bd,vertList,mode,n_Floors,envelopes,weather):
    
    '''
    This function allows to identify buildings and their attributes when
    mode == cityjson
    
    Parameters
    ----------
    name : str
        name of the district's buildings 
    bd : dict
        dictionary for each building containing attributes and info
    vertList : list
        list of all buildings vertices
    mode : str
        cityjson or geojson mode calculation
    n_Floors : int
        initialization of the building's number of floors
    envelopes : dict
        envelope information for each age class category
    weather : RC_classes.WeatherData.Weather obj
        object of the class weather WeatherData module
    
    Returns
    -------
    Building : class
    '''

    # Check input data type

    if not isinstance(name, str):
        raise TypeError(f'Ops... createBuilding function, name is not a str: name {name}')
    if not isinstance(bd, dict):
        raise TypeError(f'Ops... createBuilding function, bd is not a dict: bd {bd}')
    if not isinstance(vertList, list):
        raise TypeError(f'Ops... createBuilding function, vertList is not a list: vertList {vertList}')
    if not isinstance(mode, str):
        raise TypeError(f'Ops... createBuilding function, mode is not a str: mode {mode}')
    if not isinstance(n_Floors, int):
        raise TypeError(f'Ops... createBuilding function, n_Floors is not a int: n_Floors {n_Floors}')
    if not isinstance(envelopes, dict):
        raise TypeError(f'Ops... createBuilding function, envelopes is not a dict: envelopes {envelopes}')
    if not isinstance(weather, Weather):
        raise TypeError(f'Ops... createBuilding function, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
            
    # Check input data quality

    for vtx in vertList:
        if not isinstance(vtx, list):
            raise TypeError(f'Ops... createBuilding function, an input is not a list: input {vtx}')
        if len(vtx) != 3:
            raise TypeError(f'Ops... createBuilding function, a vertex is not a list of 3 components: input {vtx}')
        try:
            vtx[0] = float(vtx[0])
            vtx[1] = float(vtx[1])
            vtx[2] = float(vtx[2])
        except ValueError:
            raise ValueError(f'Ops... createBuilding function, a coordinate is not a float: input {vtx}')     
    if not bool(envelopes):
        wrn(f"\n\createBuilding function, the envelopes dictionary is empty..... envelopes {envelopes}\n")
    
    
    # Setting the attributes of the building
    age = bd['attributes']['Age']                                              # Age-class of the building
    use = bd['attributes']['Use']                                              # End-use of the building
    try: 
        H_plant = bd['attributes']['H_Plant']                                  # Heating plant of the building
    except KeyError:
        H_plant = 'IdealLoad'
    try:
        C_plant = bd['attributes']['C_Plant']                                  # Cooling plant of the building
    except KeyError:
        C_plant = 'IdealLoad'
    
    for geo in bd['geometry']:
        if geo['type']=='MultiSurface':
            boundaries = geo['boundaries']
            surfaces=[]
            for surface in boundaries:
                for subsurface in surface:
                    surf=[]
                    for vert_id in subsurface:
                        surf.append(vertList[vert_id])    
                    surfaces.append(surf)
        
        # Building class initialization
        return Building(name,mode,surfaces,n_Floors,use,envelopes,age,1,1,H_plant,C_plant,weather)


#%% ---------------------------------------------------------------------------------------------------
#%% JsonCity class

class JsonCity():
    
    '''
    This class manages the city via json input file: geojson or cityjson
    
    Methods
        init
        printInfo
        create_urban_canyon
        shading_effect
        paramsandloads
        designdays
        cityplants
        citysim
        complexmerge
        
    '''
    
    # class variables
    T_w_0 = 15.                                                                # Starting average temperature of the walls [°C]
    T_out_inf = 15.                                                            # Starting average temperature outgoing the buildings [°C]
    T_out_AHU = 15.                                                            # Starting average temperature outgoing the Air Handling units [°C]
    V_0_inf = 0.                                                               # Starting average volumetric flow rate outgoing the buildings due to the inflitrations[m3/s]
    V_0_vent = 0.                                                              # Starting average volumetric flow rate outgoing the Air Handling units [m3/s]
    H_waste_0 = 0.                                                             # Starting waste heating rejected by external condensers [kW]
    
    def __init__(self,json_path,envelopes,weather,model = '1C', mode='cityjson'):
        
        '''
        This method creates the city from the json file
        
        Parameters
        ----------
        json_path : str
            path of the cityJSON file
        envelopes : dict
            envelope information for each age class category
        weather : RC_classes.WeatherData.Weather obj
            object of the class weather WeatherData module
        model : str
            Resistance-Capacitance model used
        mode : str
            cityjson or geojson mode calculation
        
        Returns
        -------
        None
        
        '''
        
        # Check input data type
    
        if not isinstance(json_path, str):
            raise TypeError(f'Ops... JsonCity class, json_path is not a str: json_path {json_path}')
        if not isinstance(model, str):
            raise TypeError(f'Ops... JsonCity class, model is not a str: model {model}')
        if not isinstance(envelopes, dict):
            raise TypeError(f'Ops... JsonCity class, envelopes is not a dict: envelopes {envelopes}')
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
        if not isinstance(mode, str):
            raise TypeError(f'Ops... JsonCity class, mode is not a str: mode {mode}')

        # Check input data quality
        
        if model != '1C' and  model != '2C':
            wrn(f"\n\JsonCity class, the model doesn't exist..... model {model}\n")
        if not bool(envelopes):
            wrn(f"\n\JsonCity class, the envelopes dictionary is empty..... envelopes {envelopes}\n")

        # Creation of the city
        self.model = model
        
        # Case of cityJSON file availability:
        if mode=='cityjson':
            self.mode = mode
            self.n_Floors = 0
            self.city = cityobj(json_path)
            self.jsonBuildings= {}
            [(self.jsonBuildings.update({i:self.city.j['CityObjects'][i]})) for i in self.city.j['CityObjects'] if self.city.j['CityObjects'][i]['type']=='Building']
            self.buildings = {}
            self.complexes = {}
            for bd in self.jsonBuildings:
                if not(isinstance(self.city.j['vertices'], list)):
                    raise ValueError('json file vertices are not a list')
                if not(isinstance(self.n_Floors, int)):
                    raise ValueError('The floor number is not an integer')
                self.buildings[bd]=createBuilding(bd,self.jsonBuildings[bd],self.city.j['vertices'],self.mode,self.n_Floors,envelopes,weather)
        
        # Case of GeoJSON file availability:
        elif mode=='geojson':
            self.mode = mode
            self.city = gpd.read_file(json_path)
            self.jsonBuildings= {}
            self.buildings = {}
            self.complexes = {}
            
            # Extrusion from the footprint operation
            for i in self.city.index:
                self.jsonBuildings[self.city.loc[i]['id']] = self.city.loc[i].to_dict()
                # https://gis.stackexchange.com/questions/287306/list-all-polygon-vertices-coordinates-using-geopandas
                g = [i for i in self.city.loc[i].geometry]
                x,y = g[0].exterior.coords.xy
                coords = np.dstack((x,y)).tolist()
                coords=coords[0]
                coords.pop()
                build_surf=[]
                pavimento = []
                soffitto = []
                z_pav = 0
                z_soff = self.city.loc[i].altezza
                for n in range(len(coords)):
                    pavimento.append(coords[n]+[z_pav])
                    soffitto.append(coords[-n]+[z_soff])  
                for n in range(len(coords)):
                    build_surf.append([coords[n-1]+[z_soff],\
                                        coords[n]+[z_soff],\
                                        coords[n]+[z_pav],\
                                        coords[n-1]+[z_pav]])\
                        
                build_surf.append(pavimento)
                build_surf.append(soffitto)
                self.rh_net = self.city.loc[i]['rh_net']
                self.rh_gross = self.city.loc[i]['rh_gross']
                
                try:
                    self.heating_plant = self.city.loc[i]['Heating_plant']
                except KeyError:
                    self.heating_plant = 'IdealLoad'
                try:
                    self.cooling_plant = self.city.loc[i]['Cooling_plant']
                except KeyError:
                    self.cooling_plant = 'IdealLoad'
                    
                self.buildings[self.city.loc[i]['id']]=Building(self.city.loc[i]['nome'], 
                                                                self.mode, 
                                                                build_surf,
                                                                self.city.loc[i]['n_piani'],
                                                                self.city.loc[i]['Uso'],                                                                                                                             
                                                                envelopes,
                                                                self.city.loc[i]['età'],
                                                                self.rh_net,
                                                                self.rh_gross,
                                                                self.heating_plant,
                                                                self.cooling_plant,
                                                                weather)
                
        else:
            sys.exit('Set a proper mode')


    def printInfo(self):
        
        '''
        This method prints some variable of the district
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        '''
        
        for i in self.buildings.values():
            i.printBuildingInfo()


    def create_urban_canyon(self,sim_time,calc,data):
        
        '''
        This method allows to evaluate the Urban Heat Island effect
        
        Parameters
        ----------
        sim_time : list
            list containing timestep per hour and hours of simulation info  
        calc : bool
            True if the calculation is performed, False in the opposite case
        data : dict
            dictionary containing general and specific district information
        
        Returns
        -------
        None
        
        '''
        
        # Check input data type
        
        if not isinstance(sim_time, list):
            raise TypeError(f'Ops... JsonCity class - create_urban_canyon, sim_time is not a list: sim_time {sim_time}')
        if not isinstance(calc, bool):
            raise TypeError(f'Ops... JsonCity class - create_urban_canyon, calc is not a bool: calc {calc}')
        if not isinstance(data, dict):
            raise TypeError(f'Ops... JsonCity class - create_urban_canyon, data is not a dict: data {data}')
        
        # Check input data quality
        
        if sim_time[0] > 6:
            wrn(f"\n\JsonCity class - create_urban_canyon, more than 6 timestper per hours - it would be better to reduce ts {sim_time[0]}\n")
        
        # Urban Heat Island evaluation
        if calc:
            self.urban_canyon_calc = calc
            self.bd_ext_walls = np.array([],dtype = float)
            tot_wall_glazed_area = .0
            tot_wall_opaque_area = .0
            h_builidngs = np.array([], dtype = float)       
            building_footprints = np.array([], dtype = float)
            for bd in self.buildings.values():
                self.bd_ext_walls = np.append(self.bd_ext_walls, bd.extWallOpaqueArea + bd.extWallWinArea)
                tot_wall_opaque_area += bd.extWallOpaqueArea
                tot_wall_glazed_area += bd.extWallWinArea
                building_footprints = np.append(building_footprints, bd.footprint) 
                h_builidngs = np.append(h_builidngs, bd.buildingHeight)
            
            data['tot_wall_opaque_area'] = tot_wall_opaque_area
            data['tot_wall_glazed_area'] = tot_wall_glazed_area
            data['tot_footprint'] = np.sum(building_footprints)
            data['h_builidng_average'] = np.average(h_builidngs, weights = building_footprints)
            data['VH_urb'] = (tot_wall_opaque_area + tot_wall_glazed_area)/data['Area']
            data['Buildings_density'] = data['tot_footprint'] /data['Area']
            
            # Check data quality
            
            if ((data['Buildings_density'] + data['Vegetation_density']) > 0.9) :
                print('WARNING: building density with vegetation density higher than 0.9: check vegetation density')
                if ((data['Buildings_density'] + data['Vegetation_density']) >= .99):
                    data['Vegetation_density'] = 0.99 - data['Buildings_density'] 
            if data['tot_footprint'] > data['Area']:
                sys.exit('Error: the sum of building footprints is higher than City area. Correct city area')
            
            self.urban_canyon_data = data
            self.urban_canyon = UrbanCanyon(sim_time,self.urban_canyon_data)
        else:
            self.urban_canyon_calc = calc
            self.urban_canyon_data = None 
   
    
    def shading_effect(self,Solar_position,mode = 'cityjson',
                       toll_az = 80.,
                       toll_dist = 100. ,
                       toll_theta = 80.,
                       R_f = 0.):
        
        '''
        This method takes into account the shading effect between buildings surfaces
        
        Parameters
        ----------
        Solar_position : WeatherData.SolarPosition
            Object containing information about the position of the sun over the time of simulation
        mode : str
            cityjson or geojson mode calculation
        toll_az : float
            semi-tollerance on azimuth of shading surface [°]
        toll_dist : float
            tollerance on distance of shading surface [m]
        toll_theta : float
            semi-tollerance on position of the shading surfaces [°]
        R_f : float
            Reduction factor of the direct solar radiation due to the shading effect [0-1]
        
        Returns
        -------
        None
        
        '''
        
        # Check input data type
        
        if not isinstance(mode, str):
            raise TypeError(f'Ops... JsonCity class - shading_effect, mode is not a str: mode {mode}')
        if not isinstance(toll_az, float):
            raise TypeError(f'Ops... JsonCity class - shading_effect, toll_az is not a float: toll_az {toll_az}')
        if not isinstance(toll_dist, float):
            raise TypeError(f'Ops... JsonCity class - shading_effect, toll_dist is not a float: toll_dist {toll_dist}')
        if not isinstance(toll_theta, float):
            raise TypeError(f'Ops... JsonCity class - shading_effect, toll_theta is not a float: toll_theta {toll_theta}')
        if not isinstance(Solar_position, SolarPosition):
            raise TypeError(f'Ops... JsonCity class - shading_effect, Solar_position is not a SolarPosition object: Solar_position {Solar_position}')
        if not isinstance(R_f, float):
            raise TypeError(f'Ops... JsonCity class - shading_effect, R_f is not a float: R_f {R_f}')
        
        # Check input data quality
        
        if mode != 'geojson' and  mode != 'cityjson':
            wrn(f"\n\JsonCity class - shading_effect, the mode doesn't exist..... mode {mode}\n")
        if toll_az < 0.0 or toll_az > 90.0:
            wrn(f"\n\JsonCity class - shading_effect, toll_az is out of range [0-90]..... toll_az {toll_az}\n")
        if toll_dist < 0.0 or toll_dist > 200.0:
            wrn(f"\n\JsonCity class - shading_effect, toll_dist is out of range [0-200]..... toll_dist {toll_dist}\n")
        if toll_theta < 0.0 or toll_theta > 90.0:
            wrn(f"\n\JsonCity class - shading_effect, toll_theta is out of range [0-90]..... toll_theta {toll_theta}\n")
        if R_f < 0.0 or R_f > 1.0:
            wrn(f"\n\JsonCity class - shading_effect, R_f is out of range [0-1]..... R_f {R_f}\n")
            

        # SECTION 1: All surfaces are compared and potentially shading surfaces are stored
        self.all_Vertsurf = []
        
        if mode == 'cityjson':
            for bd in self.buildings.keys():
                self.all_Vertsurf.extend(self.buildings[str(bd)].Vertsurf)
        if mode == 'geojson':
            for i in self.city.index:
                self.all_Vertsurf.extend(self.buildings[i+1].Vertsurf)
       
        # Each surface is compared with all the others
        for x in range(len(self.all_Vertsurf)):
            for y in range(len(self.all_Vertsurf)):
                if y > x:
                    
                    # Calculation of the distance between the centroids of the two surfaces under examination
                    dist = math.sqrt((self.all_Vertsurf[x][0].centroid_coord[0]-self.all_Vertsurf[y][0].centroid_coord[0])**2+(self.all_Vertsurf[x][0].centroid_coord[1]-self.all_Vertsurf[y][0].centroid_coord[1])**2)
                    if dist == 0.0:
                        pass
                    else:
                        
                        # Calculation of the vector direction between the centroids of the two surfaces under examination
                        theta_xy = np.degrees(np.arccos((self.all_Vertsurf[y][0].centroid_coord[0]-self.all_Vertsurf[x][0].centroid_coord[0])/dist))
                        theta = -(theta_xy + 90)
                        if self.all_Vertsurf[y][0].centroid_coord[1] < self.all_Vertsurf[x][0].centroid_coord[1]:
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
                        if self.all_Vertsurf[x][0].azimuth < 0:
                            azimuth_opp = self.all_Vertsurf[x][0].azimuth + 180
                            azimuth_opp_max = azimuth_opp + toll_az
                            azimuth_opp_min = azimuth_opp - toll_az
                            theta_max = self.all_Vertsurf[x][0].azimuth + toll_theta
                            theta_min = self.all_Vertsurf[x][0].azimuth - toll_theta
                            if theta_min < -180:
                                theta_min = theta_min + 360
                                if theta_min < theta < 180 or -180 < theta < theta_max:
                                    if azimuth_opp_max > 180:
                                        azimuth_opp_max = azimuth_opp_max - 360
                                        if azimuth_opp_min < self.all_Vertsurf[y][0].azimuth <= 180 or -180 <= self.all_Vertsurf[y][0].azimuth < azimuth_opp_max:
                                            self.all_Vertsurf[x][1].extend([[dist,y]])
                                            self.all_Vertsurf[y][1].extend([[dist,x]])
                                    else:
                                        if azimuth_opp_min < self.all_Vertsurf[y][0].azimuth < azimuth_opp_max:
                                            self.all_Vertsurf[x][1].extend([[dist,y]])
                                            self.all_Vertsurf[y][1].extend([[dist,x]])
                            else:
                                if theta_min < theta < theta_max:
                                    if azimuth_opp_max > 180:
                                        azimuth_opp_max = azimuth_opp_max - 360
                                        if azimuth_opp_min < self.all_Vertsurf[y][0].azimuth <= 180 or -180 <= self.all_Vertsurf[y][0].azimuth < azimuth_opp_max:
                                            self.all_Vertsurf[x][1].extend([[dist,y]])
                                            self.all_Vertsurf[y][1].extend([[dist,x]])
                                    else:
                                        if azimuth_opp_min < self.all_Vertsurf[y][0].azimuth < azimuth_opp_max:
                                            self.all_Vertsurf[x][1].extend([[dist,y]])
                                            self.all_Vertsurf[y][1].extend([[dist,x]])
                        else:
                            azimuth_opp = self.all_Vertsurf[x][0].azimuth - 180
                            azimuth_opp_max = azimuth_opp + toll_az
                            azimuth_opp_min = azimuth_opp - toll_az
                            theta_max = self.all_Vertsurf[x][0].azimuth + toll_theta
                            theta_min = self.all_Vertsurf[x][0].azimuth - toll_theta
                            if theta_max > 180:
                                theta_max = theta_max - 360
                                if theta_min < theta < 180 or -180 <= theta < theta_max:
                                    if azimuth_opp_min < -180:
                                        azimuth_opp_min = azimuth_opp_min + 360
                                        if -180 <= self.all_Vertsurf[y][0].azimuth < azimuth_opp_max or azimuth_opp_min < self.all_Vertsurf[y][0].azimuth <= 180:
                                            self.all_Vertsurf[x][1].extend([[dist,y]])
                                            self.all_Vertsurf[y][1].extend([[dist,x]])
                                    else:
                                        if azimuth_opp_min < self.all_Vertsurf[y][0].azimuth < azimuth_opp_max:
                                            self.all_Vertsurf[x][1].extend([[dist,y]])
                                            self.all_Vertsurf[y][1].extend([[dist,x]])
                            else:
                                if theta_min < theta < theta_max:
                                    if azimuth_opp_min < -180:
                                        azimuth_opp_min = azimuth_opp_min + 360
                                        if -180 <= self.all_Vertsurf[y][0].azimuth < azimuth_opp_max or azimuth_opp_min < self.all_Vertsurf[y][0].azimuth <= 180:
                                            self.all_Vertsurf[x][1].extend([[dist,y]])
                                            self.all_Vertsurf[y][1].extend([[dist,x]])
                                    else:
                                        if azimuth_opp_min < self.all_Vertsurf[y][0].azimuth < azimuth_opp_max:
                                            self.all_Vertsurf[x][1].extend([[dist,y]])
                                            self.all_Vertsurf[y][1].extend([[dist,x]])

        # SECTION 2: Calculation of the shading effect
        for x in range(len(self.all_Vertsurf)):
            if self.all_Vertsurf[x][1] != []:
                self.all_Vertsurf[x][0].OnOff_shading = 'On'
                shading = [0]*len(self.all_Vertsurf[x][1])
                for y in range(len(self.all_Vertsurf[x][1])):
                    
                    # Calculation of the solar height limit
                    if self.all_Vertsurf[x][1][y][0] == 0:
                        sol_h_lim = 90.
                    else:
                        sol_h_lim = np.degrees(np.arctan((self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].surfHeight - self.all_Vertsurf[x][0].centroid_coord[2])/self.all_Vertsurf[x][1][y][0]))
                    self.all_Vertsurf[x][1][y].append(sol_h_lim)
                    
                    # Calculation of the solar azimuth limits
                    sol_az_lim1_xy = np.degrees(np.arccos((self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].vertList[0][0] - self.all_Vertsurf[x][0].centroid_coord[0])/math.sqrt((self.all_Vertsurf[x][0].centroid_coord[0] - self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].vertList[0][0])**2 + (self.all_Vertsurf[x][0].centroid_coord[1] - self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].vertList[0][1])**2)))
                    sol_az_lim2_xy = np.degrees(np.arccos((self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].vertList[1][0] - self.all_Vertsurf[x][0].centroid_coord[0])/math.sqrt((self.all_Vertsurf[x][0].centroid_coord[0] - self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].vertList[1][0])**2 + (self.all_Vertsurf[x][0].centroid_coord[1] - self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].vertList[1][1])**2)))
                    sol_az_lim1 = -(sol_az_lim1_xy + 90)
                    sol_az_lim2 = -(sol_az_lim2_xy + 90)
                    if self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].vertList[0][1] < self.all_Vertsurf[x][0].centroid_coord[1]:
                        sol_az_lim1 = sol_az_lim1 + 2*sol_az_lim1_xy
                    if sol_az_lim1 < -180:
                        sol_az_lim1 = sol_az_lim1 + 360
                    if sol_az_lim1 > 180:
                        sol_az_lim1 = sol_az_lim1 - 360
                    if self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].vertList[1][1] < self.all_Vertsurf[x][0].centroid_coord[1]:
                        sol_az_lim2 = sol_az_lim2 + 2*sol_az_lim2_xy
                    if sol_az_lim2 < -180:
                        sol_az_lim2 = sol_az_lim2 + 360
                    if sol_az_lim2 > 180:
                        sol_az_lim2 = sol_az_lim2 - 360
                    
                    # Necessary conditions:
                    #    1. solar height less than the solar height limit
                    #    2. solar azimuth between the solar azimuth limits
                    
                    shading_sol_h = np.less(Solar_position.elevation,self.all_Vertsurf[x][1][y][2])
                    sol_az_lim_inf = min(sol_az_lim1,sol_az_lim2)
                    sol_az_lim_sup = max(sol_az_lim1,sol_az_lim2)
                    if abs(sol_az_lim_inf - sol_az_lim_sup) < 180:
                        self.all_Vertsurf[x][1][y].append([sol_az_lim_inf,sol_az_lim_sup])
                        shading_sol_az = [np.less(Solar_position.azimuth,sol_az_lim_sup),np.greater(Solar_position.azimuth,sol_az_lim_inf)]
                        shading_tot = shading_sol_h & shading_sol_az[0] & shading_sol_az[1]
                    else:
                        sol_az_lim_inf = max(sol_az_lim1,sol_az_lim2)
                        sol_az_lim_sup = min(sol_az_lim1,sol_az_lim2)
                        self.all_Vertsurf[x][1][y].append([sol_az_lim_inf,sol_az_lim_sup])
                        shading_sol_az1 = [np.less_equal(Solar_position.azimuth,180),np.greater(Solar_position.azimuth,sol_az_lim_inf)]
                        shading_sol_az2 = [np.less(Solar_position.azimuth,sol_az_lim_sup),np.greater_equal(Solar_position.azimuth,-180)]
                        shading_az1 = shading_sol_az1[0] & shading_sol_az1[1]
                        shading_az2 = shading_sol_az2[0] & shading_sol_az2[1]
                        shading_az = shading_az1 | shading_az2
                        shading_tot = shading_sol_h & shading_az
                    shading[y] = shading_tot
                for y in range(len(self.all_Vertsurf[x][1])):
                    if y == 0:
                        shading_eff = shading[y]
                    else:
                        shading_eff = shading_eff | shading[y]
                shading_eff_01 = (1 - np.where(shading_eff==True,1,shading_eff))
                shading_eff_01 = np.where(shading_eff_01==0,R_f,shading_eff_01)
                self.all_Vertsurf[x][0].shading_effect = shading_eff_01

    
    def paramsandloads(self,envelopes,sched_db,weather,mode = 'cityjson'):
        
        '''
        This method permits to calculate the equivalent electrical network
        parameters and thermal loads
        
        Parameters
        ----------
        envelopes : dict
            envelope information for each age class category   
        sched_db : dict
            dictionary containing the operational schedules for each end-use category
        weather : RC_classes.WeatherData.Weather obj
            object of the class weather WeatherData module
        mode : str
            cityjson or geojson mode calculation
        
        Returns
        -------
        None
        
        '''
        
        # Check input data type
        
        if not isinstance(envelopes, dict):
            raise TypeError(f'Ops... JsonCity class - paramsandloads, envelopes is not a dict: envelopes {envelopes}')
        if not isinstance(sched_db, dict):
            raise TypeError(f'Ops... JsonCity class - paramsandloads, sched_db is not a dict: sched_db {sched_db}')
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
        if not isinstance(mode, str):
            raise TypeError(f'Ops... JsonCity class - paramsandloads, mode is not a str: mode {mode}')
        
        # Check input data quality
        
        if not bool(envelopes):
            wrn(f"\n\JsonCity class - paramsandloads, the envelopes dictionary is empty..... envelopes {envelopes}\n")
        if not bool(sched_db):
            wrn(f"\n\JsonCity class - paramsandloads, the envelopes dictionary is empty..... envelopes {envelopes}\n")
        if mode != 'geojson' and  mode != 'cityjson':
            wrn(f"\n\JsonCity class - paramsandloads, the mode doesn't exist..... mode {mode}\n")

        # Parameters and thermal loads calculation 
        if mode == 'cityjson':
            for bd in self.buildings.values():
                self.archId = 1
                bd.BDParamsandLoads(self.model,envelopes,sched_db,weather)
        
        elif mode == 'geojson':
            for i in self.city.index:
                self.buildings[self.city.loc[i]['id']].BDParamsandLoads(self.model,envelopes,sched_db,weather)
    
    
    '''Design Power calculation during design days'''
    def designdays(self,Plant_calc,Time_to_regime,design_days,weather):
        
        '''
        Design Power calculation during design days
        
        Parameters
        ----------
        Plant_calc : str
            YES if plant calculation is performed, NO viceversa
        Time_to_regime : int
            time needed to reach a regime condition for Design Days Calculation
        design_days : list
            period of design days calculation
        weather : RC_classes.WeatherData.Weather obj
            object of the class weather WeatherData module

        Returns
        -------
        None
        
        '''
        
        # Check input data type
        
        if not isinstance(Plant_calc, str):
            raise TypeError(f'Ops... JsonCity class - designdays, Plant_calc is not a str: Plant_calc {Plant_calc}')
        if not isinstance(Time_to_regime, int):
            raise TypeError(f'Ops... JsonCity class - designdays, Time_to_regime is not a int: Time_to_regime {Time_to_regime}')
        if not isinstance(design_days, list):
            raise TypeError(f'Ops... JsonCity class - designdays, design_days is not a list: design_days {design_days}')
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
        
       
        # Calculation in Heating and Cooling seasons
        for bd in self.buildings.values():
            bd.BDdesigndays_Heating(Plant_calc)
                 
        for t in design_days[1]:
            if t == design_days[1][0]:
                x = t
                while x < design_days[1][Time_to_regime - 1]:
                    T_e = weather.Text[t]
                    RH_e = weather.RHext[t]
        
                    if T_e < 0:
                        p_extsat = 610.5*np.exp((21.875*T_e)/(265.5+T_e))
                    else:
                        p_extsat = 610.5*np.exp((17.269*T_e)/(237.3+T_e))
        
                    for bd in self.buildings.values():
                        bd.BDdesigndays_Cooling(t,T_e,RH_e,p_extsat,weather.tau,Plant_calc,self.model)
                        
                    x = x + 1
            else:
                T_e = weather.Text[t]
                RH_e = weather.RHext[t]
            
                if T_e < 0:
                    p_extsat = 610.5*np.exp((21.875*T_e)/(265.5+T_e))
                else:
                    p_extsat = 610.5*np.exp((17.269*T_e)/(237.3+T_e))
            
                for bd in self.buildings.values():
                    bd.BDdesigndays_Cooling(t,T_e,RH_e,p_extsat,weather.tau,Plant_calc,self.model)
        
        # Reset starting values
        for bd in self.buildings.values():
            for z in bd.zones.values():
                z.reset_init_values()
    
    
    ''' Setting plant of each building and checking plant efficiency'''
    def cityplants(self,Plants_list,weather):
        
        '''
        Setting plant of each building and checking plant efficiency
        
        Parameters
        ----------
        Plants_list : dict
            dictionary contaning all the implemented plants
        weather : RC_classes.WeatherData.Weather obj
            object of the class weather WeatherData module
        
        Returns
        -------
        None
        '''
        
        # Check input data type

        if not isinstance(Plants_list, dict):
            raise TypeError(f'Ops... JsonCity class - cityplants, Plants_list is not a dict: Plants_list {Plants_list}')
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')

        # Setting plant
        for bd in self.buildings.values():
            bd.BDplants(Plants_list,weather)
    
    
    '''Energy simulation of the city'''
    def citysim(self,time,weather,Plant_calc,Plants_list):
        
        '''
        This method allows the energy simulation of the city
        
        Parameters
        ----------
        time : simulation timestep
            array of int32
        weather : RC_classes.WeatherData.Weather obj
            object of the class weather WeatherData module
        Plant_calc : str
            YES if plant calculation is performed, NO viceversa
        Plants_list : dict
            dictionary contaning all the implemented plants
        
        Returns
        -------
        None
        '''
        
        # Check input data type
        
        if not isinstance(time, np.ndarray):
            raise TypeError(f'Ops... JsonCity class - citysim, time is not a np.ndarray: time {time}')
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
        if not isinstance(Plant_calc, str):
            raise TypeError(f'Ops... JsonCity class - citysim, Plant_calc is not a str: Plant_calc {Plant_calc}')
        if not isinstance(Plants_list, dict):
            raise TypeError(f'Ops... JsonCity class - citysim, Plants_list is not a dict: Plants_list {Plants_list}')

        # Check input data quality
        
        for t in time:
            if not isinstance(time[t], np.int32):
                wrn(f"\n\nJsonCity class - citysim, at least a component of the vector time is not a np.int32: time[t] {time[t]}\n")
       
        # Energy simulation of the city
        for t in time:
            T_e = weather.Text[t]
            RH_e = weather.RHext[t]
            w = weather.w[t]
            zenith = weather.SolarPosition.zenith[t]
            if T_e < 0:
                p_extsat = 610.5*np.exp((21.875*T_e)/(265.5+T_e))
            else:
                p_extsat = 610.5*np.exp((17.269*T_e)/(237.3+T_e))
                
            # Urban Heat Island evaluation
            if self.urban_canyon_calc:
                radiation = [weather.SolarGains['0.0','0.0','global'].iloc[t],
                             weather.SolarGains['0.0','0.0','direct'].iloc[t]]
                self.urban_canyon.solve_canyon(t,
                                               T_e,
                                               w,
                                               radiation,
                                               self.T_w_0,
                                               T_e-weather.dT_er,
                                               self.H_waste_0*1000,
                                               zenith,
                                               [self.T_out_inf,self.T_out_AHU],
                                               [self.V_0_inf,self.V_0_vent])
                T_e = self.urban_canyon.T_urb[t] - 273.15
                
            # Vectors initialization
            T_w_bd = np.array([],dtype = float)
            T_out_inf_bd = np.array([],dtype = float)
            T_out_vent_bd = np.array([],dtype = float)
            V_0_inf_bd = np.array([],dtype = float)
            V_0_vent_bd = np.array([],dtype = float)
            H_waste_0 = np.array([],dtype = float)
            
            # Linear system resolution
            for bd in self.buildings.values():
                bd.solve(t,
                         T_e,
                         RH_e,
                         p_extsat,
                         weather.tau,
                         Plants_list,
                         Plant_calc,
                         self.model)
                T_w_bd = np.append(T_w_bd, bd.T_wall_0)
                T_out_inf_bd = np.append(T_out_inf_bd, bd.T_out_inf)
                T_out_vent_bd =np.append(T_out_vent_bd, bd.T_out_AHU)
                V_0_inf_bd = np.append(V_0_inf_bd, bd.G_inf_0)
                V_0_vent_bd = np.append(V_0_vent_bd, bd.G_vent_0)
                H_waste_0 = np.append(H_waste_0, bd.H_waste)
                
            # Output post-prcessing in case of Urban Heat Island evaluation
            if self.urban_canyon_calc:
                self.T_w_0 = np.average(T_w_bd, weights = self.bd_ext_walls)
                self.V_0_inf = np.sum(V_0_inf_bd)
                self.V_0_vent = np.sum(V_0_vent_bd)
                self.H_waste_0 = np.sum(H_waste_0)
                try:
                    self.T_out_inf = np.average(T_out_inf_bd, weights = V_0_inf_bd)
                except ZeroDivisionError:
                    self.T_out_inf = 20.
                try:
                    self.T_out_AHU = np.average(T_out_vent_bd, weights = V_0_vent_bd)
                except ZeroDivisionError:
                    self.T_out_AHU = 20.
        
    
    def complexmerge(self):
        
        '''
        Post-processing: output vectors in case of geojson mode.
        ATTENTION: Buildings with the same name must be listed in the
           object city() one after the other
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        '''


        nb = len(self.buildings.values())
        for y in range(nb-1,-1,-1):
            self.complexes[self.city.loc[y]['nome']] = Complex(self.city.loc[y]['nome'])
            self.complexes[self.city.loc[y]['nome']].heatFlowC = self.buildings[self.city.loc[y]['id']].heatFlowBD 
            self.complexes[self.city.loc[y]['nome']].latentFlowC = self.buildings[self.city.loc[y]['id']].latentFlowBD
            self.complexes[self.city.loc[y]['nome']].AHUDemandC = self.buildings[self.city.loc[y]['id']].AHUDemandBD
            self.complexes[self.city.loc[y]['nome']].AHUDemand_latC = self.buildings[self.city.loc[y]['id']].AHUDemand_latBD
            self.complexes[self.city.loc[y]['nome']].AHUDemand_sensC = self.buildings[self.city.loc[y]['id']].AHUDemand_sensBD
        i = 0
        x = 1
        while i < nb:
            if i < nb - 1:
                if self.city.loc[i]['nome'] == self.city.loc[x]['nome']:
                    self.complexes[self.city.loc[i]['nome']].heatFlowC += self.buildings[self.city.loc[x]['id']].heatFlowBD
                    self.complexes[self.city.loc[i]['nome']].latentFlowC += self.buildings[self.city.loc[x]['id']].latentFlowBD
                    self.complexes[self.city.loc[i]['nome']].AHUDemandC += self.buildings[self.city.loc[x]['id']].AHUDemandBD
                    self.complexes[self.city.loc[i]['nome']].AHUDemand_latC += self.buildings[self.city.loc[x]['id']].AHUDemand_latBD
                    self.complexes[self.city.loc[i]['nome']].AHUDemand_sensC += self.buildings[self.city.loc[x]['id']].AHUDemand_sensBD
            i = i + 1
            x = x + 1
             
        
#%%--------------------------------------------------------------------------------------------------- 
#%% some test methods

'''
TEST METHOD
'''
     
if __name__=='__main__':  
    env_path = os.path.join('..','Input','buildings_envelope_V02_EP.xlsx')
    envelopes = loadEnvelopes(env_path)
    '''
    path_p = os.path.join('..','Input', 'ridotto.json')
    padua = JsonCity(path_p,envelopes)
    '''    
    path_v = os.path.join('..','Input', 'verona_geojson.geojson')
    verona = JsonCity(path_v,envelopes,mode='geojson')
    '''
    for i in padua.buildings['RLB-67013'].buildingSurfaces.values():
        if (isinstance(i, Surface) and [726494.291,5032396.487,56.018] in i.vertList and [726495.621,5032407.056,56.018] in i.vertList):
            if i.type=='ExtWall':
                i.printInfo()
    '''
    '''  
    padua.buildings['RLB-67013'].printBuildingInfo()
    padua.buildings['RLB-67001'].printBuildingInfo()
    padua.buildings['RLB-67002'].printBuildingInfo()
    padua.buildings['RLB-67003'].printBuildingInfo()
    padua.buildings['RLB-66818'].printBuildingInfo()
    #[(print(s.name), print(s.area), print(s.vertList))for s in padua.buildings['RLB-67002'].buildingSurfaces.values() if s.type=='ExtWall']
    #padua.buildings['RLB-67013'].plotBuilding(padua.buildings['RLB-67013'].buildingSurfaces['Building Surface 27'])
    #padua.buildings['RLB-67013'].buildingSurfaces['Building Surface 17'].checkSurfaceCoincidence(padua.buildings['RLB-67013'].buildingSurfaces['Building Surface 27'])
    '''        