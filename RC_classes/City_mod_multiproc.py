'''IMPORTING MODULES'''

import sys
import math
import os
import geopandas as gpd
import numpy as np
import time as tm
from cjio import cityjson
from RC_classes.thermalZone import Building, Complex
from RC_classes.Envelope import loadEnvelopes
from RC_classes.Geometry import Surface
from RC_classes.Climate import UrbanCanyon
import time
import concurrent.futures

#%% Useful fubctions


''' City object creation if mode == cityjson''' 
def cityobj(p):  
    with open(p, 'r') as f:
        c = cityjson.CityJSON(file=f)
        return c

''' Buildings and surfaces creation if mode == cityjson'''    
def createBuilding(name,bd,vertList,mode,n_Floors,azSubdiv,hSubdiv,envelopes,sim_time):
    age = bd['attributes']['Age']
    use = bd['attributes']['Use']
    try: 
        H_plant = bd['attributes']['H_Plant']
    except KeyError:
        plant = 'IdealLoad'
    try:
        C_plant = bd['attributes']['C_Plant']
    except KeyError:
        plant = 'IdealLoad'
    
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
                    
        return Building(name,surfaces,use,mode,n_Floors,azSubdiv,hSubdiv,envelopes,age,1,1,H_plant,C_plant,sim_time)


#%% JsonCity class

'''
This class manages city via json input file: geojson or cityjson.

The method __init__ creates the city from the json file. It takes:
    json: path of the json file
    azSubdiv, hSubdiv: the subdivision of the sky dome in sectors (azimuth and height)
    envelopes: envelopes object containing all the envelopes' archetypes
    mode: string with the type of input json file(cityjson, geojson)
    
shading_effect method creates shading vectors for all the distric surfaces.
    .
    .
    .
    .
    .
    
printinfo prints some variable of the district. No input needed.

paramsandloads calculates the RC params and internal gains:
    .
    .
    .
    .
    
citysim does the simulation:
    .
    .
    .
    .

complexmerge. No input required.

'''

class JsonCity():
    
    T_w_0 = 15.
    T_out_inf = 15.
    T_out_AHU = 15.
    V_0_inf = 0.
    V_0_vent = 0.
    H_waste_0 = 0.
    
    
    def __init__(self,json_path,model,azSubdiv,hSubdiv,envelopes,sim_time,mode='cityjson'):
        
        self.model = model
        
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
                
                self.buildings[bd]=createBuilding(bd,self.jsonBuildings[bd],self.city.j['vertices'],self.mode,self.n_Floors,azSubdiv,hSubdiv,envelopes,sim_time)
            
            
        elif mode=='geojson':
            self.mode = mode
            self.city = gpd.read_file(json_path)
            self.jsonBuildings= {}
            self.buildings = {}
            self.complexes = {}
            '''Creates a list of surfaces starting from the footprint points and extruding the footprint'''
            for i in self.city.index:
                self.jsonBuildings[self.city.loc[i]['id']] = self.city.loc[i].to_dict()
                # https://gis.stackexchange.com/questions/287306/list-all-polygon-vertices-coordinates-using-geopandas
                g = [i for i in self.city.loc[i].geometry]
                x,y = g[0].exterior.coords.xy
                coords = np.dstack((x,y)).tolist()
                coords=coords[0]
                coords.pop()
                '''Use coords to create building surfaces list'''
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

                self.buildings[self.city.loc[i]['id']]=Building(self.city.loc[i]['nome'],build_surf,
                                                                self.city.loc[i]['Uso'],
                                                                self.mode,
                                                                self.city.loc[i]['n_piani'],
                                                                azSubdiv,hSubdiv,envelopes,
                                                                self.city.loc[i]['età'],
                                                                self.rh_net,
                                                                self.rh_gross,
                                                                self.heating_plant,
                                                                self.cooling_plant,
                                                                sim_time)
                
        else:
            sys.exit('Set a proper mode')


    def printInfo(self):
        for i in self.buildings.values():
            i.printBuildingInfo()
    
    def create_urban_canyon(self,sim_time,calc,data):
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
   
    def shading_effect(self,mode,toll_az,toll_dist,toll_theta,Solar_position,R_f):
        
        'This method takes into account the shading effect between buildings surfaces'
        
        'SECTION 1: All surfaces are compared and potentially shading surfaces are stored'
        start = tm.time()
        self.all_Vertsurf = []
        
        if mode == 'cityjson':
            for bd in self.buildings.keys():
                self.all_Vertsurf.extend(self.buildings[str(bd)].Vertsurf)
        if mode == 'geojson':
            for i in self.city.index:
                self.all_Vertsurf.extend(self.buildings[i+1].Vertsurf)
                
        end = tm.time()
        print('Shading effect Part1: ',end - start)
        
        start = tm.time()        
        'Each surface is compared with all the others'
        for x in range(len(self.all_Vertsurf)):
            for y in range(len(self.all_Vertsurf)):
                if y > x:
                    'Calculation of the distance between the centroids of the two surfaces under examination'
                    dist = math.sqrt((self.all_Vertsurf[x][0].centroid_coord[0]-self.all_Vertsurf[y][0].centroid_coord[0])**2+(self.all_Vertsurf[x][0].centroid_coord[1]-self.all_Vertsurf[y][0].centroid_coord[1])**2)
                    if dist == 0.0:
                        pass
                    else:
                        'Calculation of the vector direction between the centroids of the two surfaces under examination'
                        theta_xy = np.degrees(np.arccos((self.all_Vertsurf[y][0].centroid_coord[0]-self.all_Vertsurf[x][0].centroid_coord[0])/dist))
                        theta = -(theta_xy + 90)
                        if self.all_Vertsurf[y][0].centroid_coord[1] < self.all_Vertsurf[x][0].centroid_coord[1]:
                            theta = theta + 2*theta_xy
                        if theta < -180:
                            theta = theta + 360
                        if theta > 180:
                            theta = theta - 360
                    '''
                    Conditions:
                        1. the distance between surfaces must be less than toll_dist
                        2. the theta angle between surfaces must be within the range
                        3. the azimuth angle of the second surface must be within the range
                    '''
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
        end = tm.time()
        print('Shading effect Part2: ',end - start)
        
        start = tm.time()
        
        'SECTION 2: Calculation of the shading effect'
        for x in range(len(self.all_Vertsurf)):
            if self.all_Vertsurf[x][1] != []:
                self.all_Vertsurf[x][0].OnOff_shading = 'On'
                shading = [0]*len(self.all_Vertsurf[x][1])
                for y in range(len(self.all_Vertsurf[x][1])):
                    'Calculation of the solar height limit'
                    if self.all_Vertsurf[x][1][y][0] == 0:
                        sol_h_lim = 90.
                    else:
                        sol_h_lim = np.degrees(np.arctan((self.all_Vertsurf[self.all_Vertsurf[x][1][y][1]][0].surfHeight - self.all_Vertsurf[x][0].centroid_coord[2])/self.all_Vertsurf[x][1][y][0]))
                    self.all_Vertsurf[x][1][y].append(sol_h_lim)
                    'Calculation of the solar azimuth limits'
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
                    '''
                    Necessary conditions:
                        1. solar height less than the solar height limit
                        2. solar azimuth between the solar azimuth limits
                    '''
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
        
        end = tm.time()
        print('Shading effect Part3: ',end - start)
    
    def paramsandloads(self,envelopes,sched_db,Solar_gain,w,T_ext,dT_er,mode):
        
        if mode == 'cityjson':
            for bd in self.buildings.values():
                self.archId = 1
                bd.BDParamsandLoads(self.model,envelopes,sched_db,Solar_gain,w,T_ext,dT_er)
        
        elif mode == 'geojson':
            for i in self.city.index:
                self.buildings[self.city.loc[i]['id']].BDParamsandLoads(self.model,envelopes,sched_db,Solar_gain,w,T_ext,dT_er)
    
    
    '''Design Power calculation during design days'''
    def designdays(self,Plant_calc,Time_to_regime,design_days,T_ext,RH_ext,tau):
        
        '''Calculation in Heating and Cooling seasons'''

        for bd in self.buildings.values():
            bd.BDdesigndays_Heating(Plant_calc)
                 
        for t in design_days[1]:
            if t == design_days[1][0]:
                x = t
                while x < design_days[1][Time_to_regime - 1]:
                    T_e = T_ext[t]
                    RH_e = RH_ext[t]
        
                    if T_e < 0:
                        p_extsat = 610.5*np.exp((21.875*T_e)/(265.5+T_e))
                    else:
                        p_extsat = 610.5*np.exp((17.269*T_e)/(237.3+T_e))
        
                    for bd in self.buildings.values():
                        bd.BDdesigndays_Cooling(t,T_e,RH_e,p_extsat,tau,Plant_calc,self.model)
                        
                    x = x + 1
            else:
                T_e = T_ext[t]
                RH_e = RH_ext[t]
            
                if T_e < 0:
                    p_extsat = 610.5*np.exp((21.875*T_e)/(265.5+T_e))
                else:
                    p_extsat = 610.5*np.exp((17.269*T_e)/(237.3+T_e))
            
                for bd in self.buildings.values():
                
                    bd.BDdesigndays_Cooling(t,T_e,RH_e,p_extsat,tau,Plant_calc,self.model)
        
        '''Reset starting values'''
        for bd in self.buildings.values():
            for z in bd.zones.values():
                z.reset_init_values()
    
    
    ''' Setting plant of each building and checking plant efficiency'''
    def cityplants(self,Plants_list,T_ext_H_avg):
        
        for bd in self.buildings.values():
            bd.BDplants(Plants_list,T_ext_H_avg)
    
    
    '''Energy simulation of the city'''
    def citysim(self,time,T_ext,RH_ext,w,Solar_Gains,Solar_position,T_sky,tau,Plant_calc,Plants_list):
        
        for t in time:
            T_e = T_ext[t]
            RH_e = RH_ext[t]
            
            if T_e < 0:
                p_extsat = 610.5*np.exp((21.875*T_e)/(265.5+T_e))
            else:
                p_extsat = 610.5*np.exp((17.269*T_e)/(237.3+T_e))
            
            if self.urban_canyon_calc:
                radiation = [Solar_Gains['0.0','0.0','global'].iloc[t],
                             Solar_Gains['0.0','0.0','direct'].iloc[t]]
                self.urban_canyon.solve_canyon(t,T_e,w[t],radiation,self.T_w_0,T_sky[t],self.H_waste_0*1000,Solar_position.zenith[t],[self.T_out_inf,self.T_out_AHU],[self.V_0_inf,self.V_0_vent])
                T_e = self.urban_canyon.T_urb[t] - 273.15
            
            
            
            T_w_bd = np.array([],dtype = float)
            T_out_inf_bd = np.array([],dtype = float)
            T_out_vent_bd = np.array([],dtype = float)
            V_0_inf_bd = np.array([],dtype = float)
            V_0_vent_bd = np.array([],dtype = float)
            H_waste_0 = np.array([],dtype = float)
            '''
            for bd in self.buildings.values():
                bd.solve(t,Plants_list,T_e,RH_e,p_extsat,tau,Plant_calc,self.model)
                T_w_bd = np.append(T_w_bd, bd.T_wall_0)
                T_out_inf_bd = np.append(T_out_inf_bd, bd.T_out_inf)
                T_out_vent_bd =np.append(T_out_vent_bd, bd.T_out_AHU)
                V_0_inf_bd = np.append(V_0_inf_bd, bd.G_inf_0)
                V_0_vent_bd = np.append(V_0_vent_bd, bd.G_vent_0)
                H_waste_0 = np.append(H_waste_0, bd.H_waste)
            '''
            self.bd_parallel_list = [[bd,t,Plants_list,T_e,RH_e,p_extsat,tau,Plant_calc,self.model] for bd in self.buildings.values()]

            with concurrent.futures.ThreadPoolExecutor() as executor:
                bd_executor = executor.map(bd_parallel_solve, self.bd_parallel_list)
            
                        
            
            if self.urban_canyon_calc:
                for bd in self.buildings.values():
                    #bd.solve(t,Plants_list,T_e,RH_e,p_extsat,tau,Plant_calc,self.model)
                    T_w_bd = np.append(T_w_bd, bd.T_wall_0)
                    T_out_inf_bd = np.append(T_out_inf_bd, bd.T_out_inf)
                    T_out_vent_bd =np.append(T_out_vent_bd, bd.T_out_AHU)
                    V_0_inf_bd = np.append(V_0_inf_bd, bd.G_inf_0)
                    V_0_vent_bd = np.append(V_0_vent_bd, bd.G_vent_0)
                    H_waste_0 = np.append(H_waste_0, bd.H_waste)
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
        
    
    '''Post-processing: output vectors '''
    def complexmerge(self):
        
        '''ATTENTION: Buildings with the same name must be listed in the
           object city() one after the other'''
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
                    
#%%   

def bd_parallel_solve(x):
    bd,t,Plants_list,T_e,RH_e,p_extsat,tau,Plant_calc,model = x
    bd.solve(t,Plants_list,T_e,RH_e,p_extsat,tau,Plant_calc,model)


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