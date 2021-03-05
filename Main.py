#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 12:20:07 2020

@author: enrico
"""

''' IMPORTING MODULES '''

import sys
#import deepcopy as dpc
import pandas as pd
import os
import pvlib
import numpy as np
import time as tm
from RC_classes.IrradiancePreProcessor import PlanesIrradiances
from RC_classes.Envelope import loadEnvelopes
from RC_classes.EndUse import loadArchetype, loadSimpleArchetype
from RC_classes.CityJSON import JsonCity
from RC_classes.WeatherData import TskyCalc, SolarPosition, rescale_weather, rescale_sol_gain
from RC_classes.BuildingsPlants import loadPlants

#%%
'''Setting Input Data'''
# PREPROCESSING
iopath = os.path.join('.','Input', 'ITA_Venezia-Tessera.161050_IGDG.epw')         # Epw directory from home folder
year = 2020                                                                    # Setting an year
first_day = 7                                                                  # Setting the first day of the year (1 Monday)
tz='Europe/Rome'                                                               # Setting the thermal zone
azSubdiv = 8                                                                   # Azimuth subdivision
hSubdiv = 3                                                                    # Tilt subdivision
SolarCalc = False                                                              # Need of a Solar Calculation?
env_path = os.path.join('.','Input','Envelopes.xlsx')                          # Envelope directory from home folder
schedpath = os.path.join('.','Input','ScheduleSemp.xlsx')                     # Annual schedules
SchedMethod = str('A')                                                         # A Daily schedules, B yearly schedules
plant_path = os.path.join('.','Input','PlantsList.xlsx')                       # Plants list and properties

# SIMULATION
ts = 1                                                                         # Timestep per hour
hours = 8760                                                                   # Hours of simulation
years = 1                                                                      # Years of simulation
model = str('2C')                                                              # Select model: '1C' or '2C'
mode = str('cityjson')                                                         # Select mode: 'geojson' or 'cityjson'
jsonfile = str('PaduaRestricted.json')                      # Select .geojson or .json file
path=os.path.join('.','Input', jsonfile)
Shading_calc = str('NO')                                                       # Select 'YES' or 'NO' to take into consideration the shading effect
toll_az = float(80)                                                            # Semi-tollerance on azimuth of shading surface [°]
toll_dist = float(100)                                                         # Tollerance on distance of shading surface [m]
toll_theta = float(80)                                                          # Semi-tollerance on position of the shading surfaces [°]
R_f = float(0)                                                                 # Reduction factor of the direct solar radiation due to the shading effect [0-1]
DD_boundaries = np.array([[167,504],[4681,5017]], dtype = int)                 # Heating and Cooling Design Days Periods
Time_to_regime = 168                                                           # Time needed to reach a regime condition for Design Days Calculation
DesignDays_calc = str('YES')                                                   # Select 'YES' or 'NO' to calculate or not design days demand
Plant_calc = str('YES')                                                        # Select 'YES' or 'NO' to calculate or not buildings plant

# OUTPUT REPORT
OutRep = bool(False)

# Urban Canyon Data
UWG_calc = bool(True) 
UWG_data = {
    'Area' : float(538000),
    'Diameter': float(800),
    'Perimeter': float(3280),
    'Ortogonal_length': float(800),
    'h_layer_inf': 2,
    'h_layer_sup': 150,
    'z_i_day': 1000,
    'z_i_night': 50,
    'z_i_m': 10,
    'z_i_ref': 150,
    'ExtWall_radiative_coef': float(5),
    'ExtWall_convective_coef': float(17),
    'ExtWall_emissivity': float(0.85),
    'ExtWall_reflection': float(0.15),
    'Average_U_windows': float(2),
    'Vegetation_density': float(0.16),
    'Road_radiative_coef': float(5),
    
    'Road_layers': np.array([
        [0.05, 1.9*10**6, 273 + 10],
        [0.2 , 2*10**6,   273 + 12],
        [1   , 1.4*10**6, 273 + 13],
        ], dtype = float),                          # first column thickness [m], second capacity [J/m3K], first row asphalt, second stones, third gravel 
        
    'Soil_layers': np.array([
        [3 , 1.4*10**6, 273 + 10],
        [3 , 1.4*10**6, 273 + 12],
        [3 , 1.4*10**6, 273 + 13],
        ], dtype = float),                          # first column thickness [m], second capacity [J/m3K]
    'T_deep': float(13+273),
    
    
    'cars_number_over_1000_px': float(39545000/60360),
    'motorcycles_number_over_1000_px': float(6896000/60360),
    'othervehicle_number_over_1000_px': float(5364000/60360),
    'emission_factor_cars_over_1000_px': float(	24.74 ),
    'emission_factor_motorcycles_over_1000_px': float(13.34),
    'emission_factor_othervehicle_over_1000_px': float(104.95),
    'vehicle_multiplying_factor': float(1),
    'distance_per_hour': float(64000),
    'population_density': float(2283.2),
    'first_day': first_day,
    'Hw_week': np.array([0.9,0.5,0.9,2.1,3.3,4.2,5.8,7.1,6.3,5.8,5.4,5.4,5.8,6.3,6.7,7.1,7.5,5.8,4.2,3.3,2.5,1.7,1.3,0.9])/100,
    'Hw_end': np.array([1.7,1.2,1.,1.7,2.6,3.3,4.2,5,5.8,6.2,6.7,6.7,6.7,6.7,6.3,5.8,5,4.6,4.2,3.7,3.3,2.9,2.5,2.1])/100    
    }

'''
Reference H_traffic
Italian vehicles number: http://www.opv.aci.it/WEBDMCircolante/
emission factors
https://link.springer.com/article/10.1007/s00704-008-0086-5/tables/4
with 64 km/h

'''


#%%
'''PRE-PROCESSING'''

start = tm.time()



'''Loading Envelope and Schedule Data'''
envelopes = loadEnvelopes(env_path)                                            # Envelope file loading

if SchedMethod == 'A':
    PlantDays=[2520,3984,6192,6912]                                            # 15th April, 15th June, 15th September, 15th October
    sched = loadSimpleArchetype(schedpath,np.arange(8760),first_day,ts,PlantDays)
elif SchedMethod == 'B':
    sched = loadArchetype(schedpath,np.arange(8760,ts))
else:
    sys.exit('Set a proper schedule inporting methodology')

'''Plants List'''
Plants_list = loadPlants(plant_path)

end = tm.time()
print('Pre-processing:       ', end - start)

#%% 
''' SIMULATION '''
path=os.path.join('.','Input', jsonfile)
sim_time = np.arange(len(T_ext)) 
if DesignDays_calc == 'NO' and Plant_calc == 'YES':
    print('------------------------------------------------------------------------------'+
          '\nWARNING: To do plant calculation Design Days calculation must be done!!\nSimulation will run doing the calculation of the design days\n'+
          '------------------------------------------------------------------------------')
    DesignDays_calc = 'YES'
    
'District initialization'
start = tm.time()
Padua = JsonCity(path,model,azSubdiv,hSubdiv,envelopes,[ts,hours],mode)
Padua.create_urban_canyon([ts,hours],UWG_calc,UWG_data)
end = tm.time()
print('Jsoncity:             ', end - start)

'Mutual shading effect evaluation'
start = tm.time()
if Shading_calc == 'YES':
    Padua.shading_effect(mode,toll_az,toll_dist,toll_theta,Solar_position,R_f)
elif Shading_calc == 'NO':
    pass
else:
    sys.exit('Select an allowed value of Shading_calc')
end = tm.time()
print('Shading effect TOT:   ', end - start)

'Parameters and loads calculation'
start = tm.time()
Padua.paramsandloads(envelopes,sched,Solar_Gains,w,T_ext,dT_er,mode)
end = tm.time()
print('Paramscalc:           ', end - start)

'Design days and Buildings Plants evaluation'
start = tm.time()
design_days = [np.arange(DD_boundaries[0,0]*ts,DD_boundaries[0,1]*ts),np.arange(DD_boundaries[1,0]*ts,DD_boundaries[1,1]*ts)] 
if DesignDays_calc == 'YES':
    Padua.designdays(Plant_calc,Time_to_regime,design_days,T_ext,RH_ext,3600/ts)
elif DesignDays_calc == 'NO':
    pass
else:
    sys.exit('Select an allowed value of DesignDays_calc')

end = tm.time()
print('DesignDaysCalc:       ', end - start)
   
start = tm.time()
if Plant_calc == 'YES':
    Padua.cityplants(Plants_list,T_ext_H_avg)
else:
    pass
end = tm.time()
print('Building Plant:       ', end - start)

'Simulation'
start = tm.time()
Padua.citysim(sim_time,T_ext,RH_ext,w,Solar_Gains,Solar_position,(T_ext - dT_er),3600/ts,Plant_calc,Plants_list)
end = tm.time()
print('Simulation:           ', end - start)

#%%
''' POST-PROCESSING '''

start = tm.time()
year = (8760-1)*ts+1
nb = len(Padua.buildings.values())

Col_Bui = ['Name', 'Age class', 'End use', 'H plant', 'C plant',
            'BuiHeight','Footprint','ExtWallArea','nFloor','TotalArea','Volume','PDesH','PDesC',
            'TotOpaquaA','TotWinA','TotUA',
            'ExtWallU','RoofU','GroundU','WinU',
            'Htr_is','Htr_w','Htr_ms','Htr_em','Cm',
            'RrestAW','R1AW','RalphaStarAW','RalphsStarIL','RaplphaStarIW','R1IW','C1AW','C1IW']
BuiInfo = pd.DataFrame(index = Padua.buildings.keys() ,columns = Col_Bui)

for bd in Padua.buildings.keys():
    
    building = Padua.buildings[bd]
    z = building.zones['Zone']
    
    data = [building.name, building.archId, building.end_use, building.H_plant_type, building.C_plant_type,
            building.buildingHeight, building.footprint, building.extWallArea, building.nFloors, building.total_area, building.Volume,
            building.Pnom_H_BD, building.Pnom_C_BD,
      z.Tot_opaque_area, z.Tot_glazed_area, z.UA_tot,
      z.Strat['ExtWall'].U, z.Strat['Roof'].U, z.Strat['GroundFloor'].U, z.Strat['Window'].U,
      z.Htr_is, z.Htr_w, z.Htr_ms, z.Htr_em, z.Cm,
      z.RrestAW,z.R1AW,z.RalphaStarAW,z.RalphaStarIL,z.RalphaStarIW,z.R1IW,z.C1AW,z.C1IW]
    
    data_dic = dict(zip(Col_Bui,data))
    BuiInfo.loc[bd] = data_dic

'''Output vectors inizialization'''
HF = np.zeros([nb,year])
LF = np.zeros([nb,year])
AHUD = np.zeros([nb,year])
T = np.zeros([nb,year])
RH = np.zeros([nb,year])
ElEn = np.zeros([nb,year])
GasCon = np.zeros([nb,year])
Final = np.zeros([nb,year])

i=0
columns = []
for bd in Padua.buildings.keys():
    columns.append(bd)
    HF[i]=Padua.buildings[bd].heatFlowBD
    LF[i]=Padua.buildings[bd].latentFlowBD
    AHUD[i]=Padua.buildings[bd].AHUDemandBD
    T[i] = Padua.buildings[bd].zones['Zone'].Air_temp
    RH[i] = Padua.buildings[bd].zones['Zone'].RH_i
    ElEn[i] = Padua.buildings[bd].BDPlant.Electrical_energy_consumption
    GasCon[i] = Padua.buildings[bd].BDPlant.Gas_consumption
    Final[i] = Padua.buildings[bd].BDPlant.Final_energy_demand   
    i += 1
    
H_F = HF.transpose()
L_F = LF.transpose()
AHU_D = AHUD.transpose()
T = T.transpose()
RH = RH.transpose()
ElEn = ElEn.transpose()
GasCon = GasCon.transpose()
Final = Final.transpose()

'''Output'''
HeatFlow = H_F
LatentFlow = L_F
AHUDemand = AHU_D
Temp = T
RelHum = RH
ElectricalEnergy = ElEn
GasConsumption = GasCon
FinalEnergy = Final

if OutRep:
    pd.DataFrame(data = H_F, columns = columns).to_csv(os.path.join('OutputReport',model,'HeatFlow.csv'))
    pd.DataFrame(data = L_F, columns = columns).to_csv(os.path.join('OutputReport',model,'LantentFlow.csv'))
    pd.DataFrame(data = AHU_D, columns = columns).to_csv(os.path.join('OutputReport',model,'AHU.csv'))
    pd.DataFrame(data = T, columns = columns).to_csv(os.path.join('OutputReport',model,'Temp.csv'))
    pd.DataFrame(data = RH, columns = columns).to_csv(os.path.join('OutputReport',model,'RH.csv'))
    pd.DataFrame(data = ElEn, columns = columns).to_csv(os.path.join('OutputReport',model,'ElectricalEnergy.csv'))
    pd.DataFrame(data = GasCon, columns = columns).to_csv(os.path.join('OutputReport',model,'GasConsumption.csv'))
    # pd.DataFrame(data = Final, columns = columns).to_csv(os.path.join('OutputReport',model,'FinalEnergy.csv'))
    
    BuiInfo.to_csv(os.path.join('OutputReport',model,'BuildingsParams.csv'))

if mode == 'geojson':
    Padua.complexmerge()

end = tm.time()
print('Post-processing:      ', end - start)
