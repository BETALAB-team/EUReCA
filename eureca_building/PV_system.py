"""
PV System calculations:

Input variables from the system:
    Building:
        Thermal Zones:
            Surfaces lists
            each surface with tilt, azimuth, surface_type, and shading_coefficient
    Coverage Factor:
        the ratio of the surface that is used for PV panels
    Electrical Load
        
    
PV parameters:
    A simplified method developed by A. Driesse ^1 ^2 is used which needs 5 
    parameters describing the efficiency of the PV unit
    
Battery parameters:
    The battery has its capacity, voltage, charge and discharge efficiency
    and maximum and minimum charge as input parameters.
    
Weatherfile:
    outdoor air temperature
    direct irradiance
    global irradiance
    solar azimuth 
    solar height


1
A. Driesse and J. S. Stein, “From IEC 61853 power measurements to PV system simulations”, Sandia Report No. SAND2020-3877, 2020. DOI: 10.2172/1615179

2
A. Driesse, M. Theristis and J. S. Stein, “A New Photovoltaic Module Efficiency Model for Energy Prediction and Rating,” in IEEE Journal of Photovoltaics, vol. 11, no. 2, pp. 527-534, March 2021. DOI: 10.1109/JPHOTOV.2020.3045677


---
Created by: Mohamad
Betalab - DII, University of Padua
---

"""



import os



import pandas as pd
import numpy as np

import pvlib
from pvlib import pvsystem as pvs
## https://pvlib-python.readthedocs.io

from eureca_building.weather import WeatherFile
from eureca_building.config import CONFIG

        
class PV_system():
    '''
    initializing the PV system for each thermal zone. 
    '''
    
    def __init__(self,
                 name: str,
                 epw_path: str,
                 surface_list: list,
                 mount_surfaces=["Roof"],
                 coverage_factor=0.8):
        self.epw=epw_path
        self.name=name
        self.coverage_factor=coverage_factor
        self._surfaces=[s for s in surface_list if s.surface_type in mount_surfaces]
        self.weather,self.weather_md=pvlib.iotools.read_epw(self.epw)
        self.pv_data_adr()
        self.pv_efficiencies()
        self.pv_data_install()
        self.Battery()
    

    
    def pv_efficiencies(self):
        '''
        calculates the efficiency of the photovoltaic cells at each timestep
        currently it uses ADR method

        '''
        Weather=self.weather
        Wdf = pd.DataFrame({'ghi': Weather['ghi'], 'dhi': Weather['dhi'], 'dni': Weather['dni'],
                           'temp_air': Weather['temp_air'],
                           'wind_speed': Weather['wind_speed']
                           })
        Wdf.index=Wdf.index-pd.Timedelta(minutes=30)
        loc=pvlib.location.Location.from_epw(self.weather_md)
        # loc = Weather._site
        solpos=loc.get_solarposition(Wdf.index)
        efficiencies={}
        for s in self._surfaces:
            tilt=s._height_round
            orient=s._azimuth_round
            
            total_irrad=pvlib.irradiance.get_total_irradiance(tilt,orient,solpos.apparent_zenith,
                                                              solpos.azimuth,Wdf.dni,Wdf.ghi,Wdf.dhi)
            Wdf['poa_global']=total_irrad.poa_global
            Wdf['temp_pv']=pvlib.temperature.faiman(Wdf.poa_global,Wdf.temp_air,Wdf.wind_speed)
            Wdf['rel_eta']=pvlib.pvarray.pvefficiency_adr(Wdf['poa_global'],Wdf['temp_pv'],**self.adr_parameters)
            efficiencies[s]=Wdf
        self._pv_efficiencies=efficiencies
        
        
    
    def pv_data_adr(self):
        '''
        setting the parameters needed for the adr model

        '''
        self.adr_parameters={'k_a':0.99924,
                        'k_d':-5.49097,
                        'tc_d':0.01918,
                        'k_rs':0.06999,
                        'k_rsh':0.26144}
    def pv_data_install(self):
        '''
        setting the PV capacity modeling data:
            model power and coverage factor
        

        '''
        pv_p_mod_stc=100
        pv_A_mod=1.7
        for k,v in self._pv_efficiencies.items():
            A=k._area*self.coverage_factor
            self._pv_efficiencies[k]['power_stc']=A/pv_A_mod*pv_p_mod_stc
        self.pv_g_stc=1000 #W/m2
           
       
    def pv_production(self):
        '''
        Based on the efficiencies calculated in pv_efficiencies
        and the nominal and install parameters given in pv_data_install,
        the production is calculated for each timestep

        '''
        # ProductionWh=self.weather['ghi']
        for k,v in self._pv_efficiencies.items():
            PVDF=v
            Production_watt=PVDF['power_stc']*PVDF['rel_eta']*PVDF['poa_global']/self.pv_g_stc
            ProductionWh=Production_watt/CONFIG.ts_per_hour
            
        ProductionWh.index=ProductionWh.index+pd.Timedelta(minutes=30)
        
        return ProductionWh
    
    def Battery(self):
        '''
        sets the battery performance parameters. 

        '''
        self.battery_parameters={'charge_efficiency':0.98,
                                 'discharge_efficiency':0.98,
                                 'max_charge':0.95,
                                 'min_charge':0.2,
                                 'voltage':12,
                                 'capacity':1500}
        
    def Battery_charge(self,electricity,pv_prod):
        '''
        calculates the battery charge at each timestep. using it to calculate the energy going into the battery,
        and taken from battery. therefore it is also capable of calculating the grid 
        electricity.
        returns:
            state: battery charge level in [Ah]
            tobattery: energy stored in batteries at each timestep [Wh]
            frombattery: energy taken from the batteries at each timestep [Wh]
            togrid: energy giveb to grid at each timestep [Wh]
            fromgrid: energy taken from grid at each timestep [Wh]
            directsolar: energy consumed directly from solar production at each timestep [Wh]

        '''

        pv_prod=pv_prod.to_numpy()

        electricity=electricity.to_numpy()

        inout=pv_prod-electricity
        charger=inout.copy()
        discharger=inout.copy()
        charger[charger<0]=0
        discharger[discharger>0]=0
        state=np.zeros_like(discharger)
        statediff=np.zeros_like(state)
        maxcharge=self.battery_parameters['max_charge']*self.battery_parameters['voltage']*self.battery_parameters['capacity']
        mincharge=self.battery_parameters['min_charge']*self.battery_parameters['voltage']*self.battery_parameters['capacity']
        statechange=charger*self.battery_parameters['charge_efficiency']+discharger*self.battery_parameters['discharge_efficiency']
        
        for j, elem in np.ndenumerate(state):
            i=j[0]
            if i==0:
                state[i]=state[i]
                statediff[i]=0
            else:

                state[i]=state[i-1]+statechange[i-1]
                if state[i] > maxcharge:
                    state[i]=maxcharge
                if state[i]<mincharge:
                    if state[i]<0:
                        state[i]=0
                    if state[i-1]>=mincharge:
                        state[i]=mincharge
                statediff[i]=state[i]-state[i-1]
        
        state=state/self.battery_parameters['voltage']
        tobattery=np.minimum(statechange,statediff)
        tobattery[tobattery<0]=0
        frombattery=np.maximum(statechange,statediff)
        frombattery[frombattery>0]=0 
        frombattery=np.abs(frombattery)
        togrid=charger-tobattery
        fromgrid=np.abs(discharger)-frombattery
        directsolar=np.minimum(pv_prod,electricity)

        
        

        return [state , tobattery, frombattery, togrid, fromgrid, directsolar]
                
        
        
   