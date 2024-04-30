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



import pandas as pd
import numpy as np

import pvlib
from eureca_building.weather import WeatherFile
from eureca_building.config import CONFIG
#class
class PV_system():
    '''
    initializing the PV system for each thermal zone. 
    '''
    
    def __init__(self,
                 name: str,
                 weatherobject: WeatherFile,
                 surface_list: list,
                 mount_surfaces=["Roof"],
                 coverage_factor=0.8): #Coverage factor for the area that is covered by the PV
        self.name=name
        self.coverage_factor=coverage_factor
        self._surfaces=[s for s in surface_list if s.surface_type in mount_surfaces]
        self.weather=weatherobject._epw_hourly_data
        self.weather_md=weatherobject._epw_general_data
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
        Weather_dataframe = pd.DataFrame({'ghi': Weather['ghi'], 'dhi': Weather['dhi'], 'dni': Weather['dni'],
                           'temp_air': Weather['temp_air'],
                           'wind_speed': Weather['wind_speed']
                           })
        Weather_dataframe.index=Weather_dataframe.index-pd.Timedelta(minutes=30)
        LocationAttributes=pvlib.location.Location.from_epw(self.weather_md)
        # LocationAttributes = Weather._site
        Solar_position=LocationAttributes.get_solarposition(Weather_dataframe.index)
        efficiencies={}
        for s in self._surfaces:
            tilt=s._height_round
            orient=s._azimuth_round
            
            total_irrad=pvlib.irradiance.get_total_irradiance(tilt,orient,Solar_position.apparent_zenith,
                                                              Solar_position.azimuth,Weather_dataframe.dni,Weather_dataframe.ghi,Weather_dataframe.dhi)
            Weather_dataframe['poa_global']=total_irrad.poa_global
            Weather_dataframe['temp_pv']=pvlib.temperature.faiman(Weather_dataframe.poa_global,Weather_dataframe.temp_air,Weather_dataframe.wind_speed)
            Weather_dataframe['rel_eta']=pvlib.pvarray.pvefficiency_adr(Weather_dataframe['poa_global'],Weather_dataframe['temp_pv'],**self.adr_parameters)
            efficiencies[s]=Weather_dataframe
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
        Module_Power_STC_condition=100 #Module Power at STC
        Single_Module_Surface_Area=1.7 #Module Area
        for Surface,v in self._pv_efficiencies.items():
            Installed_PV_area=Surface._area*self.coverage_factor
            self._pv_efficiencies[Surface]['power_stc']=Installed_PV_area/Single_Module_Surface_Area*Module_Power_STC_condition
        self.pv_g_stc=1000 #W/m2
           
       
    def pv_production(self):
        '''
        Based on the efficiencies calculated in pv_efficiencies
        and the nominal and install parameters given in pv_data_install,
        the production is calculated for each timestep

        '''
        # ProductionWh=self.weather['ghi']
        for Surface,Dataframe in self._pv_efficiencies.items():
            Photovoltaic_Generation_DataFrame=Dataframe
            Production_watt=Photovoltaic_Generation_DataFrame['power_stc']*Photovoltaic_Generation_DataFrame['rel_eta']*Photovoltaic_Generation_DataFrame['poa_global']/self.pv_g_stc
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
            State: battery charge level in [Ah]
            tobattery: energy stored in batteries at each timestep [Wh]
            frombattery: energy taken from the batteries at each timestep [Wh]
            togrid: energy giveb to grid at each timestep [Wh]
            fromgrid: energy taken from grid at each timestep [Wh]
            directsolar: energy consumed directly from solar production at each timestep [Wh]

        '''

        pv_prod=pv_prod.to_numpy()

        electricity=electricity.to_numpy()

        Generation_Consumption_Balance=pv_prod-electricity
        charger=Generation_Consumption_Balance.copy()
        discharger=Generation_Consumption_Balance.copy()
        charger[charger<0]=0
        discharger[discharger>0]=0
        State=np.zeros_like(discharger)
        State_Change=np.zeros_like(State)
        maxcharge=self.battery_parameters['max_charge']*self.battery_parameters['voltage']*self.battery_parameters['capacity']
        mincharge=self.battery_parameters['min_charge']*self.battery_parameters['voltage']*self.battery_parameters['capacity']
        Momentual_Charge_Dischare=charger*self.battery_parameters['charge_efficiency']+discharger*self.battery_parameters['discharge_efficiency']
        
        for j, elem in np.ndenumerate(State):
            i=j[0]
            if i==0:
                State[i]=State[i]
                State_Change[i]=0
            else:

                State[i]=State[i-1]+Momentual_Charge_Dischare[i-1]
                if State[i] > maxcharge:
                    State[i]=maxcharge
                if State[i]<mincharge:
                    if State[i]<0:
                        State[i]=0
                    if State[i-1]>=mincharge:
                        State[i]=mincharge
                State_Change[i]=State[i]-State[i-1]
        
        State=State/self.battery_parameters['voltage']
        tobattery=np.minimum(Momentual_Charge_Dischare,State_Change)
        tobattery[tobattery<0]=0
        frombattery=np.maximum(Momentual_Charge_Dischare,State_Change)
        frombattery[frombattery>0]=0 
        frombattery=np.abs(frombattery)
        togrid=charger-tobattery
        fromgrid=np.abs(discharger)-frombattery
        directsolar=np.minimum(pv_prod,electricity)
        State_percent=State/self.battery_parameters['capacity']
        
        

        return [State_percent , tobattery, frombattery, togrid, fromgrid, directsolar]
                
        
        
   