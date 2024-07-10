# -*- coding: utf-8 -*-
"""
Thermal Storage Tank - Simplified 

The current code offers a simplified approach for a thermal storage tank 
that can be added to the urban building simulation package, EUReCA. 
The calculation is based on a one-node mass-balance and energy-balance
equations for a quasi-steady-state system. 

The main class, thermal storage tank can be initialized using the following
parameters:
    Volume: a float value giving the capacity of the tank in cubic meter
    Area: a float value for the surface area that exchanges heat with ambient,
    it is given in squared meters. if no value is given, then it is calculated 
    assuming an spherical shape and using the volume. 
    Overal_Heat_Transfer_coefficient: a float value describin in W/m2K the heat
    flux that is caused as the difference between the tank and the ambient 
    temperature exist. 
    Specific_heat_capacity and density, given in SI units, they describe the 
    working fluid. 

The tank does not work with heat exchangers, and all the inlets and outlets are
taken and added directly from and to the stored fluid. 

The inlet s and outlets need to be added to the tank using add_demand and 
add_supply methods in which the temperature and mass flow rate (in kg/s) of 
different linkages to the tank are given. 

---
Created by: Mohamad
Betalab - DII, University of Padua
---
"""



import numpy as np


from eureca_building.config import CONFIG
from eureca_building.weather import WeatherFile


class ThermalStorageTank_TES01:
    def __init__(
            self,
            Volume:float, #m3
            Area:float, #m2
            Overal_Heat_Transfer_coefficient:float, #In W/m2K
            specific_heat_capacity=4.218, #kJ/kgK
            density=1000, #kg/m3
            stratification_temperature_difference=10, #C
            charge_temperature=25, #C
            min_temperature=5,
            max_temperature=95
            ):
        self.Volume=Volume
        try: 
            self.Area=Area
        except TypeError:
            self.Area=4.836 * self.Volume ** (2/3) # assume a spherical shape
        
        self.Heat_transfer_coef=Overal_Heat_Transfer_coefficient*self.Area     # W/K                          
        self.Air_temperature=Weather._epw_hourly_data['temp_air']              # C
        self.Array_length=len(self.Air_temperature)
        self.Thermal_capacity=Volume*density*specific_heat_capacity            # J/K  
        self.specific_heat_capacity=specific_heat_capacity*1000                #J/kgK
        self.Demand=[]
        self.Supply=[]
        self.add_supply(+10000000,10,0)
        self.Temperature_stratification=10 
        self.Initial_temperature=charge_temperature                            #C
        self.Minimum_temperature=min_temperature
        self.Maximum_temperature=max_temperature
        
            
        
        
    def add_demand(self,demand_name,mass_flow_rate,temperature,has_return,Temperature_uptake=10):
        if len(mass_flow_rate)==1:
            mass_flow_rate=np.full(self.Array_length,mass_flow_rate)
        elif (len(mass_flow_rate)!=self.Array_length):
            mass_flow_rate=np.interp(np.linespace(0,1,self.Array_length),
                                  np.linspace(0,1,len(mass_flow_rate)),
                                  mass_flow_rate.flatten()).reshape(-1,1) 
        if len(temperature)==1:
            temperature=np.full(self.Array_length,temperature)
        elif (len(temperature)!=self.Array_length):
            temperature=np.interp(np.linespace(0,1,self.Array_length),
                                  np.linspace(0,1,len(temperature)),
                                  temperature.flatten()).reshape(-1,1)  
        
        if has_return==1:
            self.add_demand_mass(self,-mass_flow_rate,temperature-Temperature_uptake,0)
            
        self.Demand.append({"name":demand_name,"mass":mass_flow_rate,"temperature":temperature})
        
    def add_supply(self,supply_name,mass_flow_rate,temperature,has_return):
        if len(mass_flow_rate)==1:
            mass_flow_rate=np.full(self.Array_length,mass_flow_rate)
        elif (len(mass_flow_rate)!=self.Array_length):
            mass_flow_rate=np.interp(np.linespace(0,1,self.Array_length),
                                  np.linspace(0,1,len(mass_flow_rate)),
                                  mass_flow_rate.flatten()).reshape(-1,1) 
        if len(temperature)==1:
            temperature=np.full(self.Array_length,temperature)
        elif (len(temperature)!=self.Array_length):
            temperature=np.interp(np.linespace(0,1,self.Array_length),
                                  np.linspace(0,1,len(temperature)),
                                  temperature.flatten()).reshape(-1,1) 
        self.Supply.append({"name":supply_name,"max_mass":mass_flow_rate,"temperature":temperature,"has_return":has_return,"median_temperature":np.median(temperature)})
        
        
    def _mass_balance(self):
        outtake_mass=np.sum([item['mass'] for item in self.Demand], axis=0)
        inlet_without_return=[item for item in self.Supply if item["has_return"]==0]
        self.Supply = {k: v for k, v in self.Supply.items() if ["has_return"]==0}
        inlet_without_return = sorted(inlet_without_return, 
                                      key="median_temperature",
                                      reverse=True)
        mass_not_taken_care=outtake_mass
        for count, inlet in inlet_without_return:
            inlet["mass"]=np.clip(outtake_mass, a_min=0, a_max=inlet["max_mass"])
            mass_not_taken_care=mass_not_taken_care-inlet["mass"]
        self.Supply.append(inlet_without_return)
        
    def _energy_balance(self):
        demand_energy=np.sum([self.specific_heat_capacity*item['mass']*item['temperature'] \
                              for item in self.Demand], axis=0)
        supply_energy=np.sum([self.specific_heat_capacity*item['mass']*item['temperature'] \
                              for item in self.Supply], axis=0)
        supply_return_outtake_coefficient=np.sum([self.specific_heat_capacity*item['mass']\
                                                         for item in self.Supply\
                                                         if item["has_return"]==1], axis =0)
        fluid_heat_flow=supply_energy-demand_energy
        heat_flow_equivalent_temperature=fluid_heat_flow/self.Heat_transfer_coef
        timestep_seconds=3600/(CONFIG.ts_per_hour)
        thermal_diffusivity=(self.Heat_transfer_coef+supply_return_outtake_coefficient)\
                            /self.Thermal_capacity
        Equivalent_exchanging_temperature=self.Air_temperature+heat_flow_equivalent_temperature
        self.Temperature=np.full(self.Array_length,self.Initial_temperature)
        for i in range(self.Array_length):
            if i==0:
                pass
            else:
                Temperature=self.Temperature[i-1]+timestep_seconds*thermal_diffusivity*\
                                (Equivalent_exchanging_temperature[i-1]-self.Temperature[i-1])
                self.Temperature[i]=min(max(Temperature,self.Minimum_temperature),
                                        self.Maximum_temperature)
    
    def Solve_TES(self,iternum=10):
        for iteration in range(iternum):
            self._mass_balance()
            self._energy_balance()
            for demand_array in self.Demand:
                demand_array["mass"][demand_array["temperature"] > self.Temperature]=0

    
        
        