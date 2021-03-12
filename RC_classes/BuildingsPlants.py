'''IMPORTING MODULES'''

import pandas as pd
import numpy as np
from RC_classes.auxiliary_functions import wrn
from RC_classes.WeatherData import  Weather

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions

def loadPlants(path):
    '''
    loading Plants List
    
    Parameters
    ----------
    path : string
        Path containing the string of the PlantsList.xlsx

    Returns
    -------
    Typologies: dictionary with plant_key/Plants data    
    '''

    # Check input data type
    
    if not isinstance(path, str):
        raise TypeError(f'ERROR input path is not a string: path {path}') 
    
    try:
        Plant_list = pd.read_excel(path,sheet_name="List",header=0,index_col=[0])
    except FileNotFoundError:
        raise FileNotFoundError(f'ERROR Failed to open the plants xlsx file {path}... Insert a proper path')
        
    Typologies = dict()
    for i in Plant_list.index:
        if i != 'Plant_type':
            Typologies[i] = pd.read_excel(path,sheet_name=i,header=1)
    return Typologies

#%%--------------------------------------------------------------------------------------------------- 
#%% Plants class

class Plants():

    '''
    Class Plants provides all plants with their own properties
    This class is intended to be a container of the real plant object (see the following classes)
    It manages the interface (input output functions and attributes
    
    Methods:
        init
        setPlant
        solvePlant
    '''
    
    def __init__(self,H_plant_type,C_plant_type,l,ts):
        '''
        Initializes an Plant object
        
        Parameters
            ----------
            H_plant_type : string
                string that sets the heating plant
            C_plant_type : string
                string that sets the heating plant
            l : int
                number of simulation time steps [-]
            ts : int
                Number of time steps per hour.
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        if not isinstance(H_plant_type, str):
            raise TypeError(f'ERROR Plant inizialization, H_plant_type must be a string: H_plant_type {H_plant_type}')
        if not isinstance(C_plant_type, str):
            raise TypeError(f'ERROR Plant inizialization, C_plant_type must be a string: C_plant_type {C_plant_type}')
        if not isinstance(l, int):
            raise TypeError(f'ERROR Plant inizialization, l must be an integer: l {l}')
        if not isinstance(ts, int):
            raise TypeError(f'ERROR Plant inizialization, input ts is not an integer: ts {ts}') 
    
        # Output vectors initialization 
        
        self.ts = ts
        self.H_plant_type = H_plant_type
        self.C_plant_type = C_plant_type
        self.Total_energy_losses = np.zeros(l)                                 # [kW]
        self.Final_energy_demand = np.zeros(l)                               # [kW]
        self.Electrical_energy_consumption = np.zeros(l)                       # [kW]
        self.Gas_consumption = np.zeros(l)                                     # [Nm3]
        self.H_waste = np.zeros(l)
        
        # Heating plant 
        
        try:
            if H_plant_type == 'CondensingBoiler':
                self.heating_Plant = CondensingBoiler(H_plant_type,l)
            elif H_plant_type == 'TraditionalBoiler':
                self.heating_Plant = TraditionalBoiler(H_plant_type,l)
            elif H_plant_type == 'IdealLoad':
                self.heating_Plant = IdealLoad(H_plant_type,l)
        except KeyError:
            raise KeyError(f"ERROR Heating plant type not found: H_plant_type {H_plant_type}. Set a proper one")
            
        # Cooling plant 
        
        try:
            if C_plant_type == 'IdealLoad':
                self.cooling_Plant = IdealLoad(C_plant_type,l)
            elif C_plant_type == 'SplitAirCooler':
                self.cooling_Plant = SplitAirCooler(C_plant_type,l)
            elif C_plant_type == 'ChillerAirtoWater':
                self.cooling_Plant = ChillerAirtoWater(C_plant_type,l)
            elif C_plant_type == 'SplitAirConditioner':
                self.cooling_Plant = SplitAirConditioner(C_plant_type,l)
        except KeyError:
            raise KeyError(f"ERROR Cooling plant type not found: C_plant_type {C_plant_type}. Set a proper one")
    
    def setPlant(self,Typologies,weather,Pnom_H = 1E10,Pnom_C = -1E10):
        '''
        Sets plant design power
        
        Parameters
            ----------
            Typologies: dictionary
                dictionary with plant_key/Plants data
            weather : RC_classes.WeatherData.Weather obj
                object of the class weather WeatherData module
            Pnom_H : float
                Design Heating power  [kW]
            Pnom_C : float
                Design Cooling power  [kW]
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type

        if not isinstance(Typologies, dict):
            raise TypeError(f'ERROR Plant setPlant, Typologies must be a dictionary: Typologies {Typologies}')
        if not isinstance(Pnom_H, float):
            raise TypeError(f'ERROR Plant setPlant, Pnom_H must be a float: Pnom_H {Pnom_H}')
        if not isinstance(Pnom_C, float):
            raise TypeError(f'ERROR Plant setPlant, Pnom_C must be a float: Pnom_C {Pnom_C}')
        if not isinstance(weather, Weather):
            raise TypeError(f'ERROR Plant setPlant, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
        
        # Check data quality
        
        if Pnom_H > 1E9 or Pnom_C < -1E9:
            wrn(f"WARNING There' a building plant with a design power outside limits [-1E9:1E9]: P_non_H {Pnom_H}, P_nom_C {Pnom_C}")
        
        # Heating plant 
        
        if self.H_plant_type == 'IdealLoad':
            pass
        elif self.H_plant_type == 'CondensingBoiler':
            self.heating_Plant.setCondensingBoiler(Typologies,Pnom_H,weather.THextAv)
        elif self.H_plant_type == 'TraditionalBoiler':
            self.heating_Plant.setTraditionalBoiler(Typologies,Pnom_H,weather.THextAv)
        
        # Cooling plant 
        
        if self.C_plant_type == 'IdealLoad':
            pass
        elif self.C_plant_type == 'SplitAirCooler':
            self.cooling_Plant.setSplitAirCooler(Typologies,Pnom_C)
        elif self.C_plant_type == 'ChillerAirtoWater':
            self.cooling_Plant.setChillerAirtoWater(Typologies,Pnom_C)
        elif self.C_plant_type == 'SplitAirConditioner':
            self.cooling_Plant.setSplitAirConditioner(Typologies,Pnom_C)
    
    def solvePlant(self,Typologies,heatFlow,t,T_ext,T_int,RH_i):
        '''
        Solve plant object for a time step
        
        Parameters
            ----------
            Typologies: dictionary
                dictionary with plant_key/Plants data
            heatFlow : float
                required power  [kW]
            t : int
                Simulation time step
            T_ext : float
                External temperature
            T_int : float
                Internal temperature
            RH_i : float
                Relative Humidity
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        if not isinstance(Typologies, dict):
            raise TypeError(f'ERROR Plant solvePlant, time step {t}, Typologies must be a dictionary: Typologies {Typologies}')
        if not isinstance(t, int):            
            raise TypeError(f'ERROR Plant solvePlant, time step {t}, t must be an integer: t {t}')
        if not isinstance(heatFlow, float):
            raise TypeError(f'ERROR Plant solvePlant, time step {t}, heatFlow must be a float: heatFlow {heatFlow}')
        if not isinstance(T_ext, float):
            raise TypeError(f'ERROR Plant solvePlant, time step {t}, T_ext must be a float: T_ext {T_ext}')
        if not isinstance(T_int, float):
            raise TypeError(f'ERROR Plant solvePlant, time step {t}, T_int must be a float: T_int {T_int}')
        if not isinstance(RH_i, float):
            raise TypeError(f'ERROR Plant solvePlant, time step {t}, RH_i must be a float: RH_i {RH_i}')
               
        # Heating plant 
        
        if self.H_plant_type == 'CondensingBoiler':
            if heatFlow > 0:
                self.heating_Plant.solveCondensingBoiler(Typologies,heatFlow,t)
                self.Total_energy_losses[t] = self.heating_Plant.phi_gn_i_Px[t]
                self.Final_energy_demand[t] = heatFlow + self.Total_energy_losses[t]
                self.Electrical_energy_consumption[t] = self.heating_Plant.W_aux_Px[t]
                self.Gas_consumption[t] = self.Final_energy_demand[t]/self.heating_Plant.PCI_natural_gas/self.ts
            else:
                pass  
        
        elif self.H_plant_type == 'TraditionalBoiler':
            if heatFlow > 0:
                self.heating_Plant.solveTraditionalBoiler(Typologies,heatFlow,t)
                self.Total_energy_losses[t] = self.heating_Plant.phi_gn_i_Px[t]
                self.Final_energy_demand[t] = heatFlow + self.Total_energy_losses[t]
                self.Electrical_energy_consumption[t] = self.heating_Plant.W_aux_Px[t]
                self.Gas_consumption[t] = self.Final_energy_demand[t]/self.heating_Plant.PCI_natural_gas/self.ts
            else:
                pass 
        
        elif self.H_plant_type == 'IdealLoad':
            self.heating_Plant.solveIdealLoad(heatFlow,t)
        
        # Cooling plant 
        
        if self.C_plant_type == 'SplitAirCooler':
            if heatFlow < 0:
                self.cooling_Plant.solveSplitAirCooler(Typologies,heatFlow,t,T_ext,T_int,RH_i)
                self.Electrical_energy_consumption[t] = (self.cooling_Plant.W_el[t] + self.cooling_Plant.W_aux[t])
                self.H_waste[t] = self.cooling_Plant.H_waste[t]
            else:
                pass
        
        elif self.C_plant_type == 'ChillerAirtoWater':
            if heatFlow < 0:
                self.cooling_Plant.solveChillerAirtoWater(Typologies,heatFlow,t,T_ext)
                self.Electrical_energy_consumption[t] = (self.cooling_Plant.W_el[t] + self.cooling_Plant.W_aux[t])
                self.H_waste[t] = self.cooling_Plant.H_waste[t]
            else:
                pass
        
        elif self.C_plant_type == 'IdealLoad':
            self.cooling_Plant.solveIdealLoad(heatFlow,t)
            
        elif self.C_plant_type == "SplitAirConditioner":
            if heatFlow < 0:
                self.Plant.solveSplitAirConditioner(Typologies,heatFlow,t,T_ext)
                self.Electrical_energy_consumption[t] = self.Plant.W_el[t] 
                self.H_waste[t] = self.Plant.H_waste[t]
            else:
                pass

       

#%%--------------------------------------------------------------------------------------------------- 
#%% IdealLoad class

class IdealLoad:

    '''
    Class IdealLoad is for the ideal  zone balance. This actualy passes all methods.
    
    Methods:
        init
        solvePlant
    '''

    def __init__(self,plant_type,l):
        '''
        init. This method does nothing and passes
        
        Parameters
            ----------
            H_plant_type : string
                string that sets the heating plant
            l : int
                number of simulation time steps [-]
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
    
        pass
        
    def solveIdealLoad(self,heatFlow,t):
        '''
        solveIdealLoad. This method does nothing and passes
        This class does not fill the energy vectors arrays in the Plant class, therefore
        in case of ideal load, only heatFlow array, latentFlow array and AHU array will be 
        filled with the thermal zone sensible, latent and ventilation power.
        
        Parameters
            ----------
            heatFlow : float
                required power  [kW]
            t : int
                Simulation time step
                
        Returns
        -------
        None.
        
        '''  
        
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
    
        pass


#%%--------------------------------------------------------------------------------------------------- 
#%% CondensingBoiler class
class CondensingBoiler:
    
    '''
    Class CondensingBoiler
    This method considers a generic Condensing Boiler as the heating plant 
    of the entire building following UNI-TS 11300:2 - 2008
    
    Methods:
        init
        setCondensingBoiler
        solvePlant
    '''
    
    def __init__(self,plant_type,l):
        '''
        init. Set some attributes for the method
        
        Parameters
            ----------
            plant_type : string
                string that sets the heating plant
            l : int
                number of simulation time steps [-]
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
    
        # Input Data
        self.H_plant_type = plant_type
        self.theta_gn_test_Pn = 70                                             # [°C]
        self.theta_gn_test_Pint = 35                                           # [°C]
        self.theta_a_test = 20                                                 # [°C]
        self.PCI_natural_gas = 9.97                                            # [kWh/Nm3]
        
        # Vectors initialization
        self.phi_gn_i_Px = np.zeros(l)                                         # [kW]
        self.W_aux_Px = np.zeros(l)                                            # [kW]
        self.phi_gn_i_Px_1 = np.zeros(l)                                       # [kW]
        self.phi_gn_i_Px_2 = np.zeros(l)                                       # [kW]
        
    def setCondensingBoiler(self,plant_info,Pnom,T_ext_H_avg):
        '''
        Sets plant design power
        
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            Pnom : float
                Design Heating power  [kW]
            T_ext_H_avg : float
                Average external temperature during the heating season.
                
        Returns
        -------
        None.
        '''        
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
        
        # Choise of plant size based on estimated nominal Power 
        Size_0 = 0
        self.Pnom = Pnom                                                       # [kW]
        for i in plant_info[self.H_plant_type].index:
            if Size_0 < self.Pnom <= plant_info[self.H_plant_type]['Size'][i]:
                self.ID = i
            Size_0 = plant_info[self.H_plant_type]['Size'][i]
        if self.Pnom > plant_info[self.H_plant_type]['Size'][plant_info[self.H_plant_type].index[-1]]:
            self.ID = len(plant_info[self.H_plant_type]) - 1
        self.Pint = plant_info[self.H_plant_type]['P_int'][self.ID]*self.Pnom    # [kW]
        if plant_info[self.H_plant_type]['location'][self.ID] == 'tech_room':
            self.theta_a_gn = 15                                               # [°C]
        elif plant_info[self.H_plant_type]['location'][self.ID] == 'internal':
            self.theta_a_gn = 20                                               # [°C]
        elif plant_info[self.H_plant_type]['location'][self.ID] == 'external':
            self.theta_a_gn = T_ext_H_avg                                      # [°C]
        
        # Minimum Condensing Boiler efficiency check based on Directive 
        # 92/42/CEE'''
        if self.Pnom > 400:
            self.Pnom_check = 400                                              # [kW]
        else:
            self.Pnom_check = self.Pnom
        self.eta_gn_Pn_min = (91 + 1*np.log(self.Pnom_check))                  # [%]
        self.eta_gn_Pint_min = (97 + 1*np.log(self.Pnom_check))                # [%]
        if plant_info[self.H_plant_type]['eta_nom'][self.ID] < self.eta_gn_Pn_min or plant_info[self.H_plant_type]['eta_int'][self.ID] < self.eta_gn_Pint_min:
            print('Warning: Condensing Boiler efficiencies are lower than minimum')
        
        # No load losses estimation 
        self.phi_gn_I_P0 = (self.Pnom_check*10*4.8*self.Pnom_check**(-0.35))/1000  # [kW]
        
        # Auxiliary power input data
        self.W_aux_Pn = (45*self.Pnom**0.48)/1000                                  # [kW]
        self.W_aux_Pint = (15*self.Pnom**0.48)/1000                                # [kW]
        self.W_aux_P0 = (15)/1000                                                  # [kW]
        self.FC_Pint = self.Pint/self.Pnom
        
    def solveCondensingBoiler(self,plant_info,heatFlow,t):
        '''
        This method allows to calculate Condensing Boiler losses following
        the Standard UNI-TS 11300:2 - 2008
        
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            heatFlow : float
                required power  [kW]
            t : int
                Simulation time step
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
        
        # Corrected efficiency and losses at nominal power 
        self.eta_gn_Pn_cor = plant_info[self.H_plant_type]['eta_nom'][self.ID] + plant_info[self.H_plant_type]['f_cor_Pn'][self.ID]*(self.theta_gn_test_Pn - plant_info[self.H_plant_type]['T_gn_w'][self.ID])
        self.phi_gn_i_Pn_cor = (100 - self.eta_gn_Pn_cor)/self.eta_gn_Pn_cor*self.Pnom                                           # [kW]
        
        # Corrected efficiency and losses at intermadiate power 
        self.eta_gn_Pint_cor = plant_info[self.H_plant_type]['eta_int'][self.ID] + plant_info[self.H_plant_type]['f_cor_Pint'][self.ID]*(self.theta_gn_test_Pint - plant_info[self.H_plant_type]['T_gn_Pint'][self.ID])
        self.phi_gn_i_Pint_cor = (100 - self.eta_gn_Pint_cor)/self.eta_gn_Pint_cor*self.Pint                                     # [kW]
        
        # Corrected no-load losses 
        self.phi_gn_i_P0_cor = self.phi_gn_I_P0*((plant_info[self.H_plant_type]['T_gn_w'][self.ID] - self.theta_a_gn)/(self.theta_gn_test_Pn - self.theta_a_test))**1.25
        
        # Total losses at current power 
        if heatFlow <= 0:
            self.phi_gn_i_Px[t] = 0
        else:
            self.phi_gn_i_Px_1[t] = (heatFlow**2)*((self.Pint*(self.phi_gn_i_Pn_cor - self.phi_gn_i_P0_cor) - self.Pnom*(self.phi_gn_i_Pint_cor - self.phi_gn_i_P0_cor))/(self.Pnom*self.Pint*(self.Pnom - self.Pint)))
            self.phi_gn_i_Px_2[t] = heatFlow*(((self.Pnom**2)*(self.phi_gn_i_Pint_cor - self.phi_gn_i_P0_cor) - (self.Pint**2)*(self.phi_gn_i_Pn_cor - self.phi_gn_i_P0_cor))/(self.Pnom*self.Pint*(self.Pnom - self.Pint)))
            self.phi_gn_i_Px[t] = self.phi_gn_i_Px_1[t] + self.phi_gn_i_Px_2[t] + self.phi_gn_i_P0_cor                                               # [kW]
        
        # Auxiliary power estimation 
        if heatFlow <= 0:
            self.FC_Px = 0
            self.W_aux_Px[t] = 0
        else:
            self.FC_Px = heatFlow/self.Pnom
            if 0 < self.FC_Px <= self.FC_Pint:
                self.W_aux_Px[t] = self.W_aux_P0 + self.FC_Px/self.FC_Pint*(self.W_aux_Pint - self.W_aux_P0)                         # [kW]
            elif self.FC_Pint < self.FC_Px <= 1:
                self.W_aux_Px[t] = self.W_aux_Pint + (self.FC_Px - self.FC_Pint)*(self.W_aux_Pn - self.W_aux_Pint)/(1-self.FC_Pint)  # [kW]
            elif self.FC_Px > 1:
                self.FC_Px = 1
                self.W_aux_Px[t] = self.W_aux_Pint + (self.FC_Px - self.FC_Pint)*(self.W_aux_Pn - self.W_aux_Pint)/(1-self.FC_Pint)  # [kW]
        

#%%--------------------------------------------------------------------------------------------------- 
#%% TraditionalBoiler class
class TraditionalBoiler:
    
    '''
    Class TraditionalBoiler
    This method considers a generic Traditional Boiler as the heating plant 
    of the entire building following UNI-TS 11300:2 - 2008    
        
    Methods:
        init
        setCondensingBoiler
        solvePlant
    '''
    
    def __init__(self,plant_type,l):
        '''
        init. Set some attributes for the method
        
        Parameters
            ----------
            plant_type : string
                string that sets the heating plant
            l : int
                number of simulation time steps [-]
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary        
        
        # Input Data 
        self.H_plant_type = plant_type
        self.theta_gn_test_Pn = 70                                             # [°C]
        self.theta_gn_test_Pint = 50                                           # [°C]
        self.theta_a_test = 20                                                 # [°C]
        self.PCI_natural_gas = 9.97                                            # [kWh/Nm3]
        
        # Vectors initialization 
        self.phi_gn_i_Px = np.zeros(l)                                         # [kW]
        self.W_aux_Px = np.zeros(l)                                            # [kW]
        self.phi_gn_i_Px_1 = np.zeros(l)
        self.phi_gn_i_Px_2 = np.zeros(l)
        
    def setTraditionalBoiler(self,plant_info,Pnom,T_ext_H_avg):
        '''
        Sets plant design power
        
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            Pnom : float
                Design Heating power  [kW]
            T_ext_H_avg : float
                Average external temperature during the heating season.
                
        Returns
        -------
        None.
        '''        
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
        
        # Choise of plant size based on estimated nominal Power 
        Size_0 = 0
        self.Pnom = Pnom                                                       # [kW]
        for i in plant_info[self.H_plant_type].index:
            if Size_0 < self.Pnom <= plant_info[self.H_plant_type]['Size'][i]:
                self.ID = i
            Size_0 = plant_info[self.H_plant_type]['Size'][i]
        if self.Pnom > plant_info[self.H_plant_type]['Size'][plant_info[self.H_plant_type].index[-1]]:
            self.ID = len(plant_info[self.H_plant_type]) - 1
        self.Pint = plant_info[self.H_plant_type]['P_int'][self.ID]*self.Pnom    # [kW]
        if plant_info[self.H_plant_type]['location'][self.ID] == 'tech_room':
            self.theta_a_gn = 15                                               # [°C]
        elif plant_info[self.H_plant_type]['location'][self.ID] == 'internal':
            self.theta_a_gn = 20                                               # [°C]
        elif plant_info[self.H_plant_type]['location'][self.ID] == 'external':
            self.theta_a_gn = T_ext_H_avg                                      # [°C]
        
        # Minimum Traditional Boiler efficiency check based on Directive 
        #    92/42/CEE 
        if self.Pnom > 400:
            self.Pnom_check = 400                                              # [kW]
        else:
            self.Pnom_check = self.Pnom
        self.eta_gn_Pn_min = (84 + 2*np.log(self.Pnom_check))                  # [%]
        self.eta_gn_Pint_min = (80 + 3*np.log(self.Pnom_check))                # [%]
        if plant_info[self.H_plant_type]['eta_nom'][self.ID] < self.eta_gn_Pn_min or plant_info[self.H_plant_type]['eta_int'][self.ID] < self.eta_gn_Pint_min:
            print('Warning: Traditional Boiler efficiencies are lower than minimum')
        
        # No load losses estimation 
        self.phi_gn_I_P0 = (self.Pnom_check*10*8.5*self.Pnom_check**(-0.4))/1000   # [kW]
        
        # Auxiliary power input data 
        self.W_aux_Pn = (45*self.Pnom**0.48)/1000                                  # [kW]
        self.W_aux_Pint = (15*self.Pnom**0.48)/1000                                # [kW]
        self.W_aux_P0 = (15)/1000                                                  # [kW]
        self.FC_Pint = self.Pint/self.Pnom
        
    def solveTraditionalBoiler(self,plant_info,heatFlow,t):
        '''
        This method allows to calculate Traditional Boiler losses following
        the Standard UNI-TS 11300:2 - 2008
        
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            heatFlow : float
                required power  [kW]
            t : int
                Simulation time step
                
        Returns
        -------
        None.
        '''  
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
        
        # Corrected efficiency and losses at nominal power 
        self.eta_gn_Pn_cor = plant_info[self.H_plant_type]['eta_nom'][self.ID] + plant_info[self.H_plant_type]['f_cor_Pn'][self.ID]*(self.theta_gn_test_Pn - plant_info[self.H_plant_type]['T_gn_w'][self.ID])
        self.phi_gn_i_Pn_cor = (100 - self.eta_gn_Pn_cor)/self.eta_gn_Pn_cor*self.Pnom                                           # [kW]
        
        # Corrected efficiency and losses at intermadiate power 
        self.eta_gn_Pint_cor = plant_info[self.H_plant_type]['eta_int'][self.ID] + plant_info[self.H_plant_type]['f_cor_Pint'][self.ID]*(self.theta_gn_test_Pint - plant_info[self.H_plant_type]['T_gn_Pint'][self.ID])
        self.phi_gn_i_Pint_cor = (100 - self.eta_gn_Pint_cor)/self.eta_gn_Pint_cor*self.Pint                                     # [kW]
        
        # Corrected no-load losses 
        self.phi_gn_i_P0_cor = self.phi_gn_I_P0*((plant_info[self.H_plant_type]['T_gn_w'][self.ID] - self.theta_a_gn)/(self.theta_gn_test_Pn - self.theta_a_test))**1.25
        
        # Total losses at current power 
        if heatFlow <= 0:
            self.phi_gn_i_Px[t] = 0
        else:
            self.phi_gn_i_Px_1[t] = (heatFlow**2)*((self.Pint*(self.phi_gn_i_Pn_cor - self.phi_gn_i_P0_cor) - self.Pnom*(self.phi_gn_i_Pint_cor - self.phi_gn_i_P0_cor))/(self.Pnom*self.Pint*(self.Pnom - self.Pint)))
            self.phi_gn_i_Px_2[t] = heatFlow*(((self.Pnom**2)*(self.phi_gn_i_Pint_cor - self.phi_gn_i_P0_cor) - (self.Pint**2)*(self.phi_gn_i_Pn_cor - self.phi_gn_i_P0_cor))/(self.Pnom*self.Pint*(self.Pnom - self.Pint)))
            self.phi_gn_i_Px[t] = self.phi_gn_i_Px_1[t] + self.phi_gn_i_Px_2[t] + self.phi_gn_i_P0_cor                                               # [kW]
        
        # Auxiliary power estimation 
        if heatFlow <= 0:
            self.FC_Px = 0
            self.W_aux_Px[t] = 0
        else:
            self.FC_Px = heatFlow/self.Pnom
            if 0 < self.FC_Px <= self.FC_Pint:
                self.W_aux_Px[t] = self.W_aux_P0 + self.FC_Px/self.FC_Pint*(self.W_aux_Pint - self.W_aux_P0)                         # [kW]
            elif self.FC_Pint < self.FC_Px <= 1:
                self.W_aux_Px[t] = self.W_aux_Pint + (self.FC_Px - self.FC_Pint)*(self.W_aux_Pn - self.W_aux_Pint)/(1-self.FC_Pint)  # [kW]
            elif self.FC_Px > 1:
                self.FC_Px = 1
                self.W_aux_Px[t] = self.W_aux_Pint + (self.FC_Px - self.FC_Pint)*(self.W_aux_Pn - self.W_aux_Pint)/(1-self.FC_Pint)  # [kW]

#%%--------------------------------------------------------------------------------------------------- 
#%% SplitAirCooler class

class SplitAirCooler:
    '''
    SplitAirCooler class
    
    This method considers a generic Split Air Conditioner as the cooling plant 
    of the entire building following UNI-TS 11300:3 - 2010
        
    Methods:
        init
        setCondensingBoiler
        solvePlant
    '''
    
    def __init__(self,plant_type,l):
        '''
        init. Set some attributes for the method
        
        Parameters
            ----------
            plant_type : string
                string that sets the cooling plant
            l : int
                number of simulation time steps [-]
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary   
        
        # Input Data
        self.l = l                                                             # Timestep
        self.C_plant_type = plant_type
        self.T_ext_rif_100 = 35                                                # [°C]
        self.T_ext_rif_75 = 30                                                 # [°C]
        self.T_ext_rif_50 = 25                                                 # [°C]
        self.T_ext_rif_25 = 20                                                 # [°C]
        self.T_int_rif = 27                                                    # [°C]
        self.W_aux_gn = 0.04                                                   # [kW/kW_cond]
        self.Conversion_fac = 2.17
        
        # Vectors initialization 
        self.W_el = np.zeros(l)                                                # [kWe]
        self.W_aux = np.zeros(l)                                               # [kWe]
        self.H_waste = np.zeros(l)                                             # [kW]
        self.eta_mm = np.zeros(l)
        self.eta_1 = np.zeros(l)
        
    def setSplitAirCooler(self,plant_info,Pnom):
        '''
        Choise of plant size based on estimated nominal Power
        
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            Pnom : float
                Design Cooling power  [kW]
                
        Returns
        -------
        None.
        '''        
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
        
        Size_0 = 0
        self.Pnom = - Pnom                                                     # [kW]
        for i in plant_info[self.C_plant_type].index:
            if Size_0 < self.Pnom <= plant_info[self.C_plant_type]['Size'][i]:
                self.ID = i
            Size_0 = plant_info[self.C_plant_type]['Size'][i]
        if self.Pnom > plant_info[self.C_plant_type]['Size'][plant_info[self.C_plant_type].index[-1]]:
            self.ID = len(plant_info[self.C_plant_type]) - 1
        self.Q_max = Pnom                                                      # [kW]
        
        # EER curve in case of Air-to-Air unit 
        self.EER_100 = plant_info[self.C_plant_type]['EER_100'][self.ID]
        self.EER_75 = plant_info[self.C_plant_type]['EER_75'][self.ID]
        self.EER_50 = plant_info[self.C_plant_type]['EER_50'][self.ID]
        self.EER_25 = plant_info[self.C_plant_type]['EER_25'][self.ID]
        self.EER_20 = self.EER_25*0.94
        self.EER_15 = self.EER_25*0.85
        self.EER_10 = self.EER_25*0.73
        self.EER_5 = self.EER_25*0.50
        self.EER_2 = self.EER_25*0.26
        self.EER_1 = self.EER_25*0.14
        self.x = np.array([1,0.75,0.5,0.25,0.2,0.15,0.1,0.05,0.02,0.01,0.0])
        self.y = np.array([self.EER_100,self.EER_75,self.EER_50,self.EER_25,self.EER_20,self.EER_15,self.EER_10,self.EER_5,self.EER_2,self.EER_1,0.0])

        self.prosp_C_100 = np.array([[1.634, 1.457, 1.281, 1.105, 0.928, 0.807, 0.685],\
                            [1.720, 1.518, 1.327, 1.148, 0.979, 0.856, 0.736],\
                            [1.756, 1.534, 1.332, 1.155, 1.000, 0.871, 0.750],\
                            [1.782, 1.569, 1.370, 1.187, 1.018, 0.890, 0.767],\
                            [1.834, 1.639, 1.444, 1.249, 1.054, 0.928, 0.802]])
        self.prosp_C_75 =  np.array([[1.457, 1.281, 1.105, 0.928, 0.807, 0.685, 0.585],\
                            [1.518, 1.327, 1.148, 0.979, 0.856, 0.736, 0.636],\
                            [1.534, 1.332, 1.155, 1.000, 0.871, 0.750, 0.650],\
                            [1.569, 1.370, 1.187, 1.018, 0.890, 0.767, 0.667],\
                            [1.639, 1.444, 1.249, 1.054, 0.928, 0.802, 0.700]])
        self.prosp_C_50 =  np.array([[1.281, 1.105, 0.928, 0.807, 0.685, 0.585, 0.505],\
                            [1.327, 1.148, 0.979, 0.856, 0.736, 0.636, 0.556],\
                            [1.332, 1.155, 1.000, 0.871, 0.750, 0.650, 0.672],\
                            [1.370, 1.187, 1.018, 0.890, 0.767, 0.667, 0.587],\
                            [1.444, 1.249, 1.054, 0.928, 0.802, 0.700, 0.698]])
        self.prosp_C_25 =  np.array([[1.062, 0.962, 0.871, 0.788, 0.714, 0.646, 0.585],\
                            [1.083, 0.981, 0.888, 0.804, 0.728, 0.659, 0.596],\
                            [1.105, 1.000, 0.905, 0.820, 0.742, 0.672, 0.608],\
                            [1.126, 1.020, 0.923, 0.836, 0.757, 0.685, 0.620],\
                            [1.149, 1.040, 0.941, 0.852, 0.771, 0.698, 0.632]])
        self.prosp_C_0 = self.prosp_C_25
        self.prosp_C = np.array([self.prosp_C_100,self.prosp_C_75,self.prosp_C_50,self.prosp_C_25,self.prosp_C_0])
        self.fx = np.array([1,0.75,0.5,0.25,0])
        self.pcfx = np.array([15,20,25,30,35,40,45])
        self.pcfy = np.array([16,18,19,20,22])
    
    def solveSplitAirCooler(self,plant_info,heatFlow,t,T_ext,T_int,RH_i):
        '''
        This method allows to calculate Split Air Cooler electrical power
        following the Standard UNI-TS 11300:3 - 2010
        
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            heatFlow : float
                required power  [kW]
            t : int
                Simulation time step
            T_ext : float
                External temperature
            T_int : float
                Internal temperature
            RH_i : float
                Relative Humidity
                
        Returns
        -------
        None.        
        '''
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
        
        self.T_ext_solve = T_ext
        if self.T_ext_solve <= 15:
            self.T_ext_solve = 15.1
        elif self.T_ext_solve > 45:
            self.T_ext_solve = 45
        
        self.prosp_C_F = np.zeros([5,7])
        
        # Wet bulb temperature 
        self.T_wb = T_int*np.arctan(0.151977*(RH_i+8.313659)**(0.5))+np.arctan(T_int+RH_i)-np.arctan(RH_i-1.676331)+0.00391838*(RH_i)**(3/2)*np.arctan(0.023101*RH_i)-4.686035
        if self.T_wb <= 16:
            self.T_wb = 16.1
        elif self.T_wb > 22:
            self.T_wb = 22
            
        # Load Factor and EER in real conditions estimation 
        if heatFlow >= 0:
            pass
        else:
            self.F = heatFlow/(self.Q_max)            
            for i in range(len(self.x)-1):
                xa = self.x[i+1]
                xb = self.x[i]
                if xa < self.F <= xb:
                    self.EER = (self.F-xb)/(xa-xb)*self.y[i+1] - (self.F-xa)/(xa-xb)*self.y[i]
                elif self.F > 1:
                    self.EER = self.EER_100
                else:
                    pass
            for i in range(len(self.fx)-1):
                n_prosp_a = i+1
                n_prosp_b = i
                fxa = self.fx[i+1]
                fxb = self.fx[i]
                if fxa < self.F <= fxb:
                    for n in range(5):
                        for m in range(7):
                            self.prosp_C_F[n][m] = (self.F-fxb)/(fxa-fxb)*self.prosp_C[n_prosp_a][n][m] - (self.F-fxa)/(fxa-fxb)*self.prosp_C[n_prosp_b][n][m]
                elif self.F > 1:
                    self.prosp_C_F = self.prosp_C_100
                else:
                    pass
            for n in range(4):
                y1 = self.pcfy[n]
                y2 = self.pcfy[n+1]
                if y1 < self.T_wb <= y2:
                    for m in range(6):
                        x1 = self.pcfx[m]
                        x2 = self.pcfx[m+1]
                        if x1 < self.T_ext_solve <= x2:
                            self.eta_1[t]=1/((x2-x1)*(y2-y1))*(self.prosp_C_F[n+1][m+1]*(x2-self.T_ext_solve)*(y2-self.T_wb) + self.prosp_C_F[n][m+1]*(self.T_ext_solve-x1)*(y2-self.T_wb) + self.prosp_C_F[n+1][m]*(x2-self.T_ext_solve)*(self.T_wb-y1) + self.prosp_C_F[n][m]*(self.T_ext_solve-x1)*(self.T_wb-y1))
            
            # Electrical power required in real conditions 
            self.eta_mm[t] = self.EER*self.eta_1[t]
            self.W_el[t] = abs(heatFlow)/self.eta_mm[t]                        # [kWe]
            self.H_waste[t] = abs(heatFlow)+self.W_el[t]                       # [kW]
            self.W_aux[t] = self.W_aux_gn*(self.H_waste[t])                    # [kWe]

        
#%%--------------------------------------------------------------------------------------------------- 
#%% ChillerAirtoWater class

class ChillerAirtoWater:
    
    '''
    ChillerAirtoWater class
    
    This method considers a generic Air-to-water Chiller as the cooling plant 
    of the entire building following UNI-TS 11300:3 - 2010
    
        
    Methods:
        init
        setCondensingBoiler
        solvePlant
    '''

    
    def __init__(self,plant_type,l):
    
        '''
        init. Set some attributes for the method
        
        Parameters
            ----------
            plant_type : string
                string that sets the cooling plant
            l : int
                number of simulation time steps [-]
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary   
        
        # Input Data 
        self.l = l                                                             # Timestep
        self.C_plant_type = plant_type
        self.T_ext_rif_100 = 35                                                # [°C]
        self.T_ext_rif_75 = 30                                                 # [°C]
        self.T_ext_rif_50 = 25                                                 # [°C]
        self.T_ext_rif_25 = 20                                                 # [°C]
        self.T_w_out_rif = 7                                                   # [°C]
        self.W_aux_gn = 0.04                                                   # [kW/kW_cond]
        self.Conversion_fac = 2.17
        
        # Vectors initialization 
        self.W_el = np.zeros(l)                                                # [kWe]
        self.W_aux = np.zeros(l)                                               # [kWe]
        self.H_waste = np.zeros(l)                                             # [kW]
        self.eta_mm = np.zeros(l)
        self.eta_1 = np.zeros(l)      
    
    def setChillerAirtoWater(self,plant_info,Pnom):
        '''
        Choise of plant size based on estimated nominal Power
        
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            Pnom : float
                Design Cooling power  [kW]
                
        Returns
        -------
        None.
        '''        
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
        
        # Choise of plant size based on estimated nominal Power 
        Size_0 = 0
        self.Pnom = abs(Pnom)                                                  # [kW]
        for i in plant_info[self.C_plant_type].index:
            if Size_0 < self.Pnom <= plant_info[self.C_plant_type]['Size'][i]:
                self.ID = i
            Size_0 = plant_info[self.C_plant_type]['Size'][i]
        if self.Pnom > plant_info[self.C_plant_type]['Size'][plant_info[self.C_plant_type].index[-1]]:
            self.ID = len(plant_info[self.C_plant_type]) - 1
        self.Q_max = Pnom                                                      # [kW]
        
        # EER curve in case of Air-to-Water chiller 
        self.EER_100 = plant_info[self.C_plant_type]['EER_100'][self.ID]
        self.EER_75 = plant_info[self.C_plant_type]['EER_75'][self.ID]
        self.EER_50 = plant_info[self.C_plant_type]['EER_50'][self.ID]
        self.EER_25 = plant_info[self.C_plant_type]['EER_25'][self.ID]
        self.EER_20 = self.EER_25*0.95
        self.EER_15 = self.EER_25*0.94
        self.EER_10 = self.EER_25*0.87
        self.EER_5 = self.EER_25*0.71
        self.EER_2 = self.EER_25*0.46
        self.EER_1 = self.EER_25*0.29
        self.x = np.array([1,0.75,0.5,0.25,0.2,0.15,0.1,0.05,0.02,0.01,0.0])
        self.y = np.array([self.EER_100,self.EER_75,self.EER_50,self.EER_25,self.EER_20,self.EER_15,self.EER_10,self.EER_5,self.EER_2,self.EER_1,0.0])
        
        # Assumption: T_water_out = 7°C and delta(T) = 5°C 
        self.prosp_C_100 = np.array([1.756, 1.534, 1.332, 1.155, 1.000, 0.871, 0.750])
        self.prosp_C_75 =  np.array([1.534, 1.332, 1.155, 1.000, 0.871, 0.750, 0.650])
        self.prosp_C_50 =  np.array([1.332, 1.155, 1.000, 0.871, 0.750, 0.650, 0.570])
        self.prosp_C_25 =  np.array([1.155, 1.000, 0.871, 0.750, 0.650, 0.570, 0.500])
        self.prosp_C_0 = self.prosp_C_25
        self.prosp_C = np.array([self.prosp_C_100,self.prosp_C_75,self.prosp_C_50,self.prosp_C_25,self.prosp_C_0])
        self.fx = np.array([1,0.75,0.5,0.25,0])
        self.pcfx = np.array([15,20,25,30,35,40,45])

    def solveChillerAirtoWater(self,plant_info,heatFlow,t,T_ext):
        '''
        This method allows to calculate Air to Water Chiller electrical power
        following the Standard UNI-TS 11300:3 - 2010
        
               
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            heatFlow : float
                required power  [kW]
            t : int
                Simulation time step
            T_ext : float
                External temperature
                
        Returns
        -------
        None.        
        '''
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
        
        self.T_ext_solve = T_ext
        if self.T_ext_solve <= 15:
            self.T_ext_solve = 15.1
        elif self.T_ext_solve > 45:
            self.T_ext_solve = 45
        
        self.prosp_C_F = np.zeros(7)
            
        # Load Factor and EER in real conditions estimation 
        if heatFlow >= 0:
            pass
        else:
            self.F = heatFlow/(self.Q_max)            
            for i in range(len(self.x)-1):
                xa = self.x[i+1]
                xb = self.x[i]
                if xa < self.F <= xb:
                    self.EER = (self.F-xb)/(xa-xb)*self.y[i+1] - (self.F-xa)/(xa-xb)*self.y[i]
                elif self.F > 1:
                    self.EER = self.EER_100
                else:
                    pass
            for i in range(len(self.fx)-1):
                fxa = self.fx[i+1]
                fxb = self.fx[i]
                if fxa < self.F <= fxb:
                    for m in range(7):
                        self.prosp_C_F[m] = (self.F-fxb)/(fxa-fxb)*self.prosp_C[i+1][m] - (self.F-fxa)/(fxa-fxb)*self.prosp_C[i][m]
                elif self.F > 1:
                    self.prosp_C_F = self.prosp_C_100
                else:
                    pass
            for m in range(6):
                x1 = self.pcfx[m]
                x2 = self.pcfx[m+1]
                if x1 < self.T_ext_solve <= x2:
                    self.eta_1[t] = (self.T_ext_solve-x1)/(x2-x1)*self.prosp_C_F[m+1] - (self.T_ext_solve-x2)/(x2-x1)*self.prosp_C_F[m]
            
            # Electrical power required in real conditions 
            self.eta_mm[t] = self.EER*self.eta_1[t]
            self.W_el[t] = abs(heatFlow)/self.eta_mm[t]                        # [kWe]
            self.H_waste[t] = abs(heatFlow)+self.W_el[t]                       # [kW]
            self.W_aux[t] = self.W_aux_gn*(self.H_waste[t])                    # [kWe]
 
        
#%%--------------------------------------------------------------------------------------------------- 
#%% SplitAirConditioner class
class SplitAirConditioner:
    
    '''
    SplitAirConditioner class
    
    This method considers a generic Split Air Conditioner as the cooling plant 
    of the entire building following two documents of literature
    
    Berend Jan Christiaan van Putten, Nariman Mahdavi, Julio H. Braslavsky,
    An Analytical Model for Demand Response of Variable-Speed Air Conditioners,
    IFAC-PapersOnLine,
    Volume 51, Issue 28,
    2018,
    Pages 426-431,
    ISSN 2405-8963,
    https://doi.org/10.1016/j.ifacol.2018.11.740.
    (https://www.sciencedirect.com/science/article/pii/S2405896318334608)
    
    https://ieeexplore.ieee.org/document/7515217
    N. Mahdavi, J. H. Braslavsky and C. Perfumo, 
    "Mapping the Effect of Ambient Temperature on the Power Demand of Populations of Air Conditioners," 
    in IEEE Transactions on Smart Grid, vol. 9, no. 3, pp. 1540-1550, May 2018, 
    doi: 10.1109/TSG.2016.2592522.
    
            
    Methods:
        init
        setCondensingBoiler
        solvePlant
    '''
    
    def __init__(self,plant_type,l):
        '''
        init. Set some attributes for the method
        
        Berend Jan Christiaan van Putten, Nariman Mahdavi, Julio H. Braslavsky,
        An Analytical Model for Demand Response of Variable-Speed Air Conditioners,
        IFAC-PapersOnLine,
        Volume 51, Issue 28,
        2018,
        Pages 426-431,
        ISSN 2405-8963,
        https://doi.org/10.1016/j.ifacol.2018.11.740.
        (https://www.sciencedirect.com/science/article/pii/S2405896318334608)
        
        https://ieeexplore.ieee.org/document/7515217
        N. Mahdavi, J. H. Braslavsky and C. Perfumo, 
        "Mapping the Effect of Ambient Temperature on the Power Demand of Populations of Air Conditioners," 
        in IEEE Transactions on Smart Grid, vol. 9, no. 3, pp. 1540-1550, May 2018, 
        doi: 10.1109/TSG.2016.2592522.
        
        Parameters
            ----------
            plant_type : string
                string that sets the cooling plant
            l : int
                number of simulation time steps [-]
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary   
        
        # Input Data 
        self.l = l
        self.plant_type = plant_type #plant_type è una stringa
        self.T_ext_rif = 35                                                    # [°C]
        self.T_int_rif = 27                                                    # [°C]
        self.W_aux_gn = 0.04                                                   # [kW/kW_cond]
        self.W_el = np.zeros(l)                                                # [kWe]
        self.W_aux = np.zeros(l)                                               # [kWe]
        self.H_waste = np.zeros(l)                                             # [kW]
        self.eta_mm = np.zeros(l)
        self.EER = np.zeros(l)
        self.F = np.zeros(l)
        
        
    def setSplitAirConditioner(self,plant_info,Pnom):        #da modificare dentro il file excel di impianti
        
        '''
        Choise of plant size based on estimated nominal Power for air conditioner WITHOUT inverter
        
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            Pnom : float
                Design Cooling power  [kW]
                
        Returns
        -------
        None.
        '''        
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
   
        Size_0 = 0
        self.Pnom = - Pnom                                                     # [kW]
        for i in plant_info[self.plant_type].index:
            if Size_0 < self.Pnom <= plant_info[self.plant_type]['Size'][i]:
                self.ID = i
            Size_0 = plant_info[self.plant_type]['Size'][i]
        if self.Pnom > plant_info[self.plant_type]['Size'][plant_info[self.plant_type].index[-1]]:
            self.ID = len(plant_info[self.plant_type]) - 1
        self.Q_max = Pnom                                                      # [kW]
        
        if plant_info[self.plant_type]['inverter'][self.ID] == 'no':
            self.EER_targa = plant_info[self.plant_type]['EER_medio'][self.ID]  
        elif plant_info[self.plant_type]['inverter'][self.ID] == 'yes':                                                                  #Choise of plant size based on estimated nominal Power for air conditioner WITH inverter
            self.a = 11.2
            self.b = -0.363
            self.c = 2.6*pow(10,-3)
            self.d = 0.829
            self.e = 6.9*pow(10,-3)
            self.f = -0.159
        
    def solveSplitAirConditioner(self,plant_info,heatFlow,t,T_ext):
        '''
        This method allows to calculate from literature: Split Air Cooler condensing power
        
        Berend Jan Christiaan van Putten, Nariman Mahdavi, Julio H. Braslavsky,
        An Analytical Model for Demand Response of Variable-Speed Air Conditioners,
        IFAC-PapersOnLine,
        Volume 51, Issue 28,
        2018,
        Pages 426-431,
        ISSN 2405-8963,
        https://doi.org/10.1016/j.ifacol.2018.11.740.
        (https://www.sciencedirect.com/science/article/pii/S2405896318334608)

        https://ieeexplore.ieee.org/document/7515217
        N. Mahdavi, J. H. Braslavsky and C. Perfumo, 
        "Mapping the Effect of Ambient Temperature on the Power Demand of Populations of Air Conditioners," 
        in IEEE Transactions on Smart Grid, vol. 9, no. 3, pp. 1540-1550, May 2018, 
        doi: 10.1109/TSG.2016.2592522.
        
        Parameters
            ----------
            plant_info: dictionary
                dictionary with plant_key/Plants data
            heatFlow : float
                required power  [kW]
            t : int
                Simulation time step
            T_ext : float
                External temperature
                
        Returns
        -------
        None.        
        '''
        
        # Check input data type
        
        # The quality control is implemented in the Plant class.. Here is not necessary
        
        if heatFlow >= 0:
           pass 
        else:   
            self.F[t] = heatFlow/(self.Q_max)            
            if plant_info[self.plant_type]['inverter'][self.ID] == 'no':
                self.T_ext_solve = T_ext
                if self.T_ext_solve <= 20:
                    self.T_ext_solve = 20.1
                elif self.T_ext_solve > 35:
                    self.T_ext_solve = 35 
                self.EER[t] = -0.14*(self.T_ext_solve-35) + self.EER_targa
                
            elif plant_info[self.plant_type]['inverter'][self.ID] == 'yes':       #Choise of plant size based on estimated nominal Power for air conditioner WITH inverter'
                self.Pnom_cond = 7 #[kW] P_nom conditioner 
                self.number_cond = heatFlow/self.Pnom_cond #number of air conditioners suppose in each thermal zone
                self.EER[t]= self.a + self.b*T_ext + self.c*T_ext**2 + self.d*(heatFlow/self.number_cond) + self.e*T_ext*(heatFlow/self.number_cond) + self.f*(heatFlow/self.number_cond)**2
            
            self.H_waste[t] = (abs(heatFlow)*(1+self.EER[t]))/self.EER[t]
            self.W_el[t] = self.H_waste[t]-abs(heatFlow)            