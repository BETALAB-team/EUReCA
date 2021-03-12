'''
Urban Canyon class and other 
'''

'''IMPORTING MODULES'''
import numpy as np
# import os
import pandas as pd
from RC_classes.RSMDef import RSMDef
from RC_classes.param import Param
import os
from RC_classes.auxiliary_functions import wrn

# import pvlib
# import math
# import matplotlib.pyplot as plt

#%% ---------------------------------------------------------------------------------------------------
#%% Road class

class Road(object):
    '''
    Road class
    
    Methods  
        init
        solve_road
        update_temp_r
    '''
    
    # class variables
    
    road_albedo = 0.08                                                         # albedo strada
    alfa_road = 1 - road_albedo                                                # coefficiente di assorbimento della strada []  
    lamda_soil = 0.4                                                           # [W/m*K] conducibilità termica suolo

    def __init__(self,sim_time,suolo,strada,temp_layer_r,h_rd_sky):
        '''
        Initializes an Road object and sets some parameters
        
        Parameters
            ----------
            sim_time : list of int
                List of int with the number of time steps per hour (first), the number of hours (second) 
            suolo : np.array
                variable with soil layers properties (3 layers)
                example:
                np.array([
                    [3 , 1.4*10**6, 273 + 10],  
                    [3 , 1.4*10**6, 273 + 12],
                    [3 , 1.4*10**6, 273 + 13],
                    ], dtype = float), 
                    first column thickness [m]
                    second column thermal capacity [J/(m3 K)]
                    third column temperature [K]
            strada : np.array
                variable with road layers properties (3 layers)
                example:
                np.array([
                    [0.05, 1.9*10**6, 273 + 10],
                    [0.2 , 2*10**6,   273 + 12],
                    [1   , 1.4*10**6, 273 + 13],
                    ], dtype = float), 
                    first column thickness [m]
                    second column thermal capacity [J/(m3 K)]
                    third column temperature [K]
            temp_layer_r : np.array (4 components)
                array with road layers init temperatures [K]            
            h_rd_sky : float
                radiative heat transfer coefficient    [W/(m2 K)]
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        if not isinstance(sim_time, list) or not isinstance(sim_time[0], int) or not isinstance(sim_time[1], int):
            raise TypeError(f'ERROR road class climate_rev_v1, sim_time is not a list of int: sim_time {sim_time}') 
        if not isinstance(suolo, np.ndarray):
            raise TypeError(f'ERROR road class climate_rev_v1, input suolo is not a np.array: suolo {suolo}') 
        if not isinstance(strada, np.ndarray):
            raise TypeError(f'ERROR road class climate_rev_v1, input strada is not a np.array: strada {strada}') 
        if not isinstance(temp_layer_r, np.ndarray):
            raise TypeError(f'ERROR road class climate_rev_v1, input temp_layer_r is not a np.array: temp_layer_r {temp_layer_r}') 
        if not isinstance(h_rd_sky, float):
            raise TypeError(f'ERROR road class climate_rev_v1, input h_rd_sky is not a float: h_rd_sky {h_rd_sky}') 
        
        # Control input data quality
        
        if sim_time[0] > 4:
            wrn(f"WARNING Road class climate_rev_v1, input ts is higher than 4, this means more than 4 time steps per hours were set: ts {sim_time[0]}")
        if len(sim_time) > 2:
            wrn(f"WARNING Road class climate_rev_v1, sim_time is longer than 2, look the input comment: sim_time {sim_time}")
        if not suolo.shape == (3,3):
            wrn(f"WARNING Road class climate_rev_v1, suolo must be a 3,3 np.array, look the input comment: suolo {suolo}")
        if not strada.shape == (3,3):
            wrn(f"WARNING Road class climate_rev_v1, strada must be a 3,3 np.array, look the input comment: strada {strada}")
        if not temp_layer_r.shape == (4,):
            wrn(f"WARNING Road class climate_rev_v1, temp_layer_r must be a 4 np.array, look the input comment: temp_layer_r {temp_layer_r}")
        if h_rd_sky < 2 or h_rd_sky > 20:
            wrn(f"WARNING Road class climate_rev_v1, are you sure about the radiative heat tranfer coeff?? h_rd_sky {h_rd_sky}")
  
        ts = sim_time[0]
        hours = sim_time[1]
        self.l = ts*(hours-1) + 1
        self.h_rd_sky = h_rd_sky
        self.T_deep = temp_layer_r[0]
        self.T_r_0 = temp_layer_r[1]
        self.T_layer2_r_0 = temp_layer_r[2]
        self.T_layer3_r_0 = temp_layer_r[3]
        self.d_soil = suolo[:,0].sum()                                                   # [m] Profondità del suolo
        self.c_v_soil = suolo[2,1]                                               # [J/(m3*K)] Capacità termica volumetrica del suolo 
        self.C_soil = self.lamda_soil/self.d_soil 
        self.d_asphalt = strada[0,0]                                             # [m] profondità comunque da rivedere perchè in qquesto caso son 9.25m
        self.d_stones = strada[1,0]
        self.d_gravel = strada[2,0]
        self.c_v_gravel = strada[0,1]
        self.c_v_stones = strada[1,1]
        self.c_v_asphalt = strada[2,1]
        self.C_asphalt = 0.74/self.d_asphalt                                   # [W/(m2*K)] conduttanza asfalto 
        self.C_stones = 2.1/self.d_stones
        self.C_gravel = 0.4/self.d_gravel
        self.T_r_prova = np.zeros(self.l)
        self.T_layer3_r = np.zeros(self.l)
        self.T_layer2_r = np.zeros(self.l)
        
    def solve_road(self,t,T_urb_0,T_sky,h_r,S_r):
        '''
        solve the road system
        
        Parameters
            ----------
            t : int
                Simulation time step
            T_urb_0 : float
                Urban initial temperature [K]
            T_sky : float
                Sky vault temperature [K]
            h_r : float
                convective/radiative heat transfer coefficient with urban air volume [W/(m2 K)]
            S_r : float
                road surface [m2]
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        if not isinstance(t, np.int32):
            raise TypeError(f'ERROR Road class climate_rev_v1, t must be an integer: t {t}')
        if not isinstance(T_urb_0, float):
            raise TypeError(f'ERROR Road class climate_rev_v1, T_urb_0 must be an float: T_urb_0 {T_urb_0}')
        if not isinstance(T_sky, float):
            raise TypeError(f'ERROR Road class climate_rev_v1, T_sky must be an float: T_sky {T_sky}')
        if not isinstance(h_r, float):
            raise TypeError(f'ERROR Road class climate_rev_v1, h_r must be an float: h_r {h_r}')
        if not isinstance(S_r, float):
            raise TypeError(f'ERROR Road class climate_rev_v1, S_r must be an float: S_r {S_r}')
        
        # Control input data quality
        
        if T_urb_0 > 350. or T_urb_0 < 250. :
            wrn(f"WARNING Road class climate_rev_v1, T_urb_0 is outside the boudary limits [250,350] K: T_urb_0 {T_urb_0}")
        if T_sky > 350. or T_sky < 230. :
            wrn(f"WARNING Road class climate_rev_v1, T_sky is outside the boudary limits [230,350] K: T_sky {T_sky}")
        if h_r > 60.:
            wrn(f"WARNING Road class climate_rev_v1, h_r is outside the boudary limits (>50 W/(m2 K)): h_r {h_r}")
        
        # Road system solution, surface temperature calculation
        self.T_layer3_r[t] = (self.T_layer3_r_0*(self.d_gravel*self.c_v_gravel)/3600 + self.C_gravel*self.T_deep)/((self.d_gravel*self.c_v_gravel/3600) + self.C_gravel) #[K]
        self.T_layer2_r[t] = ((self.C_gravel*self.T_layer3_r[t]+self.C_asphalt*self.T_r_0)+(self.d_stones*self.c_v_stones)/3600*self.T_layer2_r_0)/((self.d_stones*self.c_v_stones)/3600+self.C_asphalt+self.C_soil)
        self.T_r = (self.T_r_0*(self.d_asphalt*self.c_v_asphalt)/3600 + self.C_asphalt*self.T_layer2_r[t] + S_r*self.alfa_road +self.h_rd_sky*T_sky + h_r*T_urb_0)/((self.d_asphalt*self.c_v_asphalt/3600) + self.C_asphalt + h_r + self.h_rd_sky) #[K]
        self.T_r_prova[t] = self.T_r
        
    def update_temp_r(self,t):
        '''
        update road temperature
        
        Parameters
            ----------
            t : int
                Simulation time step
                
        Returns
        -------
        None.
        
        '''  
        
        # Check input data type
        
        if not isinstance(t, np.int32):
            raise TypeError(f'ERROR Road class climate_rev_v1, t must be an integer: t {t}')             
        
        self.T_r_0 = self.T_r
        self.T_layer3_r_0 = self.T_layer3_r[t]
        self.T_layer2_r_0 = self.T_layer2_r[t]
        
        
#%% ---------------------------------------------------------------------------------------------------
#%% Soil class

class Soil():
    '''
    Soil class
    
    Methods  
        init
        solve_soil
        update_temp
    '''
    
    # class variables
    
    alfa = 0.7                                                                 # coefficiente di assorbimento del suolo []                                                       #[J/kg*K] cp/cv = 1.4
    lamda_soil = 0.4                                                           # [W/m*K] conducibilità termica suolo

    def __init__(self,sim_time,d_soil,c_v_soil,temp_layer,h_rd_sky):
        '''
        Initializes an Soil object and sets some parameters
        
        Parameters
            ----------
            sim_time : list of int
                List of int with the number of time steps per hour (first), the number of hours (second) 
            d_soil : np.array
                variable with soil layers thickness (3 layers)
                example:
                np.array([0.05, 0.2 , 1 ], dtype = float), 
                     thickness [m]
            c_v_soil : np.array
                variable with soil layers thickness (3 layers)
                example:
                np.array([1.9*10**6,2*10**6, 1.4*10**6,], dtype = float), 
                     thermal capacity [J/(m3 K)]
            temp_layer : np.array (4 components)
                array with road layers init temperatures [K]            
            h_rd_sky : float
                radiative heat transfer coefficient    [W/(m2 K)]
                
        Returns
        -------
        None.        
        '''  
        
        # Check input data type
        
        if not isinstance(sim_time, list) or not isinstance(sim_time[0], int) or not isinstance(sim_time[1], int):
            raise TypeError(f'ERROR road class climate_rev_v1, sim_time is not a list of int: sim_time {sim_time}') 
        if not isinstance(d_soil, np.ndarray):
            raise TypeError(f'ERROR road class climate_rev_v1, input suolo is not a np.array: d_soil {d_soil}') 
        if not isinstance(c_v_soil, np.ndarray):
            raise TypeError(f'ERROR road class climate_rev_v1, input strada is not a np.array: c_v_soil {c_v_soil}') 
        if not isinstance(temp_layer, np.ndarray):
            raise TypeError(f'ERROR road class climate_rev_v1, input temp_layer_r is not a np.array: temp_layer {temp_layer}') 
        if not isinstance(h_rd_sky, float):
            raise TypeError(f'ERROR road class climate_rev_v1, input h_rd_sky is not a float: h_rd_sky {h_rd_sky}') 
        
        # Control input data quality
        
        if sim_time[0] > 6:
            wrn(f"WARNING Soil class climate_rev_v1, input ts is higher than 6, this means more than 6 time steps per hours were set: sim_time {sim_time[0]}")
        if len(sim_time) > 2:
            wrn(f"WARNING Soil class climate_rev_v1, sim_time is longer than 2, look the input comment: sim_time {sim_time}")
        if not d_soil.shape == (3,):
            wrn(f"WARNING Soil class climate_rev_v1, d_soil must be a 3 np.array, look the input comment: d_soil {d_soil}")
        if not c_v_soil.shape == (3,):
            wrn(f"WARNING Soil class climate_rev_v1, c_v_soil must be a 3 np.array, look the input comment: c_v_soil {c_v_soil}")
        if not temp_layer.shape == (4,):
            wrn(f"WARNING Soil class climate_rev_v1, temp_layer must be a 4 np.array, look the input comment: temp_layer {temp_layer}")
        if h_rd_sky < 2 or h_rd_sky > 20:
            wrn(f"WARNING Soil class climate_rev_v1, are you sure about the radiative heat tranfer coeff?? h_rd_sky {h_rd_sky}")
  
        ts = sim_time[0]
        hours = sim_time[1]
        self.l = ts*(hours-1) + 1
        self.h_rd_sky = h_rd_sky
        self.T_deep = temp_layer[0]
        self.T_suolo_0 = temp_layer[1]
        self.T_layer2_0 = temp_layer[2]
        self.T_layer3_0 = temp_layer[3]
        self.d_soil = d_soil                                                   # [m] Profondità del suolo
        self.c_v_soil = c_v_soil                                               # [J/m^3*K] Capacità termica volumetrica del suolo 
        self.C_soil = self.lamda_soil/self.d_soil                              # [W/m^2*K] Trasmittanza suolo
        self.T_suolo_prova = np.zeros(self.l)
        
        
    def solve_soil(self,t,h_conv_rsl,T_air_rural,T_sky,Q_rad):
        '''
        Soil temperature calculation
        
        Parameters
            ----------
            t : int
                Simulation time step
            h_conv_rsl : float
                Convective heat transfer coefficient of the rural layers [W/(m2 K)]
            T_air_rural : float
                Rural air temperature [K]
            T_sky : float
                Sky vault temperature [K]
            Q_rad : float
                Irradiance on the horizontal plane rural [W/m2]
                
        Returns
        -------
        None.
        '''  
        
        # Check input data type
        
        if not isinstance(t, np.int32):
            raise TypeError(f'ERROR Soil class climate_rev_v1, t must be an integer: t {t}')
        if not isinstance(T_air_rural, float):
            raise TypeError(f'ERROR Soil class climate_rev_v1, T_air_rural must be an float: T_air_rural {T_air_rural}')
        if not isinstance(T_sky, float):
            raise TypeError(f'ERROR Soil class climate_rev_v1, T_sky must be an float: T_sky {T_sky}')
        if not isinstance(h_conv_rsl, float):
            raise TypeError(f'ERROR Soil class climate_rev_v1, h_conv_rsl must be a float: h_conv_rsl {h_conv_rsl}')
        if not isinstance(Q_rad, float):
            raise TypeError(f'ERROR Soil class climate_rev_v1, Q_rad must be a float: Q_rad {Q_rad}')
        
        # Control input data quality
        
        if T_air_rural > 350. or T_air_rural < 250. :
            wrn(f"WARNING Soil class climate_rev_v1, T_air_rural is outside the boudary limits [250,350] K: T_air_rural {T_air_rural}")
        if T_sky > 350. or T_sky < 230. :
            wrn(f"WARNING Soil class climate_rev_v1, T_sky is outside the boudary limits [230,350] K: T_sky {T_sky}")
        if h_conv_rsl > 60.:
            wrn(f"WARNING Soil class climate_rev_v1, h_conv_rsl is outside the boudary limits (>50 W/(m2 K)): h_conv_rsl {h_conv_rsl}")
        if Q_rad < 0. or Q_rad > 2000.:
            wrn(f"WARNING Soil class climate_rev_v1, Q_rad is outside the boudary limits (>50 W/(m2 K)): Q_rad {Q_rad}")        
        
        # System solution
        
        self.T_layer3 = (self.T_layer3_0*(self.d_soil[2]*self.c_v_soil[2])/3600 + self.C_soil[2]*self.T_deep)/((self.d_soil[2]*self.c_v_soil[2]/3600) + self.C_soil[2]) #[K]
        self.T_layer2 = (self.C_soil[1]*(self.T_layer2_0+self.T_suolo_0)+self.T_layer2_0*(self.d_soil[1]*self.c_v_soil[1])/3600)/((self.d_soil[1]*self.c_v_soil[1])/3600+2*self.C_soil[1])
        self.T_suolo = (self.T_suolo_0*(self.d_soil[0]*self.c_v_soil[0])/3600 + self.C_soil[0]*self.T_layer2+ Q_rad*self.alfa + h_conv_rsl*T_air_rural + self.h_rd_sky*T_sky)/(self.d_soil[0]*self.c_v_soil[0]/3600 + self.C_soil[0] + h_conv_rsl + self.h_rd_sky) #[K]
        self.T_suolo_prova[t] = self.T_suolo
    
    def update_temp(self):
        '''
        update soil temperature
        
        Parameters
            ----------
            None.
            
        Returns
        -------
        None.
        
        '''  

        self.T_suolo_0 = self.T_suolo
        self.T_layer2_0 = self.T_layer2
        self.T_layer3_0 = self.T_layer3


#%% ---------------------------------------------------------------------------------------------------
#%% Urban canyon class

class UrbanCanyon(object):
    '''
    Urban canyon class
    Look at:
    
    B. Bueno, L. Norford, J. Hidalgo, and G. Pigeon, “The urban weather generator,” 
    J. Build. Perform. Simul., vol. 6, no. 4, pp. 269–281, 2013, 
    doi: 10.1080/19401493.2012.718797.
    
    B. Bueno, L. Norford, G. Pigeon, and R. Britter, 
    “A resistance-capacitance network model for the analysis of the interactions 
    between the energy performance of buildings and the urban climate,” 
    Build. Environ., vol. 54, pp. 116–125, Aug. 2012, 
    doi: 10.1016/J.BUILDENV.2012.01.023.
    
    Methods  
        init
        solve_canyon
        variante_u_exc
    '''

    # Class variables

    # Input data
    ro_air = 1.225                                                             # [kg/m^3]
    cp_air = 1005.                                                            # [J/kg*K]
    cv_air = cp_air*ro_air                                                     # [J/m3*K]
    h_rd_sky = 5.                                                              # [w/m^2*K] =4*emissivita_vapore*1* 5.67*10^-8*(T_medio_can,j)^3
    
    H_urb_0 = 10.                                                              # [W/m^2]
    T_urb_0 = 24+273.15                                                        # [K]
    T_boundary_0 = 18+273.15                                                   # [K
    
    def __init__(self,sim_time,data,meso_path = os.path.join('.','RC_classes')):
        '''
        Initializes an Urban Canyon and sets some parameters
        
        Parameters
            ----------
            sim_time : list of int
                List of int with the number of time steps per hour (first), the number of hours (second) 
            data : dictionary
                this dictionary contains all urban canyon input data. Example:
                UWG_data = {
                'Area' : float(100000),
                'Diameter': float(200),
                'Perimeter': float(800),
                'Ortogonal_length': float(200),
                'h_layer_inf': int(2),
                'h_layer_sup': int(150),
                'z_i_day': int(1000),
                'z_i_night': int(50),
                'z_i_m': int(10),
                'z_i_ref': int(150),
                'ExtWall_radiative_coef': float(5),
                'ExtWall_convective_coef': float(17),
                'ExtWall_emissivity': float(0.85),
                'ExtWall_reflection': float(0.15),
                'Average_U_wondows': float(2),
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
                'T_deep': float(13+273)    
    
                'cars_number_over_1000_px': float(39545000/60360),
                'motorcycles_number_over_1000_px': float(6896000/60360),
                'othervehicle_number_over_1000_px': float(5364000/60360),
                'emission_factor_cars_over_1000_px': float(	24.74 ),
                'emission_factor_motorcycles_over_1000_px': float(13.34),
                'emission_factor_othervehicle_over_1000_px': float(104.95),
                'vehicle_multiplying_factor': float(1),
                'distance_per_hour': float(64000),
                'population_density': float(2283.2),
                'first_day': Sim_input['first_day'],
                'Hw_week': np.array([0.9,0.5,0.9,2.1,3.3,4.2,5.8,7.1,6.3,5.8,5.4,5.4,5.8,6.3,6.7,7.1,7.5,5.8,4.2,3.3,2.5,1.7,1.3,0.9])/100,
                'Hw_end': np.array([1.7,1.2,1.,1.7,2.6,3.3,4.2,5,5.8,6.2,6.7,6.7,6.7,6.7,6.3,5.8,5,4.6,4.2,3.7,3.3,2.9,2.5,2.1])/100                    
                }        
                
            meso_path : string
                path containing the meso file
                
        Returns
        -------
        None.        
        '''  
        
        # Check input data type
        
        if not isinstance(sim_time, list) or not isinstance(sim_time[0], int) or not isinstance(sim_time[1], int):
            raise TypeError(f'ERROR Urban Canyon class climate_rev_v1, sim_time is not a list of int: sim_time {sim_time}') 
        if not isinstance(data, dict):
            raise TypeError(f'ERROR Urban Canyon class data input must be a dictionary : data {data}') 
        if not isinstance(meso_path, str):
            raise TypeError(f'ERROR Urban Canyon class meso_path input must be a string : meso_path {meso_path}') 
        
        ts = sim_time[0]
        hours = sim_time[1]
        self.data = data
        self.l = ts*(hours-1) + 1
        
        """
        self.nome = 'Padova'
        self.vettore = np.zeros(25)
        self.vettore2 = np.arange(10)
        self.area = 100
        self.height = 10
        self.volume = self.area * self.height
        """
        
        try:
            temp_layer = np.array([data['T_deep'],
                                   data['Soil_layers'][2,2],
                                   data['Soil_layers'][1,2],
                                   data['Soil_layers'][0,2]])
            self.soil = Soil(sim_time,data['Soil_layers'][:,0],data['Soil_layers'][:,1],temp_layer,self.h_rd_sky)
            
            temp_layer_r = np.array([data['T_deep'],
                           data['Road_layers'][2,2],
                           data['Road_layers'][1,2],
                           data['Road_layers'][0,2]])
            self.road = Road(sim_time,data['Soil_layers'],data['Road_layers'],temp_layer_r,self.h_rd_sky)

            self.A_w = self.data['tot_wall_glazed_area']\
                        + self.data['tot_wall_opaque_area']                         # [m^2] area totale pareti esterne
            self.emissivity_w = self.data['ExtWall_emissivity']                                        # emissività pareti estrerne
            self.riflex_w = self.data['ExtWall_reflection']                                               # riflessività pareti estrerne
            self.h_rad_w = self.data['ExtWall_radiative_coef']                                                 # [W/m^2*K] coefficiente scambio radiativo pareti
            self.h_rad_r = self.h_rd_sky
            self.A_win = self.data['tot_wall_glazed_area']                                                # [m^2] Area totale finestre
            self.U_win = self.data['Average_U_windows']                                                     # [W/m^2*K] Trasmittanza finestre
            
            self.VH_urb = self.data['VH_urb']                                                   # rapporto tra: area verticale degli edifici/area del sito urbano -> valori presi dai dati di Basilea di Bueno
            self.horiz_buil_densi = self.data['Buildings_density']                                # percentuale densità edifici 
            self.horiz_veg_densi = self.data['Vegetation_density']                                  # percentuale densità di vegetazione(alberi)
            self.h_aver_buil = self.data['h_builidng_average']                                          # [m] altezza media edifici
            self.D_city = self.data['Diameter']                                                     # [m] Diametro della città
            self.P_city = self.data['Perimeter']                                        # [m] Perimetro della città
            self.A_city = self.data['Area']                                   # [m^2] Area della città
            self.W = self.data['Ortogonal_length']                                 # [m] lunghezza della città ortogonale alla direzione del vento 
            self.A_canyon = self.A_city-self.horiz_buil_densi*self.A_city          # [m^2] Area del canyon
            self.V_canyon = self.A_canyon * self.h_aver_buil                       # [m^3] Volume del canyon
            self.w_aver_r = 2*self.h_aver_buil*(1-self.horiz_buil_densi)/self.VH_urb # [m] lunghezza media strade
            self.A_r = self.A_city *(1-self.horiz_veg_densi-self.horiz_buil_densi) # [m^2] Area della strada

            self.altezza_layer_inf = self.data['h_layer_inf']                              # [m] 2 -> la temperatura rurale è misurata alla stazione meteorologica a 2m
            self.altezza_layer_sup = self.data['h_layer_sup']                              # [m] 150 -> il gradiente velocità arriva fino a 150m
            self.altezza = np.arange(self.altezza_layer_inf,self.altezza_layer_sup+1,1)      # [m] altezza del vettore: gradiente velocità
            self.zi_day = self.data['z_i_day']                                                    # [m] 1000 -> valori presi dal documento di Bueno, daytime boundary layer height
            self.zi_night = self.data['z_i_night']                                                # [m] 50 -> nighttime boundary layer height
            self.z_m = self.data['z_i_m']                                                          # [m] 10 -> quando siamo in condizione di vento geostrofico, la velocità di riferimento è misurata alla stazione meteorologica a 10m
            self.z_ref = self.data['z_i_ref']                                                      # [m] 150-> altezza di riferimento nel quale il profilo di temperatura è uniforme
        except KeyError:
            raise KeyError(f"""ERROR There's something wrong with the creation of the urban canyon. I cannot find all the necessary input data.
            See the __init__ help to look at the required variables           
            """)
        
        
        # Vectors initialization 
        self.C_adv = np.zeros(self.l)
        self.arg_ucirc = np.zeros(self.l)
        self.u_circ = np.zeros(self.l)
        self.H_rur = np.zeros(self.l)
        self.diff = np.zeros(self.l)
        self.H_urb = np.zeros(self.l)
        self.T_urb = np.zeros(self.l)
        self.T_boundary = np.zeros(self.l)
        self.err_H_urb = np.zeros(self.l)
        self.Ri = np.zeros(self.l)
        self.h_r = np.zeros(self.l)
        self.S_r = np.zeros(self.l)
        self.f_m = np.zeros(self.l)
        self.H_waste = np.zeros(self.l)
        self.T_w = np.zeros(self.l)
        
        # Traffic calculation
        try:
            emission = (data['cars_number_over_1000_px']*data['emission_factor_cars_over_1000_px']
                    + data['motorcycles_number_over_1000_px']*data['emission_factor_motorcycles_over_1000_px']
                    + data['othervehicle_number_over_1000_px']*data['emission_factor_othervehicle_over_1000_px'])/1000 #over 1 px conversion
            
            H_traffic_nom = emission*data['vehicle_multiplying_factor']*24*data['distance_per_hour']*data['population_density']/3600000000
            week = np.hstack((np.tile(data['Hw_week'],5), np.tile(data['Hw_end'],2)))
            first_week = week[(data['first_day']-1)*24 :]
        except KeyError:
            raise KeyError(f"""ERROR There's something wrong with the creation of the urban canyon. I cannot find all the necessary input data.
            See the __init__ help to look at the required variables           
            """)        
        
        year = np.hstack((first_week, np.tile(week,53)))[:8760]
        m = str(60/ts) + 'min'
        timeIndex = pd.date_range('2020-01-01', periods=8760, freq='1h')
        year_df = pd.DataFrame(year).set_index(timeIndex)         
        year_df = year_df.resample(m).interpolate(method='linear')              # Steps interpolation Resampling
        self.H_traffic = year_df.to_numpy()[:,0] * H_traffic_nom * data['Area']   # [W]
        
        self.RSMParams = Param()
        self.RSM = RSMDef(45.41,9.81,1,2.,4.,101026,self.RSMParams,meso_path)
        
        self.vettore_prova = np.zeros(self.l)
        self.vettore_prova_2 = np.zeros(self.l)
        self.vettore_prova_3 = np.zeros(self.l)
        self.vettore_prova_4 = np.zeros(self.l)
        
        
    def solve_canyon(self,t,T_air_rural,u_atm_input,radiation,T_w,T_sky,H_waste,lamda_zenith,T_exf,V_exf):             
        '''
        Solve the canyon balance
        
        B. Bueno, L. Norford, J. Hidalgo, and G. Pigeon, “The urban weather generator,” 
        J. Build. Perform. Simul., vol. 6, no. 4, pp. 269–281, 2013, 
        doi: 10.1080/19401493.2012.718797.
        
        B. Bueno, L. Norford, G. Pigeon, and R. Britter, 
        “A resistance-capacitance network model for the analysis of the interactions 
        between the energy performance of buildings and the urban climate,” 
        Build. Environ., vol. 54, pp. 116–125, Aug. 2012, 
        doi: 10.1016/J.BUILDENV.2012.01.023.        
        
        Parameters
            ----------
            t : int
                Simulation time step
            T_air_rural : float
                Rural air temperature [°C]
            u_atm_input : float
                Wind velocity [m/s]
            radiation : list of floats
                Global and direct irradiance     
            T_w : float
                Wall temperature [°C]
            T_sky : float
                Sky vault temperature [K]
            H_waste : float
                Heat flux from buildings [W]
            lamda_zenith : SolarPosition.zenith[t] float
                Solar zenith from SolarPosition.zenith[t] (WeatherData.SolarPosition object)
            T_exf : list of floats
                List of floats with infiltration and ventilation temperature
            V_exf : list of floats
                List of floats with infiltration and ventilation volume flow rate
                
        Returns
        -------
        None.
        '''  
        
        # Check input data type
        
        if not isinstance(t, np.int32):
            raise TypeError(f'ERROR Urban canyon class climate_rev_v1, t must be an integer: t {t}')
        if not isinstance(T_air_rural, float):
            raise TypeError(f'ERROR Urban canyon class climate_rev_v1, T_air_rural must be an float: T_air_rural {T_air_rural}')
        if not isinstance(u_atm_input, float):
            raise TypeError(f'ERROR Urban canyon class climate_rev_v1, u_atm_input must be an float: u_atm_input {u_atm_input}')
        if not isinstance(T_w, float):
            raise TypeError(f'ERROR Urban canyon class climate_rev_v1, T_w must be an float: T_w {T_w}')
        if not isinstance(T_sky, float):
            raise TypeError(f'ERROR Urban canyon class climate_rev_v1, T_sky must be an float: T_sky {T_sky}')
        if not isinstance(T_exf, list):
            raise TypeError(f'ERROR Urban canyon class climate_rev_v1, T_exf must be a list: T_exf {T_exf}')
        if not isinstance(V_exf, list):
            raise TypeError(f'ERROR Urban canyon class climate_rev_v1, V_exf must be a list: V_exf {V_exf}')
        if not isinstance(radiation, list) or not isinstance(radiation[0], float) or not isinstance(radiation[1], float):
            raise TypeError(f'ERROR Urban canyon class climate_rev_v1, radiation must be a list of floats: radiation {radiation}')
        if not isinstance(lamda_zenith, float):
            raise TypeError(f'ERROR Urban canyon class climate_rev_v1, lamda_zenith must be a float: lamda_zenith {lamda_zenith}')

        # Check input data quality 
        
        if T_air_rural < -20. or T_air_rural > 50.:
            wrn(f"WARNING Urban canyon class solve method, time step {t}: T_air_rural outside limit [-20.,50.]: T_air_rural, {T_air_rural}")
        if u_atm_input < 0. or u_atm_input > 15.:
            wrn(f"WARNING Urban canyon class solve method, time step {t}: u_atm_input outside limit [0.,15.]: u_atm_input, {u_atm_input}")
        if T_w < -20. or T_w > 50.:
            wrn(f"WARNING Urban canyon class solve method, time step {t}: T_w outside limit [-20.,50.]: T_w, {T_w}")
        if T_sky < -30. or T_sky > 50.:
            wrn(f"WARNING Urban canyon class solve method, time step {t}: T_sky outside limit [-20.,50.]: T_sky, {T_sky}")
        if T_exf[0] < -20. or T_exf[0] > 50.:
            wrn(f"WARNING Urban canyon class solve method, time step {t}: T_exf outside limit [-20.,50.]: T_exf, {T_exf}")
        if T_exf[1] < -20. or T_exf[1] > 50.:
            wrn(f"WARNING Urban canyon class solve method, time step {t}: T_exf outside limit [-20.,50.]: T_exf, {T_exf}")
        if V_exf[0] < 0.:
            wrn(f"WARNING Urban canyon class solve method, time step {t}: V_exf negative value: V_exf, {V_exf}")
        if radiation[0] < 0. or radiation[0] > 2000.:
            wrn(f"WARNING Urban canyon class solve method, time step {t}: radiation outside limit [0.,2000.]: radiation, {radiation}")
        if radiation[1] < 0. or radiation[1] > 2000.:
            wrn(f"WARNING Urban canyon class solve method, time step {t}: radiation outside limit [0.,2000.]: radiation, {radiation}")

        # Set some input parameter

        T_w += 273.15
        T_air_rural += 273.15
        T_sky += 273.15
        T_exf_inf, T_exf_vent = T_exf
        T_exf_inf += 273.15
        T_exf_vent += 273.15
        V_exf_inf, V_exf_vent = V_exf
        
        Q_rad, Q_dir = radiation
        Q_diff = Q_rad - Q_dir
        
        flag_H_urb = True
        it = 0
        
        while flag_H_urb:                                                      # Fintanto che flag_H_urb sia True viene ripetuto il calcolo nello stesso timestep!!
        
            # Wind velocity control
            if u_atm_input > 0.1:
               self.u_atm = u_atm_input
            else: 
               self.u_atm = 0.1  
            # self.u_atm = .5
            self.h_conv_rsl = 5.8+3.7*self.u_atm                               # [W/(m2*K)]
           
            # 1.Equation RSM 
            self.soil.solve_soil(t,self.h_conv_rsl,T_air_rural,T_sky,Q_rad)
            self.H_rur[t] = self.h_conv_rsl*(T_air_rural-self.soil.T_suolo_prova[t]) # [W/m2] perchè uscente dal nodo T_rur
        
            # 2.Equation VDM 
            self.T_altezza = (T_air_rural-0.0065*self.altezza)                 # [K] una matrice di 24 colonne e 150 righe(perchè si hanno 150m di altitudine)
            
            
            # This is taken from the original UWG model....
            # The vertical diffusion model is quite important
            self.RSM.VDM(T_air_rural, 101325, self.u_atm, self.H_rur[t], self.RSMParams)
            
            
            self.T_altezza_mod =  np.interp(np.arange(0, len(self.altezza)), 
                      np.arange(0, len(np.array(self.RSM.tempProf)))*len(self.altezza)/len(self.RSM.tempProf), self.RSM.tempProf)
            
            self.T_altezza = self.T_altezza_mod
            
            # 3. Equation UBL 
            self.kw = 1
            self.beta = 9.81/self.T_boundary_0
            self.diff[t] = self.H_urb_0-self.H_rur[t]
            if self.diff[t] > 0:
                if Q_rad < 10:
                    self.arg_ucirc[t] = (self.beta*self.zi_night*(self.H_urb_0-self.H_rur[t]))/(self.ro_air*self.cp_air)
                    self.u_circ[t] = self.kw*self.arg_ucirc[t]**(1/3)
                    if self.u_circ[t] < self.u_atm:
                        self.C_surf = (self.H_urb_0*3600)/(self.zi_night*self.ro_air*self.cv_air)                            # [K]
                        self.C_adv[t]  = (self.u_atm*self.zi_night*3600*self.cp_air)/(2*self.z_m*self.D_city*self.cv_air)  # []
                        self.teta_eq = 2/3*self.T_altezza[self.zi_night]+1/3*self.T_altezza[self.altezza_layer_inf]          # [K]
                        self.T_boundary[t] = (self.T_boundary_0+self.C_surf+self.C_adv[t]*self.teta_eq)/(self.C_adv[t]+1)    # [K]
                    else:
                        self.C_surf = (self.H_urb_0*3600)/(self.zi_night*self.ro_air*self.cv_air)                            # [K]
                        self.C_adv[t] = (self.P_city*self.u_circ[t]*3600*self.cp_air)/(self.A_city*self.cv_air)              # []
                        self.teta_eq = 1/2*self.T_altezza[self.zi_night]+1/2*self.T_altezza[self.altezza_layer_inf]          # [K]
                        self.T_boundary[t] = (self.T_boundary_0+self.C_surf +self.C_adv[t]*self.teta_eq)/(self.C_adv[t]+1)   # [K]
                else:
                    self.arg_ucirc[t] = (self.beta*self.zi_day*(self.H_urb_0-self.H_rur[t]))/(self.ro_air*self.cp_air)
                    self.u_circ[t] = self.kw*pow(self.arg_ucirc[t],1/3)
                    if self.u_circ[t] < self.u_atm:
                        self.C_surf = (self.H_urb_0*3600)/(self.zi_day*self.ro_air*self.cv_air)                              # [K]
                        self.C_adv[t] = (self.W*self.u_atm*3600*self.cp_air)/(self.A_city*self.cv_air)                       # []
                        self.teta_eq = self.T_altezza[self.altezza_layer_sup-self.altezza_layer_inf]                         # [K] 
                        self.T_boundary[t] = (self.T_boundary_0+self.C_surf+self.C_adv[t]*self.teta_eq)/(self.C_adv[t]+1)    # [K]
                    else:
                        self.C_surf = (self.H_urb_0*3600)/(self.zi_day*self.ro_air*self.cv_air)                              # [K]
                        self.C_adv[t] = (self.P_city*self.u_circ[t]*3600*self.cp_air)/(self.A_city*self.cv_air)              # []           
                        self.teta_eq = self.T_altezza[self.altezza_layer_sup-self.altezza_layer_inf]                         # [K]
                        self.T_boundary[t] = (self.T_boundary_0+self.C_surf+self.C_adv[t]*self.teta_eq)/(self.C_adv[t]+1)    # [K]
            else:
                self.u_circ[t] = 0
                if Q_rad < 1:
                    if self.u_circ[t] < self.u_atm:
                        self.C_surf = (self.H_urb_0*3600)/(self.zi_night*self.ro_air*self.cv_air)                            # [K]
                        self.C_adv[t]  = (self.u_atm*self.zi_night*3600*self.cp_air)/(2*self.z_m*self.D_city/4*self.cv_air)  # []
                        self.teta_eq = 2/3*self.T_altezza[self.zi_night]+1/3*self.T_altezza[self.altezza_layer_inf]          # [K]
                        self.T_boundary[t] = (self.T_boundary_0+self.C_surf+self.C_adv[t]*self.teta_eq)/(self.C_adv[t]+1)    # [K]
                    else:
                        self.C_surf = (self.H_urb_0*3600)/(self.zi_night*self.ro_air*self.cv_air)                            # [K]
                        self.C_adv[t] = (self.P_city*self.u_circ[t]*3600*self.cp_air)/(self.A_city*self.cv_air)              # [] 
                        self.teta_eq = 1/2*self.T_altezza[self.zi_night]+1/2*self.T_altezza[self.altezza_layer_inf]          # [K]
                        self.T_boundary[t] = (self.T_boundary_0+self.C_surf +self.C_adv[t]*self.teta_eq)/(self.C_adv[t]+1)   # [K]
                else:
                    
                    if self.u_circ[t] < self.u_atm:
                        self.C_surf = (self.H_urb_0*3600)/(self.zi_day*self.ro_air*self.cv_air)                              # [K]
                        self.C_adv[t] = (self.W*self.u_atm*3600*self.cp_air)/(self.A_city*self.cv_air)                       # []
                        self.teta_eq = self.T_altezza[self.altezza_layer_sup-self.altezza_layer_inf]                         # [K]
                        self.T_boundary[t] = (self.T_boundary_0+self.C_surf+self.C_adv[t]*self.teta_eq)/(self.C_adv[t]+1)    # [K]
                    else:
                        self.C_surf = (self.H_urb_0*3600)/(self.zi_day*self.ro_air*self.cv_air)                              # [K]
                        self.C_adv[t] = (self.P_city*self.u_circ[t]*3600*self.cp_air)/(self.A_city*self.cv_air)              # []
                        self.teta_eq = self.T_altezza[self.altezza_layer_sup-self.altezza_layer_inf]                         # [K]
                        self.T_boundary[t] = (self.T_boundary_0+self.C_surf+self.C_adv[t]*self.teta_eq)/(self.C_adv[t]+1)    # [K]
    
            # 4. Equation UC 
            '''
            self.C_vk = 0.4                                                    # Von-Karman constant
            self.zr = 1.5*self.h_aver_buil                                     # blending height
            self.z0 = 0.1*self.h_aver_buil                                     # roughness length
            
            self.Ri[t] = (9.81*self.zr*(self.T_boundary[t]-self.T_urb_0))/(self.T_boundary[t]*pow(self.u_atm,2)) # Numero di Richardson (vedere in letterattura se u_atm=0)
            self.a = self.C_vk/(np.log(self.zr/self.z0))                       # Drag coefficient
            if self.Ri[t] >= 0:
                self.f_m[t] = 1/(1+4.7*self.Ri[t])**2                          # Coefficient for atmosphere stability if Ri>=0 (non si mescolano molto i 2 fluidi)
            else:
                self.c = 69.56*pow(self.a,2)*pow(self.zr/self.z0,0.5)
                self.f_m[t] = (1-9.4*self.Ri[t])/(1+self.c*(-self.Ri[t])**2)   # Coefficient for atmosphere unstability Ri<0
            
            self.u_friction = self.a*self.u_atm*pow(self.f_m[t],0.5)
            self.u_can = self.u_friction*pow(8/self.VH_urb,0.5)
            #self.u_ex = self.u_friction/(self.u_atm/self.u_friction-pow(8/self.VH_urb,0.5))  # velocità di scambio
            '''
            self.u_ex, self.u_friction = self.variante_u_exc(self.u_atm,self.H_urb_0,self.T_urb_0) 
            self.u_can = self.u_friction*pow(8/self.VH_urb,0.5)
            
            self.h_w = 5.8+3.7*self.u_atm+self.h_rad_w                         # [W/(m2*K)]
            self.h_r[t] = 5.8+3.7*self.u_can+self.h_rad_r                      # [W/(m2*K)]
            # zenith = np.min([lamda_zenith,float(89)])*np.pi/180
            zenith = np.min([lamda_zenith,89.])*np.pi/180
            self.argomento = self.w_aver_r/(self.h_aver_buil*np.tan(zenith))
            # print(self.w_aver_r,self.h_aver_buil,self.argomento, zenith,lamda_zenith)
            if self.argomento<1:
                self.teta0 = np.arcsin(self.argomento)
            else:
                self.teta0 = np.arcsin(1)
                    
            # Radiation on the road
            self.F_r = pow(pow(self.h_aver_buil/self.w_aver_r,2)+1,0.5)-(self.h_aver_buil/self.w_aver_r)
            self.F_w = 0.5*(self.h_aver_buil/self.w_aver_r+1-pow(pow(self.h_aver_buil/self.w_aver_r,2)+1,0.5))/(self.h_aver_buil/self.w_aver_r)
            self.K_r = Q_dir*(2*self.teta0/np.pi+2/np.pi*self.h_aver_buil/self.w_aver_r*np.tan(zenith)*(1-np.cos(self.teta0)))+self.F_r*Q_diff
            self.K_w = Q_dir*(self.w_aver_r/self.h_aver_buil*(0.5-self.teta0/np.pi)+ 1/np.pi*np.tan(zenith)*(1-np.cos(self.teta0)))+self.F_w*Q_diff
            self.R_w = self.riflex_w*self.K_w
            self.R_r = self.road.road_albedo*self.K_r
            self.M_w = (self.R_w+self.F_w*self.riflex_w*self.R_r)/(1-(1-2*self.F_w)*self.riflex_w+(1-self.F_r)*self.F_w*self.road.road_albedo*self.riflex_w)
            self.M_r = (self.R_r+(1-self.F_r)*self.road.road_albedo*(self.R_w+self.F_w*self.riflex_w*self.R_r))/(1-(1-2*self.F_w)*self.riflex_w+(1-self.F_r)*self.F_w*self.road.road_albedo*self.riflex_w)
            self.S_r[t] = self.K_r+(1-self.F_r)*self.M_w
            self.S_w = self.K_w+(1-2*self.F_w)*self.M_w + self.F_w*self.M_r
    
            self.road.solve_road(t,self.T_urb_0,T_sky,self.h_r[t],self.S_r[t])
        
            # Urban temeprature calculation (Air node thermal balance)
            self.T_urb[t] = (self.A_w*self.h_w*T_w +
                             self.A_r*self.h_r[t]*self.road.T_r_prova[t] + 
                             self.A_r*self.h_rd_sky*T_sky + 
                             self.A_win*self.U_win*T_exf_inf + 
                             V_exf_inf*self.ro_air*self.cp_air*T_exf_inf + 
                             V_exf_vent*self.ro_air*self.cp_air*T_exf_vent + 
                             self.A_canyon*self.u_ex*self.ro_air*self.cp_air*self.T_boundary[t]+ 
                             H_waste + 
                             self.H_traffic[t] + 
                             self.T_urb_0*self.V_canyon*self.ro_air*self.cv_air/3600)\
                             /(self.V_canyon*self.ro_air*self.cv_air/3600 + 
                              self.A_w*self.h_w + 
                              self.A_r*self.h_r[t] + 
                              self.A_r*self.h_rd_sky + 
                              self.A_win*self.U_win + 
                              (V_exf_inf + V_exf_vent)*self.ro_air*self.cp_air + 
                              self.A_canyon*self.u_ex*self.ro_air*self.cp_air) #[K]
                            
                            
            self.H_urb[t] = self.u_ex*self.ro_air*self.cp_air*(self.T_urb[t]-self.T_boundary[t])#+ self.h_w*(T_w-self.T_boundary[t])+ (H_waste + self.H_traffic[t])/self.A_canyon #[W/m^2] flusso convettivo tra T_urb-T_ubl => uscente dal canyon verso urban boundary layer
            
            self.H_waste[t] = H_waste
            self.T_w[t] = T_w
            
            # H_urb control 
            self.err_H_urb[t] = abs(self.H_urb_0 - self.H_urb[t])
            
            if self.err_H_urb[t] > 10 and it < 10:
                self.H_urb_0 = self.H_urb[t]
                it += 1
            else:
                self.soil.update_temp()
                self.road.update_temp_r(t)
                self.T_urb_0 = self.T_urb[t]
                self.H_urb_0 = self.H_urb[t]
                self.T_boundary_0 = self.T_boundary[t]
                flag_H_urb = False                                             # Se l'errore è minore della tolleranza allora flag_H_urb diventa False e si passa al timestep successivo!
            '''    
            self.vettore_prova[t] = (
                             -1*self.A_canyon*self.u_ex*self.ro_air*self.cp_air*(self.T_boundary[t] -self.T_urb[t]))
            
            self.vettore_prova_2[t] = (self.A_w*self.h_w*(T_w -self.T_urb[t])+
                             self.A_r*self.h_r[t]*(self.road.T_r_prova[t] -self.T_urb[t])+ 
                             self.A_r*self.h_rd_sky*(T_sky -self.T_urb[t]) + 
                             self.A_win*self.U_win*(T_exf_inf -self.T_urb[t]) + 
                             V_exf_inf*self.ro_air*self.cp_air*(T_exf_inf -self.T_urb[t])+ 
                             V_exf_vent*self.ro_air*self.cp_air*(T_exf_vent -self.T_urb[t])+ 
                             H_waste + 
                             self.H_traffic[t])
            
            self.vettore_prova_3[t] = (self.A_w*self.h_w*(T_w -self.T_urb[t])+
                             self.A_r*self.h_r[t]*(self.road.T_r_prova[t] -self.T_urb[t])+ 
                             self.A_r*self.h_rd_sky*(T_sky -self.T_urb[t]) + 
                             self.A_win*self.U_win*(T_exf_inf -self.T_urb[t]) + 
                             V_exf_inf*self.ro_air*self.cp_air*(T_exf_inf -self.T_urb[t])+ 
                             V_exf_vent*self.ro_air*self.cp_air*(T_exf_vent -self.T_urb[t])+ 
                             self.A_canyon*self.u_ex*self.ro_air*self.cp_air*(self.T_boundary[t] -self.T_urb[t])+ 
                             H_waste + 
                             self.H_traffic[t])
            '''
            self.vettore_prova[t] = (self.C_adv[t]*(self.teta_eq-self.T_boundary[t])\
                             )
            
            self.vettore_prova_2[t] = (
                             self.C_surf
                             )
            
            self.vettore_prova_3[t] = (self.C_adv[t])
            self.vettore_prova_4[t] = (self.teta_eq)
            
            
     
    def variante_u_exc(self,wind_w,H_u,T_u):
        """
        Internal method for the u_exchange calculation
        Method from UWG
        
        Parameters
            ----------
            wind_w : float
                wind velocity [m/s]
            H_u : float
                Urban heat flux [W]
            T_u : float
                Urban Temeprature [K]
            
                
        Returns
        -------
        float : u exchange [m/s]
        float : u friction [m/s]        
        """
    
    
        z0r = 1.
        zb = self.h_aver_buil*2.
        z_ref = 150.
        zm = 10.
        lambda_f = self.VH_urb/4.
        if lambda_f <0.15:
            zou = lambda_f*self.h_aver_buil
        else:
            zou = 0.15*self.h_aver_buil
        
        if lambda_f <0.05:
            du = 3.*lambda_f*self.h_aver_buil
        elif lambda_f>=0.15:
            du =(0.7 + 0.35*(lambda_f-0.15))*self.h_aver_buil
        else:
            du =(0.15 + 5.5*(lambda_f-0.05))*self.h_aver_buil
                
        
        
        A = np.log(z_ref / z0r)
        B = np.log(zm / z0r)
        C = np.log(zb / zou)
        D = np.log(z_ref / zou)
        
        u_b = wind_w*A*C/(B*D)
        u_star_u = 0.4*u_b/np.log((zb - du)/zou)
        
        w_star = np.float64((9.8*H_u*z_ref)/(self.ro_air*self.cp_air*T_u*1.))**1/3
        
        u_star_f = np.max([u_star_u,w_star])
        
        return u_star_f*0.3, u_star_f 
            
            
             
#%% ---------------------------------------------------------------------------------------------------
#%% Some functions to test the urban climate


if False:
    """ Main Climate """
    
    """PRE-PROCESSING"""
    
    "Weather input data"
    Q_rad = np.array([0,0,0,0,0,29.164,101.527,200.029,324.238,463.25,593.808,683.925,714.774,675.527,577.158,441.681,299.915,174.034,75.6694,9.57349,0,0,0,0])  # [W/m2] 
    Q_dir = np.array([0,0,0,0,0,0.164042658,13.52728513,54.02875899,125.2376782,219.2499998,315.8084683,386.9251599,412.773546,384.5272029,312.1579894,215.6806761,122.9145198,53.03385169,13.6694451,0.57349043,0,0,0,0])
    Q_diff = Q_rad-Q_dir
    T_air_rural = 273.15 + np.array([23.4,23,22.5,21.8,21.7,21.7,22,23.2,24.9,27.2,28.6,29.8,30.8,31,30.8,30.1,29.6,29,28.2,27.4,26.5,25.6,24.9,24.3]) # [K]
    T_sky = 273.15 + np.array([18.422,16.1685,18.1232,16.3289,13.1646,16.8236,15.9892,20.144,20.8652,20.2225,17.6048,22.0515,20.557,28.2431,25.455,23.8089,25.0776,19.7162,22.4727,21.825,22.9978,21.5028,23.2381,15.2739])  # [K] 
    lamda_zenith = np.pi/180*np.array([116.4252579,115.8075354,112.3967515,106.6031521,98.97569925,90.03934276,80.23230418,69.91648077,59.42468648,49.13937191,39.63628388,31.95942971,27.89448333,29.09130397,34.99576857,43.65019608,53.59634162,64.03904829,74.51213152,84.66814494,94.16988148,102.6267237,109.5599514,114.4175936]) #[rad]
    u_atm_input = np.array([1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,2.0,3.0,4.0,5.0,5.0,5.0,4.0,3.0,1.0,1.0,1.0,1.0,1.0,0.0]) # [m/s]
    
    "Soil input data"
    d_soil = 3                                                                     # [m]
    c_v_soil = 1.4*10**6                                                           # [J/m3*K]
    
    "Road input data"
    d_asphalt = 0.05                                                               # [m]                                                  
    d_stones = 0.2                                                                 # [m]
    d_gravel = 1                                                                   # [m]
    c_v_asphalt = 1.9*10**6                                                        # [J/m3*K]
    c_v_stones = 2*10**6                                                           # [J/m3*K]
    c_v_gravel = 1.4*10**6                                                         # [J/m3*K]
    
    "City input data"
    H_waste = 0.1                                                                  # [W/K]
    H_traffic = 5                                                                  # [W/K]
    H_traffic_weekday = H_traffic*np.array([0.2,0.2,0.2,0.2,0.2,0.4,0.7,0.9,0.9,0.6,0.6,0.6,0.6,0.6,0.7,0.8,0.9,0.9,0.8,0.8,0.7,0.3,0.2,0.2])
    VH_urb = 0.48                                                   
    horiz_buil_densi = 0.54
    horiz_veg_densi = 0.16
    h_aver_buil = 14.6                                                             # [m]                 
    D_city = 1000                                                                  # [m]
    A_w = 10**6  
    T_w = T_air_rural                                                                            
    T_in = 26+273.15                                                               # [K]
    emissivity_w = 0.85                                       
    riflex_w = 0.15                                              
    A_win = A_w/7                                                                  # [m2]
    U_win = 2                                                                      # [W/(m2*K)]
    h_rad_w = 5                                                                    # [W/(m2*K)]
    h_rad_r = 5                                                                    # [W/(m2*K)] 
    #h_rd_sky = 5                                                                   # [w/(m2*K)] =4*emissivita_vapore*1* 5.67*10^-8*(T_medio_can,j)^3                                               
    V_inf = 10**(-4)                                                               # [m3/s]
    altezza_layer_inf = 2                                                          # [m]
    altezza_layer_sup = 150                                                        # [m]
    zi_day = 1000                                                                  # [m]
    zi_night = 50                                                                  # [m]
    z_m = 10                                                                       # [m]
    z_ref = 150                                                                    # [m]
    
    
    """SIMULATION"""
    
    "Padova Canyon inizialization"
    Padova_Canyon = UrbanCanyon(d_soil,c_v_soil,d_asphalt,d_stones,d_gravel,c_v_asphalt,c_v_stones,c_v_gravel,A_w,emissivity_w,riflex_w,h_rad_w,h_rad_r,U_win,V_inf,VH_urb,horiz_buil_densi,horiz_veg_densi,h_aver_buil,D_city,altezza_layer_inf,altezza_layer_sup,zi_day,zi_night,z_m,z_ref)
    
    for t in range (24):
        Padova_Canyon.solve_canyon(t,u_atm_input[t],Q_rad[t],T_w[t],T_in,T_air_rural[t],T_sky[t],H_waste,H_traffic_weekday[t],lamda_zenith[t],Q_diff[t],Q_dir[t])
    
    
    """POST-PROCESSING"""
    
    Temperatures = np.array([T_air_rural,Padova_Canyon.T_boundary,Padova_Canyon.T_urb]).transpose()
    Coeff_H = np.array([Padova_Canyon.H_rur,Padova_Canyon.H_urb]).transpose()
    Temp_surfaces = np.array([Padova_Canyon.soil.T_suolo_prova,Padova_Canyon.road.T_r_prova]).transpose()

if __name__ == '__main__':

        
    import sys
    #import deepcopy as dpc
    import pandas as pd
    import os
    import pvlib
    import numpy as np
    import time as tm
    from IrradiancePreProcessor import PlanesIrradiances
    from Envelope import loadEnvelopes
    from EndUse import loadArchetype, loadSimpleArchetype
    from CityJSON import JsonCity
    from WeatherData import TskyCalc, SolarPosition, rescale_weather, rescale_sol_gain
    from BuildingsPlants import loadPlants
    import matplotlib.pyplot as plt
    
    #%%
    '''Setting Input Data'''
    # PREPROCESSING
    iopath = os.path.join('..','Input', 'ARPAV_Padova_Legnaro_2019_v2.epw')      # Epw directory from home folder
    year = 2020                                                                    # Setting an year
    first_day = 7                                                                  # Setting the first day of the year (1 Monday)
    tz='Europe/Rome'                                                               # Setting the thermal zone
    azSubdiv = 8                                                                   # Azimuth subdivision
    hSubdiv = 3                                                                    # Tilt subdivision
    SolarCalc = False                                                              # Need of a Solar Calculation?
    env_path = os.path.join('..','Input','B_env_580.xlsx')                          # Envelope directory from home folder
    schedpath = os.path.join('..','Input','schedule580ed.xlsx')                     # Annual schedules
    SchedMethod = str('A')                                                         # A Daily schedules, B yearly schedules
    plant_path = os.path.join('..','Input','PlantsList.xlsx')                       # Plants list and properties
    
    # SIMULATION
    ts = 1                                                                         # Timestep per hour
    hours = 8760                                                                   # Hours of simulation
    years = 1                                                                      # Years of simulation
    model = str('2C')                                                              # Select model: '1C' or '2C'
    mode = str('geojson')                                                         # Select mode: 'geojson' or 'cityjson'
    jsonfile = str('PiovegoRC_hcplant1.geojson')                                                # Select .geojson or .json file
    path=os.path.join('.','Input', jsonfile)
    Shading_calc = str('YES')                                                       # Select 'YES' or 'NO' to take into consideration the shading effect
    toll_az = float(80)                                                            # Semi-tollerance on azimuth of shading surface [°]
    toll_dist = float(100)                                                         # Tollerance on distance of shading surface [m]
    toll_theta = float(80)                                                         # Semi-tollerance on position of the shading surfaces [°]
    R_f = float(0)                                                                 # Reduction factor of the direct solar radiation due to the shading effect [0-1]
    DD_boundaries = np.array([[167,504],[4681,5017]], dtype = int)                 # Heating and Cooling Design Days Periods
    Time_to_regime = 168                                                           # Time needed to reach a regime condition for Design Days Calculation
    DesignDays_calc = str('YES')                                                   # Select 'YES' or 'NO' to calculate or not design days demand
    Plant_calc = str('YES')                                                        # Select 'YES' or 'NO' to calculate or not buildings plant
    
    # OUTPUT REPORT
    OutRep = bool(True)
    
    # Urban Canyon Data
    UWG_calc = bool(True) 
    UWG_data = {
        'Area' : float(400*400*3.14),
        'Diameter': float(400*2),
        'Perimeter': float(400*3.14),
        'Ortogonal_length': float(400*2),
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
        'Vegetation_density': float(0.24),
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
        'vehicle_multiplying_factor': float(0),
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
    
    '''Importing and processing weather data from .epw'''
    epw = pvlib.iotools.read_epw(iopath, coerce_year = year)                       # Reading the epw via pvlib 
    epw_res = epw[0].reset_index(drop=True)                                        # Exporting the hourly values
    lat, lon = epw[1]['latitude'], epw[1]['longitude']                             # Extracting latitude and longitude from the epw
    site = pvlib.location.Location(lat, lon, tz = tz)                              # Creating a location variable
    time = np.arange(8760)                                                         # Time vector inizialization
    
    '''Weather Data and Average temperature difference between Text and Tsky'''
    w = epw_res['wind_speed']                                                      # [m/s]
    T_ext =epw_res['temp_air']                                                     # [°C]
    T_ext_H_avg = np.mean([T_ext[0:2160],T_ext[6599:8759]])                        # [°C]
    RH_ext = epw_res['relative_humidity']/100                                      # [0-1]
    T_dp = epw_res['temp_dew']                                                     # [°C]
    P_ = epw_res['atmospheric_pressure']                                           # [Pa]
    n_opaque = epw_res['opaque_sky_cover']                                         # [0-10]
    dT_er = TskyCalc(T_ext,T_dp,P_,n_opaque)                                       # Average temperature difference between Text and Tsky
    w,T_ext,RH_ext,T_dp,P_,n_opaque = rescale_weather(ts,w,T_ext,RH_ext,T_dp,P_,n_opaque)
    Solar_position = site.get_solarposition(times=epw[0].index).reset_index(drop=True)
    Solar_position = SolarPosition(ts,Solar_position)
    
    '''Inizialization of Solar Gain DataFrame'''
    '''Creates a dataframe with the hourly solar radiation ond angle of incidence for several directions'''
    
    if SolarCalc:
        Irradiances=PlanesIrradiances(site,epw,year,azSubdiv,hSubdiv)
        Irradiances.Irradiances.to_csv(os.path.join('..','Input','PlanesIrradiances.csv'))
    
    '''Solar Gain DataFrame'''
    Solar_Gains = pd.read_csv(os.path.join('..','Input','PlanesIrradiances.csv'),header=[0,1,2],index_col=[0])
    Solar_Gains = rescale_sol_gain(ts,Solar_Gains)
    
    # '''Loading Envelope and Schedule Data'''
    # envelopes = loadEnvelopes(env_path)                                            # Envelope file loading
    
    # if SchedMethod == 'A':
    #     PlantDays=[2520,3984,6192,6912]                                            # 15th April, 15th June, 15th September, 15th October
    #     sched = loadSimpleArchetype(schedpath,time,year,PlantDays,first_day,ts)
    # elif SchedMethod == 'B':
    #     sched = loadArchetype(schedpath,time,ts)
    # else:
    #     sys.exit('Set a proper schedule inporting methodology')
    
    # '''Plants List'''
    # Plants_list = loadPlants(plant_path)
    
    end = tm.time()
    print('Pre-processing:       ', end - start)
    
    UWG_data['tot_wall_opaque_area'] = 1000000
    UWG_data['tot_wall_glazed_area'] = 1000000/10*3
    UWG_data['tot_footprint'] = 30000
    UWG_data['h_builidng_average'] = 10.04
    UWG_data['VH_urb'] =0.48
    UWG_data['Buildings_density'] = 0.32
    urban_canyon = UrbanCanyon([ts,8760],UWG_data,os.path.join('..','RC_classes'))
    '''
    Q_rad = np.array([0,0,0,0,0,29.164,101.527,200.029,324.238,463.25,593.808,683.925,714.774,675.527,577.158,441.681,299.915,174.034,75.6694,9.57349,0,0,0,0])  # [W/m2] 
    Q_dir = np.array([0,0,0,0,0,0.164042658,13.52728513,54.02875899,125.2376782,219.2499998,315.8084683,386.9251599,412.773546,384.5272029,312.1579894,215.6806761,122.9145198,53.03385169,13.6694451,0.57349043,0,0,0,0])
    Q_diff = Q_rad-Q_dir
    T_ext = np.array([23.4,23,22.5,21.8,21.7,21.7,22,23.2,24.9,27.2,28.6,29.8,30.8,31,30.8,30.1,29.6,29,28.2,27.4,26.5,25.6,24.9,24.3]) # [K]
    T_sky = np.array([18.422,16.1685,18.1232,16.3289,13.1646,16.8236,15.9892,20.144,20.8652,20.2225,17.6048,22.0515,20.557,28.2431,25.455,23.8089,25.0776,19.7162,22.4727,21.825,22.9978,21.5028,23.2381,15.2739])  # [K] 
    lamda_zenith = np.array([116.4252579,115.8075354,112.3967515,106.6031521,98.97569925,90.03934276,80.23230418,69.91648077,59.42468648,49.13937191,39.63628388,31.95942971,27.89448333,29.09130397,34.99576857,43.65019608,53.59634162,64.03904829,74.51213152,84.66814494,94.16988148,102.6267237,109.5599514,114.4175936]) #[rad]
    u_atm_input = np.array([1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,2.0,3.0,4.0,5.0,5.0,5.0,4.0,3.0,1.0,1.0,1.0,1.0,1.0,0.0]) # [m/s]
    '''
    
    for t in time:
            T_e = T_ext[t]
            RH_e = RH_ext[t]
            
            if T_e < 0:
                p_extsat = 610.5*np.exp((21.875*T_e)/(265.5+T_e))
            else:
                p_extsat = 610.5*np.exp((17.269*T_e)/(237.3+T_e))
        
            radiation = [Solar_Gains['0.0','0.0','global'].iloc[t],
                             Solar_Gains['0.0','0.0','direct'].iloc[t]]
            
            #radiation = [Q_rad[t],Q_dir[t]]
            
            
            T_w = T_e
            
            H_waste = 0.001*UWG_data['Area']
            
            urban_canyon.solve_canyon(t,T_e,w[t],radiation,T_w,T_e-dT_er,H_waste ,Solar_position.zenith[t],[26,26],[0,0])
            #urban_canyon.solve_canyon(t,T_e,u_atm_input[t],radiation,T_w,T_sky[t],H_waste ,lamda_zenith[t],[26,26],[0,0])
            T_e = urban_canyon.T_urb[t] - 273.15
            
    A = T_ext
    B = urban_canyon.T_urb - 273.15
    
    fig1, [ax1, ax2] = plt.subplots(nrows= 2)
    
    UHI = B-A
    media2 = UHI[3627:6554]  
    

    #ax2.plot(urban_canyon.vettore_prova_4[4810:4890]-273.15)
    #ax.plot(urban_canyon.u_circ[4810:4890])
    #ax.plot(urban_canyon.T_boundary[4810:4890]-273.15)
    