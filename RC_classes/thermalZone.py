'''IMPORTING MODULES'''

import sys
import statistics
import os
from copy import deepcopy as cp
from warnings import warn as wrn
import pvlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate
from RC_classes.Geometry import  Surface, SurfaceInternalMass, SurfaceInternalAdjacent
from RC_classes.Envelope import loadEnvelopes, Envelope
from RC_classes.EndUse import loadArchetype, Archetype
from RC_classes.AHU import AirHandlingUnit
from RC_classes.PVPlant import DistrictPVGIS
from RC_classes.BuildingsPlants import Plants
from RC_classes.WeatherData import  Weather
#from mpl_toolkits.mplot3d import Axes3D

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions

def impedenceParallel(R,C,T_RA = 5.):
    '''
    equivComplexRes Given two vectors (thermal resistances r and thermal 
    capacitances c) of length m (number of walls of the same type, ie either 
    IW or AW), calculates the equivalent complex thermal resistance Zeq
    according to T_RA (period in days)
    
    Parameters
    ----------
        R : np.array
            np.array with the surfaces resistances
        C: np.array
            np.array with the surfaces capacitances
        T_RA: float
            Reference time (number of days)
        
    Returns
    -------
    tuple of floats: the equivalent resistance anc capacitance 
    '''    
    
    # Check input data type
    
    if not isinstance(R, np.ndarray):
        raise TypeError(f'Ops... impedenceParallel function, input R is not a np.array: R {R}') 
    if not isinstance(C, np.ndarray):
        raise TypeError(f'Ops... impedenceParallel function, input C is not a np.array: C {C}') 
    if not isinstance(T_RA, float):
        raise TypeError(f'Ops... impedenceParallel function, input T_RA is not a float: T_RA {T_RA}') 
        
    # T_RA = 5;       % Period = 5 days
    omega_RA = 2*np.pi/(86400*T_RA)
    z = np.zeros(len(R),complex)
    z = (R+1j/(omega_RA*C))   # vettore delle Z 
    
    Z1eq = 1/(sum(1/z))                # Z equivalente
    
    R1eq = np.real(Z1eq)
    C1eq = 1/(omega_RA*np.imag(Z1eq))
    
    # Check output quality
    
    if R1eq < 0. or C1eq <0.:
        wrn(f"\n\nimpedenceParallel funtion, There's something wrong with the calculation, R or C is negative.. R {R}, C {C}\n")
    
    return R1eq, C1eq


def tri2star(T1,T2,T3):
    '''
    #tri2star Transforms three resistances in triangular connection into
    #three resistances in star connection
    
    Parameters
    ----------
        T1 : float
            Resistance 1
        T2: float
            Resistance 2
        T3: float
            Resistance 3
        
    Returns
    -------
    tuple of floats: 3 new resistances (star connection)
    '''
    
    # Check input data type
    
    if not isinstance(T1, float):
        raise TypeError(f'Ops... tri2star function, input T1 is not a float: T1 {T1}') 
    if not isinstance(T2, float):
        raise TypeError(f'Ops... tri2star function, input T2 is not a float: T2 {T2}') 
    if not isinstance(T3, float):
        raise TypeError(f'Ops... tri2star function, input T3 is not a float: T3 {T3}')

    # Check input data quality

    if T1 < 0. or T2 < 0. or T3 < 0.:
        wrn(f"\n\ntri2star funtion, There's something wrong with the input, one of them is negative.. T1 {T1}, T2 {T2}, T3 {T3}\n")
    
    T_sum = T1+T2+T3
    S1 = T2*T3/T_sum
    S2 = T1*T3/T_sum
    S3 = T2*T1/T_sum
    return S1, S2, S3


def longWaveRadiation(theta_a,SSW = 1.):
    '''
    Estimation of sky and ground temperatures via vdi6007 model:
        theta_a outdoor air temperature [°C]
        SSW factor to count the clear non-clear sky
        
    Parameters
    ----------
        theta_a : np.array
            external temeprature [°C]
        SSW : float
            factor to count the clear non-clear sky renge 0-1
        
    Returns
    -------
    tuple of np.array: 4 arrays, irradiance from sky vault,
                                irradiance from ground,
                                ground equivalent temeprature,
                                sky equivalent temeperature
    '''

    # Check input data type

    if not isinstance(theta_a, np.ndarray):
        raise TypeError(f'Ops... longWaveRadiation function, input theta_a is not a np.array: theta_a {theta_a}') 
    if not isinstance(SSW, float):
        try:
            SSW = float(SSW)
        except ValueError:            
            raise TypeError(f'Ops... longWaveRadiation function, input SSW is not a float: SSW {SSW}') 

    # Check input data quality
    
    if not np.all(np.greater(theta_a,-50.)) or not np.all(np.less(theta_a,60.)):
        wrn(f"\n\nlongWaveRadiation function, the theta_a input is outside range [-50,60]: theta_a {theta_a}")
    if not 0. <= SSW <= 1.:
        wrn(f"\n\nlongWaveRadiation function, the SSW input is outside range [0,1]: SSW {SSW}")

    Ea_1 = 9.9*5.671*10**(-14)*(273.15+theta_a)**6
    
    alpha_L = 2.30 - 7.37*10**(-3)*(273.15+theta_a)
    alpha_M = 2.48 - 8.23*10**(-3)*(273.15+theta_a)
    alpha_H = 2.89 - 1.00*10**(-2)*(273.15+theta_a)
    
    Ea = Ea_1*(1+(alpha_L+(1-(1-SSW)/3)*alpha_M+((1-(1-SSW)/3)**2)*alpha_H)*((1-SSW)/3)**2.5)
    Ee = -(0.93*5.671*10**(-8)*(273.15+theta_a)**4+(1-0.93)*Ea)
    
    theta_erd = ((-Ee/(0.93*5.67))**0.25)*100 - 273.15  # [°C]
    theta_atm = ((Ea/(0.93*5.67))**0.25)*100 - 273.15   # [°C]  

    return Ea, Ee, theta_erd, theta_atm


def loadHK(perc_rad, perc_rad_aw, perc_altro_irr, A_aw, A_raum):
    '''
    loadHK  -  Distribution of heat load on the nodes (surface nodes and air node) based on heat emitters
    of the building
    
    Parameters
    ----------
        perc_rad (sigma_fhk): float
            percentage of heat flow by radiant floors
        perc_rad_aw (sigma_fhk_aw): float
            percentage of radiant floors installed on AW walls 
        perc_altro_irr (sigma_rad_str): float
            percentage of radiant load (out of total which is rad+conv) by other types of heat emitters (e.g.: fan-coils, radiators)
        A_aw : float
            sum of the exterior opaque building components
        A_raum : float
            sum of internal partitions (all IW components) and exterior opaque building components
    
    Returns
    tuple of float:
        sigma_hk_iw: radiant heat load  on IW building components (on surface node IW)
        sigma_hk_aw: radiant heat load  on AW building components (on surface node AW)
        sigma_hk_kon: convective heat load (on air node)
    '''
    
    # Check input data type
    
    if not isinstance(perc_rad, float):
        try:
            perc_rad = float(perc_rad)
        except ValueError:            
            raise TypeError(f'Ops... loadHK function, input perc_rad is not a float: perc_rad {perc_rad}') 
    if not isinstance(perc_rad_aw, float):
        try:
            perc_rad_aw = float(perc_rad_aw)
        except ValueError:            
            raise TypeError(f'Ops... loadHK function, input perc_rad_aw is not a float: perc_rad_aw {perc_rad_aw}') 
    if not isinstance(perc_altro_irr, float):
        try:
            perc_altro_irr = float(perc_altro_irr)
        except ValueError:            
            raise TypeError(f'Ops... loadHK function, input perc_altro_irr is not a float: perc_altro_irr {perc_altro_irr}') 
    if not isinstance(A_raum, float):
        try:
            A_raum = float(A_raum)
        except ValueError:            
            raise TypeError(f'Ops... loadHK function, input A_raum is not a float: A_raum {A_raum}') 
    if not isinstance(A_aw, float):
        try:
            A_aw = float(A_aw)
        except ValueError:            
            raise TypeError(f'Ops... loadHK function, input A_aw is not a float: A_aw {A_aw}') 
    
    # Check input data quality
    
    for p in [perc_rad,perc_rad_aw,perc_altro_irr]:
        if not 0. <= p <= 1.:
            wrn(f"\n\nloadHK function, one of the percentage inputs is outside range [0,1]: perc_rad {perc_rad}, perc_rad_aw {perc_rad_aw}, perc_altro_irr {perc_altro_irr}")
    for area in [A_raum,A_aw]:
        if not 0. <= area :
            wrn(f"\n\nloadHK function, one of the area inputs is negative: A_raum {A_raum}, A_aw {A_aw}")    
    
    # %Note: the sum of the 3 outputs must be equal to 1
    
    if perc_rad == 1:

        perc_rad_iw = 1 - perc_rad_aw
        sigma_hk_iw = perc_rad_iw
        sigma_hk_aw = perc_rad_aw
        sigma_hk_kon = 0

    elif perc_rad == 0:

        sigma_hk_iw = 0
        sigma_hk_aw = 0
        sigma_hk_kon = 1
    
    else:
        perc_rad_iw = 1 - perc_rad_aw
        perc_altro = 1 - perc_rad
        perc_altro_irr = perc_altro*perc_altro_irr
        sigma_hk_iw = perc_altro_irr*(A_raum-A_aw)/A_raum + perc_rad_iw
        sigma_hk_aw = perc_altro_irr*A_aw/A_raum + perc_rad_aw
        sigma_hk_kon = 1 - sigma_hk_iw - sigma_hk_aw
        
    return [sigma_hk_iw, sigma_hk_aw, sigma_hk_kon]

#%% ---------------------------------------------------------------------------------------------------
#%% ThermalZone class

class ThermalZone:
    '''
    thermal zone class
    
    __init__:
        zone number
        name of the building
        envelope object of the zone
        schedule archetype object of the zone
        list of the surfaces
        volume
        zone area
        
    zoneParameter13790 calculates 1C params. No input
    
    calculate_zone_loads_ISO13790 calculates 1C zone loads:
        Solar Gains vector [W/m2]
        average difference between t air and t sky [°C]
        .
        .       
        
    Sensible1C 1C system solver:
        solver typology phi set or T set
        vector with ventilation and infiltration heat transfer coefficients [W/K]
        Setpoint temperature [°C]
        supply setpoint temperature [°C]
        external temperature [°C]
        time step duration [s]
        plant power setpoint [W]
        
    zoneParameter6007 calculates 2C params. No input
    
    calculate_zone_loads_vdi6007 calculate 2C zone loads:
        external temeprature [°C]
        Solar Gains dataframe [W/m2]
        
    Sensible2C 2C system solver:
        solver typology phi set or T set
        gains vector [W]
        vector with ventilation and infiltration heat transfer coefficients [W/K]
        external and sun-air temperature [°C]
        supply setpoint temperature [°C]
        time step [s]
        Setpoint temperature [°C]
        plant power setpoint [W]

    solveTZ thermal zone balances solver:
        time step [-]
        external temperature [°C]
        external relative humidity [-]
        external saturation pressure [Pa]
        time step duration [s]
        mode 2C or 1C
        
    Methods list:
        init
        zoneParameter13790
        calculate_zone_loads_ISO13790
        Sensible1C
        zoneParameter6007
        calculate_zone_loads_vdi6007
        Sensible2C
        solveDD_Heating
        solveTZ
        reset_init_values
    '''
    
    # Class attributes
    
    htr_ms = 9.1                                                               # heat tranfer coeff. ISO 13790 [W/(m2 K)]
    his = 3.45                                                                 # heat tranfer coeff. ISO 13790 [W/(m2 K)]
    rho_air = 1.2                                                              # [kg/m3]
    cp_air = 1000.                                                              # [J/(kg K)]
    sens_frac = 0.5462                                                         # Based on Standard ISO 18523
    p_atm = 101325.                                                             # [Pa]
    r_0 = 2501000.                                                              # [J/kg]
    cpv = 1875.                                                                 # [J/(kg K)]
    x_m0 = 0.0105                                                              # [kg_v/kg_as]
    x_m0DD = 0.0105                                                            # [kg_v/kg_as]
    T_wall_0 =  15                                                              # [°C]
    
    conv_people = 1.                                                            # people convective fraction (with respect to convective and radiant)
    conv_app = 1.                                                               # appliances convective fraction (with respect to convective and radiant)
    
    sigma_fhk = 0         # percentage of heating from from radiant floor (vdi 6007)
    sigma_fhk_aw = 0      # percentage of radiant floor systems embedded in external (AW) walls (vdi 6007) 
    sigma_hk_str = 1      # radiant fraction of other terminal units (not embedded in envelopes) 50% rad, 50% conv
    
    def __init__(self,zoneNumber,bd_name,envelope,sched_db,surfList,volume,zone_area,l):
        '''
        Initializes the thermal zone
        creates some np.array which will contain the zone parameters and variables
        
        Parameters
            ----------
            zoneNumber : int
                the number id of the thermal zone
            envelope : Envelope object
                from the dictionary that returns from loadEnvelopes
            sched_db : Archetype object
                from the dictionary that returns from loadSchedules
            surfList : list or tuple
                list or tuple of Surface objects
            volume : float
                the zone volume [m3]
            zone_area : float      
                the zone area [m2]
            l : int 
                number of simulation time steps
                    
        Returns
        -------
        None.        
        ''' 

        # Check input data type

        if not isinstance(zoneNumber, int):
            try:
                zoneNumber = int(zoneNumber)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone initialization, bd {bd_name}, zoneNumber input is not an int: zoneNumber {zoneNumber}')
        if not isinstance(l, int):
            try:
                l = int(l)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone initialization, bd {bd_name}, l input is not an int: l {l}')
        if not isinstance(volume, float):
            try:
                volume = float(volume)
            except ValueError:
                raise TypeError(f'Ops... thermal zone initialization, bd {bd_name}, volume input is not an float: volume {volume}')
        if not isinstance(zone_area, float):
            try:
                zone_area = float(zone_area)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone initialization, bd {bd_name}, volume input is not an float: zone_area {zone_area}')
        if not isinstance(envelope, Envelope):
            raise TypeError(f'Ops... thermal zone initialization, bd {bd_name}, envelope input is not an Envelope object: envelope {envelope}')
        if not isinstance(sched_db, Archetype):
            raise TypeError(f'Ops... thermal zone initialization, bd {bd_name}, envelope input is not an Archetype object: envelope {envelope}')
                 
        # Check input data quality
        
        if volume < 0.:
            wrn(f"\n\n Thermal zone initialization, bd {bd_name}, the zone volume is negative: volume {volume}. It will be set to 1.")
            volume = 1.
        if zone_area < 0.:
            wrn(f"\n\n Thermal zone initialization, bd {bd_name}, the zone area is negative: zone_area {zone_area}. It will be set to 1.")
            zone_area = 1.
        
        # Sets some attributes of the zone
        
        self.bd_name = bd_name
        self.name ='Zone ' + str(zoneNumber)
        self.schedules = sched_db
        self.Strat={'ExtWall':envelope.ExtWall,\
        'IntWall':envelope.IntWall,\
        'IntCeiling':envelope.IntCeiling,\
        'GroundFloor':envelope.GroundFloor,\
        'Roof':envelope.Roof,
        'Window':envelope.Window,
        'IntFloor':envelope.IntFloor}                            
        self.V = volume
        self.Ca = self.V*self.rho_air*self.cp_air 
        self.surfaces = {}
        self.zone_area = zone_area
        self.Ta0 = 15
        self.theta_m0 = 15
        self.theta_m0_vdi = [15, 15]                                           # aw the first, iw the second

        # Calculates the external wall area and other geometrical params

        i = 0
        self.Araum = 0
        self.Aaw = 0
        self.Tot_glazed_area = 0
        self.Tot_opaque_area = 0
        for surface in surfList:
            if not (isinstance(surface,Surface) or isinstance(surface,SurfaceInternalMass) or isinstance(surface,SurfaceInternalAdjacent)):
                raise TypeError(f'Ops... thermal zone initialization, surfList should contain all Surface objects: surfList {surfList}')
            i += 1
            self.surfaces[('Surface '+str(i))]=surface
            self.Araum += self.surfaces[('Surface '+str(i))].area
            if surface.type == 'ExtWall' or surface.type == 'GroundFloor' or surface.type == 'Roof' :
                self.Aaw += self.surfaces[('Surface '+str(i))].area
            if surface.type == 'ExtWall':
                self.Tot_glazed_area += surface.glazedArea
                self.Tot_opaque_area += surface.opaqueArea
                
        # Vectors inizialization 
        
        self.heatFlow = np.zeros(l)
        self.Air_temp = np.zeros(l)
        self.theta_m0_2c = np.zeros([l,2])
        self.ext_surf_heat_gain = np.zeros(l)                                                                     
        self.RH_i = np.zeros(l)
        self.latentFlow = np.zeros(l)
        self.x_ext = np.zeros(l)
        self.G_v = np.zeros(l)
        self.G_da_inf = np.zeros(l)
        self.G_da_vent = np.zeros(l)
        self.x_int = np.zeros(l)
        self.T_sup = np.zeros(l)
        self.ZoneAHU = AirHandlingUnit(l)
        self.T_wall_0_vector = np.zeros(l)
        self.R_lim_ext_wall_2C = self.Strat['ExtWall'].R_se / (self.Tot_opaque_area + self.Tot_glazed_area)
        self.H_lim_ext_wall_1C = self.Tot_opaque_area / self.Strat['ExtWall'].R_se
        
        # Parameters init
        
        self.Htr_is = 0.
        self.Htr_w = 0.
        self.Htr_ms = 0.
        self.Htr_em = 0.
        self.Cm = 0.
        self.DenAm = 0.
        self.Atot = 0.
        self.Htr_op = 0.
        
        self.RrestAW = 0.
        self.R1AW = 0.
        self.RalphaStarAW = 0.
        self.RalphaStarIL = 0.
        self.RalphaStarIW = 0.
        self.R1IW = 0.
        self.C1AW = 0.
        self.C1IW = 0.              
        
    def zoneParameter13790(self):
        '''
        Calculates the thermal zone parameters of the ISO 13790
        it does not require input
        
        Parameters
            ----------
            None
                    
        Returns
        -------
        None.        
        ''' 
        
        self.Htr_int = []
        
        # list all surface to extract the window and opeque area and other thermo physical prop
        
        for surface in self.surfaces.values():
            if surface.type == 'ExtWall' or surface.type == 'IntWall' or surface.type == 'IntCeiling' or surface.type == 'Roof' or surface.type == 'GroundFloor':
                self.Cm += surface.opaqueArea*(self.Strat[surface.type].k_int)
                self.DenAm += surface.opaqueArea*(self.Strat[surface.type].k_int)**2
            if surface.type == 'IntFloor':
                self.Cm += surface.opaqueArea*(self.Strat[surface.type].k_est)
                self.DenAm += surface.opaqueArea*(self.Strat[surface.type].k_est)**2
            self.Atot += surface.area
            
            if surface.type == 'ExtWall' or surface.type == 'GroundFloor' or surface.type == 'Roof' :
                self.Htr_op += surface.opaqueArea*(self.Strat[surface.type].U)
                self.Htr_w += surface.glazedArea*(self.Strat['Window'].U)
                
            if type(surface) == SurfaceInternalAdjacent:
                self.Htr_int.append({'H':surface.opaqueArea*self.Strat[surface.type].U,'AdjacentZone':surface.adjacentZone})
        
        # Final calculation
        
        self.Am = self.Cm**2/self.DenAm
        self.Htr_ms = self.Am*self.htr_ms
        self.Htr_em = 1/(1/self.Htr_op - 1/self.Htr_ms)
        self.Htr_is = self.his*self.Atot
        self.UA_tot = self.Htr_op + self.Htr_w

        
    def calculate_zone_loads_ISO13790(self,weather,h_r):
        '''
        Calculates the heat gains on the three nodes of the ISO 13790 network
        Vectorial calculation
        
        Parameters
            ----------
            weather : RC_classes.WeatherData.weather
                weather obj
            h_r : float    
                external radiative heat transfer coefficient W/(m2 K)
                
        Returns
        -------
        None.        
        ''' 
        
        # Check input data type
        
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
        if not isinstance(weather.dT_er, float):
            try:
                weather.dT_er = float(weather.dT_er)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone calculate_zone_loads_ISO13790, bd {self.bd_name}, dT_er input is not an float: dT_er {weather.dT_er}')
        if not isinstance(h_r, float):
            try:
                h_r = float(h_r)
            except ValueError:
                raise TypeError(f'Ops... thermal zone calculate_zone_loads_ISO13790, bd {self.bd_name}, h_r input is not an float: h_r {h_r}')
    
        # Check input data quality
        
        if h_r < 0.:
            wrn(f"\n\n Thermal zone calculate_zone_loads_ISO13790, bd {self.bd_name}, the h_r is negative: h_r {h_r}.")
        
        if not 0. <= weather.dT_er < 30.:
            wrn(f"\n\n Thermal zone calculate_zone_loads_ISO13790, bd {self.bd_name}, the dT_er is outside limits [0,30] °C: dT_er {weather.dT_er}.")
        
        # First calculation of internal heat gains 
        
        phi_int = (self.schedules.people*self.sens_frac+self.schedules.appliances+self.schedules.lighting)*self.zone_area
        phi_sol_gl_tot = 0
        phi_sol_op_tot = 0
        
        # Solar radiation

        for i in self.surfaces.values():
            
            phi_sol_op = 0
            phi_sol_gl = 0
            
            if i.type == 'ExtWall' or i.type == 'Roof':
                irradiance = weather.SolarGain[str(float(i.azimuth_round))][str(float(i.height_round))]  
                BRV = irradiance['direct'].to_numpy()                           
                TRV = irradiance['global'].to_numpy()
                DRV = TRV -BRV
                A_ww = i.glazedArea
                self.A_ew = i.opaqueArea
                
                # Glazed surfaces 
                F_sh = self.Strat['Window'].F_sh
                F_w = self.Strat['Window'].F_w
                F_f = self.Strat['Window'].F_f
                F_so = self.Strat['Window'].F_so
                AOI = irradiance['AOI'].to_numpy()
                shgc = interpolate.splev(AOI, self.Strat['Window'].SHGC_profile, der=0)
                shgc_diffuse = interpolate.splev(70, self.Strat['Window'].SHGC_profile, der=0)
                if i.OnOff_shading == 'On':
                    phi_sol_gl = F_so*(BRV*F_sh*F_w*(1-F_f)*shgc*A_ww*i.shading_effect)+F_so*DRV*F_sh*F_w*(1-F_f)*shgc_diffuse*A_ww
                else:
                    phi_sol_gl = F_so*(BRV*F_sh*F_w*(1-F_f)*shgc*A_ww)+F_so*DRV*F_sh*F_w*(1-F_f)*shgc_diffuse*A_ww
                
                
            # Opaque surfaces 
            if i.type == 'ExtWall':
                self.F_so_op = F_so
                self.F_r = i.F_r
                self.alpha = self.Strat['ExtWall'].alpha_est
                self.sr_ew = self.Strat['ExtWall'].R_se
                self.U_ew_net = self.Strat['ExtWall'].U_net
                if i.OnOff_shading == 'On':
                    phi_sol_op = self.F_so_op*(BRV*i.shading_effect + DRV)*self.alpha*self.sr_ew*self.U_ew_net*self.A_ew-self.F_r*self.sr_ew*self.U_ew_net*self.A_ew*h_r*weather.dT_er
                else:
                    phi_sol_op = self.F_so_op*TRV*self.alpha*self.sr_ew*self.U_ew_net*self.A_ew-self.F_r*self.sr_ew*self.U_ew_net*self.A_ew*h_r*weather.dT_er
                
            if i.type == 'Roof':
                TRH = weather.SolarGain['0.0']['0.0']['global'].to_numpy()
                self.F_r = i.F_r
                self.alpha = self.Strat['Roof'].alpha_est
                self.sr_rf = self.Strat['Roof'].R_se
                self.U_rf_net = self.Strat['Roof'].U_net
                phi_sol_op = TRH*self.alpha*self.sr_rf*self.U_rf_net*self.A_ew-self.F_r*self.sr_rf*self.U_rf_net*self.A_ew*h_r*weather.dT_er
                
            # Total solar gain     
            phi_sol_gl_tot += phi_sol_gl
            phi_sol_op_tot += phi_sol_op
        phi_sol = phi_sol_gl_tot+phi_sol_op_tot
       
        # Distribute heat gains to temperature nodes 
        self.phi_ia = 0.5*phi_int
        self.phi_st = (1-self.Am/self.Atot-self.Htr_w/(9.1*self.Atot))*(0.5*phi_int + phi_sol)
        self.phi_m = self.Am/self.Atot*(0.5*phi_int + phi_sol)
        
          
    def Sensible1C(self, flag, Hve, T_set, T_sup_AHU, T_e, phi_load, tau, phi_HC_set = 0.):
        '''
        Solves ISO 13790 network for a specific time step
                
        Parameters
            ----------
            flag : string
                string 'Tset' or 'phiset'
            Hve : list of positive floats
                ventilation and infiltration heat tranfer coeff [W/K]
            T_set : float    
                setpoint temperature [°C]
            T_sup_AHU: float    
                ventilation supply temperature [°C]
            T_e: float    
                external temperature [°C]
            phi_load: list of three integers
                the load on the three nodes network (ia, sm, m respectively) [W]
            tau: int
                time constant (3600 s for hourly sim)
            phi_HC_set : float
                thermal power to the ambient [W]            
                
        Returns
        -------
        np.array
            with demand [W], T_air [°C], T_s [°C] and T_m [°C]
        ''' 
        
        # Check input data type
        
        if flag != 'Tset' and flag != 'phiset':
            raise TypeError(f"Ops... thermal zone Sensible1C, bd {self.bd_name}, flag input is not a 'Tset' or 'phiset': flag {flag}")
        if not isinstance(Hve, list) or not isinstance(Hve[0],float) or not isinstance(Hve[1],float):
            try: 
                Hve[0]=float(Hve[0])
                Hve[1]=float(Hve[1])
            except ValueError:
                raise TypeError(f"Ops... thermal zone Sensible1C, bd {self.bd_name}, Hve input is not a list of floats: Hve {Hve}")
        if not isinstance(T_set, float):
            try:
                T_set = float(T_set)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone Sensible1C, bd {self.bd_name}, T_set input is not a float: T_set {T_set}')
        if not isinstance(T_sup_AHU, float):
            try:
                T_sup_AHU = float(T_sup_AHU)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone Sensible1C, bd {self.bd_name}, T_sup_AHU input is not a float: T_sup_AHU {T_sup_AHU}')
        if not isinstance(T_e, float):
            try:
                T_e = float(T_e)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone Sensible1C, bd {self.bd_name}, T_e input is not a float: T_e {T_e}')
        if not isinstance(phi_load, list) or not isinstance(phi_load[0],float) or not isinstance(phi_load[1],float) or not isinstance(phi_load[2],float):
            try: 
                phi_load[0]=float(phi_load[0])
                phi_load[2]=float(phi_load[2])
                phi_load[1]=float(phi_load[1])
            except ValueError:
                raise TypeError(f"Ops... thermal zone Sensible1C, bd {self.bd_name}, phi_load input is not a list of floats: phi_load {phi_load}")
        if not isinstance(tau, int):
            try:
                tau = int(tau)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone Sensible1C, bd {self.bd_name}, tau input is not a float: tau {tau}')
        if not isinstance(phi_HC_set, float):
            try:
                phi_HC_set = float(phi_HC_set)
            except ValueError: 
                raise TypeError(f'Ops... thermal zone Sensible1C, bd {self.bd_name}, phi_HC_set input is not a float: phi_HC_set {phi_HC_set}')
                
        # Check input data quality
        
        if Hve[0] < 0. or Hve[1] < 0.:
            wrn(f"\n\n Thermal zone Sensible1C, bd {self.bd_name}, the Hve is negative: Hve {Hve}.")
        if not -50. <= T_set < 50.:
            wrn(f"\n\n Thermal zone Sensible1C, bd {self.bd_name}, the T_set is outside limits [-50,50] °C: T_set {T_set}.")
        if not -50. <= T_sup_AHU < 50.:
            wrn(f"\n\n Thermal zone Sensible1C, bd {self.bd_name}, the T_sup_AHU is outside limits [-50,50] °C: T_sup_AHU {T_sup_AHU}.")
        if not -50. <= T_e < 60.:
            wrn(f"\n\n Thermal zone Sensible1C, bd {self.bd_name}, the T_e is outside limits [-50,50] °C: T_e {T_e}.")
        
        # Set some data and build up the system   
        
        phi_ia = phi_load[0]
        phi_st = phi_load[1]
        phi_m = phi_load[2]
        Hve_vent = Hve[0]
        Hve_inf = Hve[1]
        Y = np.zeros((3,3))
        q = np.zeros((3))
        
        if flag == 'Tset':
            
            Y[0,0] = 1
            Y[0,1] = self.Htr_is
            Y[1,1] = -(self.Htr_is+self.Htr_w+self.Htr_ms)
            Y[1,2] = self.Htr_ms
            Y[2,1] = self.Htr_ms
            Y[2,2] = -self.Cm/tau-self.Htr_em-self.Htr_ms
        
            q[0] = Hve_inf*(T_set-T_e)+Hve_vent*(T_set-T_sup_AHU)-phi_ia+self.Htr_is*T_set
            q[1] = -self.Htr_is*T_set-phi_st-self.Htr_w*T_e
            q[2] = -self.Htr_em*T_e-phi_m-self.Cm*self.theta_m0/tau
            x = np.linalg.inv(Y).dot(q)
            return np.insert(x,1,T_set)
        
        if flag == 'phiset':
            
            Y[0,0] = -(self.Htr_is+Hve_inf+Hve_vent)
            Y[0,1] = self.Htr_is
            Y[1,0] = self.Htr_is
            Y[1,1] = -(self.Htr_is+self.Htr_w+self.Htr_ms)
            Y[1,2] = self.Htr_ms
            Y[2,1] = self.Htr_ms
            Y[2,2] = -self.Cm/tau-self.Htr_em-self.Htr_ms
            
            q[0] = -phi_HC_set-Hve_inf*T_e-Hve_vent*T_sup_AHU-phi_ia
            q[1] = -phi_st - self.Htr_w*T_e
            q[2] = -self.Htr_em*T_e-phi_m-self.Cm*self.theta_m0/tau
            y = np.linalg.inv(Y).dot(q)
            return np.insert(y,0,phi_HC_set)            
            
    def zoneParameter6007(self):      
        '''
        Calculates the thermal zone parameters of the VDI 6007
        it does not require input
        
        Parameters
            ----------
            None
                    
        Returns
        -------
        None.        
        ''' 

        # Creation of some arrys of variables and parameters
    
        R1IW_m = np.array([])
        C1IW_m = np.array([])
        R_IW = np.array([])
        R1AW_v = np.array([])
        C1AW_v = np.array([])
        R_AW = np.array([])
        R1_AF_v = np.array([])
        HAW_v = np.array([])
        HAF_v = np.array([])
        alphaKonAW = np.array([])
        alphaKonIW = np.array([])
        alphaKonAF = np.array([])
        alphaStr = 5 #vdi Value
        alphaKonA = 20 #vdi value
        RalphaStrAW = np.array([])
        RalphaStrIW = np.array([])
        RalphaStrAF = np.array([])
        AreaAW =  np.array([])
        AreaAF =  np.array([])
        AreaIW =  np.array([])
        
        # Cycling surface to calculates the Resistance and capacitance of the vdi 6007
        
        for surface in self.surfaces.values():
            if surface.type == 'ExtWall' or surface.type == 'GroundFloor' or surface.type == 'Roof':
                surface_R1, surface_C1 = self.Strat[surface.type].vdi6007surfaceParams(surface.opaqueArea,True)
                #R1AW_v.append(surface_R1)
                C1AW_v = np.append(C1AW_v,[surface_C1],axis = 0)
                #R_AW = np.append(R_AW,[sum(self.Strat[surface.type].r)],axis=0)
                # considering glazing component
                R_AF_v = (self.Strat['Window'].Rl_w/surface.glazedArea) #Eq 26
                #R1_AF_v = np.append(R1_AF_v,[R_AF_v/6],axis=0) #Jacopo utilizes a different formula, but this is what I understood from the standard
                # this part is a little different in Jacopo model,
                # However this part calculates opaque R, glazed R and insert the parallel ass wall R
                R1AW_v = np.append(R1AW_v,[1/(1/surface_R1+1/R_AF_v)],axis=0)
                #R1AW_v = np.append(R1AW_v,[1/(1/surface_R1+6/R_AF_v)],axis=0) ALTERNATIVA NORMA
                
                HAW_v = np.append(HAW_v,self.Strat[surface.type].U*surface.opaqueArea)
                HAF_v = np.append(HAF_v,self.Strat['Window'].U*surface.glazedArea)
                alphaKonAW = np.append(alphaKonAW,[surface.opaqueArea*(1/self.Strat[surface.type].R_si-alphaStr)],axis=0)
                alphaKonAF = np.append(alphaKonAF,[surface.glazedArea*(1/self.Strat['Window'].Ri_w-alphaStr)],axis=0)
                                  
                RalphaStrAW = np.append(RalphaStrAW,[1/(surface.opaqueArea*alphaStr)])
                RalphaStrAF = np.append(RalphaStrAF,[1/(surface.glazedArea*alphaStr)])
                
                AreaAW = np.append(AreaAW,surface.opaqueArea)
                AreaAF = np.append(AreaAF,surface.glazedArea)
                
            elif surface.type == 'IntCeiling' or surface.type == 'IntWall' or  surface.type=='IntFloor':
                surface_R1, surface_C1 = self.Strat[surface.type].vdi6007surfaceParams(surface.opaqueArea,False)
                R1IW_m = np.append(R1IW_m,[surface_R1],axis=0)
                C1IW_m = np.append(C1IW_m,[surface_C1],axis =0)
                R_IW = np.append(R_IW,[sum(self.Strat[surface.type].r)],axis=0)
                alphaKonIW = np.append(alphaKonIW,[surface.opaqueArea*(1/self.Strat[surface.type].R_si-alphaStr)],axis=0)
                
                #if surface.opaqueArea*alphaStr == 0:
                #    print(surface.name)
                RalphaStrIW = np.append(RalphaStrIW,[1/(surface.opaqueArea*alphaStr)])
                
                AreaIW = np.append(AreaIW,surface.area)
            else:
                print('Error.. surface type not found')
        
        # Doing the parallel of the impedances
        
        self.R1AW, self.C1AW = impedenceParallel(R1AW_v,C1AW_v)  #eq 22
        self.R1IW, self.C1IW = impedenceParallel(R1IW_m,C1IW_m)
        
        # Final params
        
        self.RgesAW = 1 / (sum(HAW_v)+sum(HAF_v))  #eq 27

        RalphaKonAW = 1/(sum(alphaKonAW)+sum(alphaKonAF))   #scalar
        RalphaKonIW = 1/sum(alphaKonIW)   #scalar 
        
        if sum(AreaAW) <= sum(AreaIW):
            RalphaStrAWIW = 1/(sum(1/RalphaStrAW)+sum(1/RalphaStrAF))  #eq 29
        else:
            RalphaStrAWIW = 1/sum(1/RalphaStrIW)  #eq 31
            
        self.RrestAW = self.RgesAW - self.R1AW - 1/(1/RalphaKonAW + 1/RalphaStrAWIW) #eq 28
        
        RalphaGesAW_A = 1/(alphaKonA*(sum(AreaAF) + sum(AreaAW)))
        
        if self.RgesAW < RalphaGesAW_A:   # this is different from Jacopo's model but equal to the standard
            self.RrestAW = RalphaGesAW_A  #eq 28a
            self.R1AW = self.RgesAW - self.RrestAW - 1/(1/RalphaKonAW + 1/RalphaStrAWIW) #eq 28b
            
            if self.R1AW < 10**(-10):
                self.R1AW = 10**(-10)   # Thresold (only numerical to avoid division by zero)  #eq 28c
                
        self.RalphaStarIL, self.RalphaStarAW, self.RalphaStarIW = tri2star(RalphaStrAWIW,RalphaKonIW,RalphaKonAW)
        self.UA_tot = sum(HAW_v)+sum(HAF_v)   
        self.Htr_op = sum(HAW_v)
        self.Htr_w = sum(HAF_v) 

    def calculate_zone_loads_vdi6007(self,weather):
        '''
        Calculates zone loads for the vdi 6007 standard
        Also the external equivalent temperature
        
        Parameters
        ----------
            weather : RC_classes.WeatherData.weather
                weather obj
               
        Returns
        -------
        None
        '''

        # Check input data type
        
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
        
        # Check input data quality
                
        '''
        Eerd = Solar_gain['0.0']['0.0']['Global']*rho_ground                          #
        Eatm = Solar_gain['0.0']['0.0']['Global']-Solar_gain['0.0']['0.0']['Direct']
        
        T_ext vettore
        
        '''
        
        Eatm, Eerd, theta_erd, theta_atm = longWaveRadiation(weather.Text)
        #fil =  (Eatm + Eerd)*(theta_erd - theta_atm)
        alpha_str_A = 5
        
        # Creates some vectors and set some parameters
        
        T_ext = weather.Text
        theta_eq = np.zeros([len(T_ext),len(self.surfaces)])
        delta_theta_eq_lw = np.zeros([len(T_ext),len(self.surfaces)])
        delta_theta_eq_kw = np.zeros([len(T_ext),len(self.surfaces)])
        theta_eq_w = np.zeros([len(T_ext),len(self.surfaces)])
        frame_factor = 1 - self.Strat['Window'].F_f
        F_sh = self.Strat['Window'].F_sh*self.Strat['Window'].F_so*self.Strat['Window'].F_w
        Q_il_str_A_iw = 0
        Q_il_str_A_aw = 0
        #self.Q_il_str_A = 0 
        
        i = -1
        
        # Lists all surfaces to calculate the irradiance on each one and creates the solar gains
        
        for surface in self.surfaces.values():
            i += 1
            if surface.type == 'ExtWall' or surface.type == 'Roof':
                if surface.OnOff_shading == 'On':
                    shading = surface.shading_effect
                else:
                    shading = 1
                irradiance = weather.SolarGains[str(float(surface.azimuth_round))][str(float(surface.height_round))]
                AOI = irradiance['AOI'].to_numpy()
                BRV = irradiance['direct'].to_numpy()                           
                TRV = irradiance['global'].to_numpy()
                DRV = TRV -BRV
                phi = surface.F_r  #eventualmente si può importare dalla surface
                alpha_a = self.Strat[surface.type].alpha_conv_est + alpha_str_A
                eps_F = 0.9
                #eps_F = self.Strat[surface.type] 
                delta_theta_eq_lw[:,i] = ((theta_erd - T_ext)*(1-phi) + (theta_atm - T_ext)*phi)*(eps_F*alpha_str_A)/(0.93*alpha_a)
                delta_theta_eq_kw[:,i] = (BRV*shading+(TRV-BRV))*self.Strat[surface.type].alpha_est/alpha_a
                theta_eq[:,i] = (T_ext + delta_theta_eq_lw[:,i] + delta_theta_eq_kw[:,i])*self.Strat[surface.type].U*surface.opaqueArea/self.UA_tot
                theta_eq_w[:,i] = (T_ext + delta_theta_eq_lw[:,i])*self.Strat['Window'].U*surface.glazedArea/self.UA_tot
                
                shgc = interpolate.splev(AOI, self.Strat['Window'].SHGC_profile, der=0)
                shgc_diffuse = interpolate.splev(70, self.Strat['Window'].SHGC_profile, der=0)                                                                                                                                                         
                # Jacopo quì usa come A_v l'area finestrata, mentre la norma parla di area finestrata + opaca per la direzione 
                Q_il_str_A_iw += frame_factor*F_sh*surface.glazedArea*(shgc*BRV*shading+shgc_diffuse*(TRV-BRV))*((self.Araum - self.Aaw)/(self.Araum -surface.glazedArea))
                Q_il_str_A_aw += frame_factor*F_sh*surface.glazedArea*(shgc*BRV*shading+shgc_diffuse*(TRV-BRV))*((self.Aaw - surface.glazedArea)/(self.Araum -surface.glazedArea))
                   
            if surface.type == 'GroundFloor':
                theta_eq[:,i] = T_ext*self.Strat[surface.type].U*surface.opaqueArea/self.UA_tot
        
        #self.Q_il_str_A = self.Q_il_str_A.to_numpy()
        #self.carichi_sol = (Q_il_str_A_iw+ Q_il_str_A_aw).to_numpy()
        
        self.theta_eq_tot = theta_eq.sum(axis = 1) + theta_eq_w.sum(axis = 1) 
        
        # Calculates internal heat gains
        
        Q_il_str_I = (self.schedules.people*self.sens_frac*(1-self.conv_people)+self.schedules.appliances*(1-self.conv_app)+self.schedules.lighting*(1-self.conv_app))*self.zone_area
        self.Q_il_kon_I = (self.schedules.people*self.sens_frac*(self.conv_people)+self.schedules.appliances*(self.conv_app)+self.schedules.lighting*(self.conv_app))*self.zone_area
        
        Q_il_str_I_iw = Q_il_str_I * (self.Araum-self.Aaw)/self.Araum
        Q_il_str_I_aw = Q_il_str_I * self.Aaw/self.Araum
        
        self.Q_il_str_iw =  Q_il_str_A_iw + Q_il_str_I_iw
        self.Q_il_str_aw =  Q_il_str_A_aw + Q_il_str_I_aw
        
        # sigma_fhk_iw = 1 - self.sigma_fhk_aw
                
        self.sigma = loadHK(self.sigma_fhk, self.sigma_fhk_aw, self.sigma_hk_str, self.Aaw, self.Araum)                                                          
    
    def Sensible2C(self, flag, phi, H_ve, theta_bound, theta_sup, tau, theta_set = 20., Q_hk = 0.):
        """
        function x = buildingLS_2C_Tset(phi, theta_bound, params, anteil, theta_set, theta_m0, tau)
        buildingLS_2C_Tset solves the linear system (Y*x = q) of VDI6007 at each
        iteration with a given setpoint temperature
        (*)Note: if  phi_HC > 0 HEATING LOAD; phi_HC < 0 COOLING LOAD.
        
        Parameters
            ----------
            flag : string
                string 'Tset' or 'phiset' 
            phi : list of floats
                internal and solar gains, three components (convective, aw and iw)[W]
            Hve : list of floats
                ventilation coefficients (ventilation and infiltration) [W/K]
            theta_bound : list of floats
                equivalent external air temperature and real outdoor air temperature [°C]
            tau: int
                time constant (3600 s for hourly sim)
            theta_set : float
                set-point of considered thermal zone [°C]
            Q_hk : float
                Thermal power entering in the system [W]
            
        Returns
        -------
        np.array
            temperature nodes of the RC model: 
            theta_m_aw thermal mass of AW building components [°C]
            theta_s_aw surface of AW building components [°C]
            theta_lu_star No physical meaning (node obtained from the delta-->star transformation) [°C]
            theta_I_lu internal air temperature [°C]  <--- WHAT WE WILL EXTRACT AS OUTPUT outside the function     
            Q_hk_ges heating/cooling load for maintaining the given setpoint temperature [W]  <--- WHAT WE WILL EXTRACT AS OUTPUT outside the function     
            theta_s_iw surface of IW building components [°C]
            theta_m_iw thermal mass of IW building components [°C]             
        """
        
        # Check input data type
        
        if flag != 'Tset' and flag != 'phiset':
            raise TypeError(f"Ops... thermal zone Sensible2C, bd {self.bd_name}, flag input is not a 'Tset' or 'phiset': flag {flag}")
        if not isinstance(phi, list) or not isinstance(phi[0],float) or not isinstance(phi[1],float) or not isinstance(phi[2],float):
            try: 
                phi[0]=float(phi[0])
                phi[2]=float(phi[2])
                phi[1]=float(phi[1])
            except ValueError:
                raise TypeError(f"Ops... thermal zone Sensible2C, bd {self.bd_name}, phi input is not a list of floats: phi {phi}")
        if not isinstance(H_ve, list) or not isinstance(H_ve[0],float) or not isinstance(H_ve[1],float):
            try: 
                H_ve[0]=float(H_ve[0])
                H_ve[1]=float(H_ve[1])
            except ValueError:
                raise TypeError(f"Ops... thermal zone Sensible2C, bd {self.bd_name}, H_ve input is not a list of floats: H_ve {H_ve}")
        if not isinstance(theta_bound, list) or not isinstance(theta_bound[0],float) or not isinstance(theta_bound[1],float):
            try: 
                theta_bound[0]=float(theta_bound[0])
                theta_bound[1]=float(theta_bound[1])
            except ValueError:
                raise TypeError(f"Ops... thermal zone Sensible2C, bd {self.bd_name}, theta_bound input is not a list of floats: theta_bound {theta_bound}")
        if not isinstance(theta_sup, float):
            try:
                theta_sup = float(theta_sup)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone Sensible2C, bd {self.bd_name}, theta_sup input is not a float: theta_sup {theta_sup}')
        if not isinstance(tau, int):
            try:
                tau = int(tau)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone Sensible2C, bd {self.bd_name}, tau input is not a float: tau {tau}')
        if not isinstance(theta_set, float):
            try:
                theta_set = float(theta_set)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone Sensible2C, bd {self.bd_name}, theta_set input is not a float: theta_set {theta_set}')
        if not isinstance(Q_hk, float):
            try:
                Q_hk = float(Q_hk)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone Sensible2C, bd {self.bd_name}, Q_hk input is not a float: Q_hk {Q_hk}')
                
        # Check input data quality
        
        if H_ve[0] < 0. or H_ve[1] < 0.:
            wrn(f"\n\n Thermal zone Sensible1C, bd {self.bd_name}, the H_ve is negative: H_ve {H_ve}.")
        if not -50. <= theta_bound[0] < 60.:
            wrn(f"\n\n Thermal zone Sensible1C, bd {self.bd_name}, the theta_bound[0] is outside limits [-50,50] °C: theta_bound[0] {theta_bound[0]}.")
        if not -50. <= theta_bound[1] < 60.:
            wrn(f"\n\n Thermal zone Sensible1C, bd {self.bd_name}, the theta_bound[1] is outside limits [-50,50] °C: theta_bound[1] {theta_bound[1]}.")
        if not -50. <= theta_sup < 50.:
            wrn(f"\n\n Thermal zone Sensible1C, bd {self.bd_name}, the theta_sup is outside limits [-50,50] °C: theta_sup {theta_sup}.")
        if not -50. <= theta_set < 50.:
            wrn(f"\n\n Thermal zone Sensible1C, bd {self.bd_name}, the theta_set is outside limits [-50,50] °C: theta_set {theta_set}.")       
        
        # % INPUT
        
        # % Resistances and capacitances of the 7R2C model
        R_lue_ve = 1e20 if H_ve[0] == 0 else 1/H_ve[0]
        R_lue_inf = 1e20 if H_ve[1] == 0 else  1/H_ve[1]
        
        Q_il_kon = phi[0]      # convective heat gains (internal) 
        Q_il_str_aw = phi[1]   # radiant heat gains on surface node AW (internal + solar) 
        Q_il_str_iw = phi[2]   # radiant heat gains on surface node IW (internal + solar) 
        
        theta_A_eq = theta_bound[0]   # equivalent outdoor temperature (sol-air temperature)
        theta_lue = theta_bound[1]     # outdoor air temperature 
       
        if flag == 'Tset':
            theta_I_lu = theta_set          # internal air setpoint
                        
            # MATRIX OF THERMAL TRANSMITTANCES
            
            Y = np.zeros([6,6])
            
            Y[0,0] = -1/self.RrestAW - 1/self.R1AW - self.C1AW/tau
            Y[0,1] = 1/self.R1AW
            
            Y[1,0] = 1/self.R1AW
            Y[1,1] = -1/self.R1AW - 1/ self.RalphaStarAW
            Y[1,2] = 1/ self.RalphaStarAW
            Y[1,3] = self.sigma[1]
            
            Y[2,1] = 1/ self.RalphaStarAW
            Y[2,2] = -1/ self.RalphaStarAW - 1/self.RalphaStarIL -1/self.RalphaStarIW
            Y[2,4] = 1/self.RalphaStarIW
            
            Y[3,2] = 1/self.RalphaStarIL
            Y[3,3] = self.sigma[2]
            
            Y[4,2] = 1/self.RalphaStarIW
            Y[4,3] = self.sigma[0]
            Y[4,4] = -1/self.RalphaStarIW - 1/self.R1IW
            Y[4,5] = 1/self.R1IW
        
            Y[5,4] = 1/self.R1IW
            Y[5,5] = -1/self.R1IW - self.C1IW/tau
            
            # VECTOR OF KNOWN VALUES
            
            q = np.zeros([6,1])
            
            q[0] = -theta_A_eq/self.RrestAW - self.C1AW*self.theta_m0_vdi[0]/tau
            q[1] = -Q_il_str_aw
            q[2] = -theta_I_lu/self.RalphaStarIL
            q[3] = theta_I_lu/self.RalphaStarIL - Q_il_kon \
                    - (theta_lue - theta_I_lu)/R_lue_inf \
                    - (theta_sup - theta_I_lu)/R_lue_ve \
                    + self.Ca*(theta_I_lu - self.Ta0)/tau
            q[4] = -Q_il_str_iw
            q[5] = -self.C1IW*self.theta_m0_vdi[1]/tau
            
            
            # OUTPUT (UNKNOWN) VARIABLES OF THE LINEAR SYSTEM
            
            y = np.linalg.inv(Y).dot(q)
            return np.insert(y,3,theta_set)
            
        elif flag == 'phiset':
            # Note: the heat load in input is already distributed on the 3 nodes
            
            Q_hk_iw = Q_hk*self.sigma[0]        #% radiant heat flow from HVAC system (on surface node IW)
            Q_hk_aw = Q_hk*self.sigma[1]        #% radiant heat flow from HVAC system (on surface node AW)
            Q_hk_kon = Q_hk*self.sigma[2]       #% convective heat flow from HVAC system (on air node)
            
            # MATRIX OF THERMAL TRANSMITTANCES

            Y = np.zeros([6,6])  
            
            Y[0,0] = -1/self.RrestAW - 1/self.R1AW - self.C1AW/tau 
            Y[0,1] = 1/self.R1AW 
            
            Y[1,0] = 1/self.R1AW 
            Y[1,1] = -1/self.R1AW - 1/ self.RalphaStarAW 
            Y[1,2] = 1/ self.RalphaStarAW 
            					  
            
            Y[2,1] = 1/ self.RalphaStarAW 
            Y[2,2] = -1/ self.RalphaStarAW - 1/self.RalphaStarIL -1/self.RalphaStarIW 
            Y[2,3] = 1/self.RalphaStarIL 
            Y[2,4] = 1/self.RalphaStarIW 
            
            Y[3,2] = 1/self.RalphaStarIL 
            Y[3,3] = -1/self.RalphaStarIL  -1/R_lue_inf  -1/R_lue_ve  -self.Ca/tau 
            
            Y[4,2] = 1/self.RalphaStarIW 
            					 
            Y[4,4] = -1/self.RalphaStarIW - 1/self.R1IW 
            Y[4,5] = 1/self.R1IW 
            
            Y[5,4] = 1/self.R1IW 
            Y[5,5] = -1/self.R1IW - self.C1IW/tau 
            
            # VECTOR OF KNOWN TERMS
            
            q = np.zeros(6) 
            q[0] = -theta_A_eq/self.RrestAW - self.C1AW*self.theta_m0_vdi[0]/tau 
            q[1] = -Q_hk_aw - Q_il_str_aw 
            q[2] = 0
            q[3] = -Q_hk_kon - Q_il_kon - theta_lue/R_lue_inf - theta_sup/R_lue_ve -self.Ca*self.Ta0/tau
            q[4] = -Q_hk_iw - Q_il_str_iw 
            q[5] = -self.C1IW*self.theta_m0_vdi[1]/tau 
            
            # OUTPUT LINEAR SYSTEM
            
            y = np.linalg.inv(Y).dot(q)
            return np.insert(y,4,Q_hk)
        
        else:
            return 'wrong flag'
    
    
    def solveDD_Heating(self):
        '''
        Solves the Heating Design power
        it does not require input
        Stationary conditions
        
        Parameters
            ----------
            None
                    
        Returns
        -------
        None.        
        ''' 
        
        T_e = -5
        if T_e < 0:
            p_extsat = 610.5*np.exp((21.875*T_e)/(265.5+T_e))
        else:
            p_extsat = 610.5*np.exp((17.269*T_e)/(237.3+T_e))
        RH_e = 0.8
        T_set_inv = 20
        RH_int_set_H = 0.6
        x_ext = 0.622*(RH_e*p_extsat/(self.p_atm-(RH_e*p_extsat)))
        G_v = 0
        G_da_inf = 0.2*self.rho_air*self.V/3600  
        Hve_inf = 0.1*self.cp_air*self.rho_air*self.V/3600
        G_da_vent = 0.5*self.rho_air*self.V/3600
        Hve_vent = 0.5*self.cp_air*self.rho_air*self.V/3600
        p_intsat = 610.5*np.exp((17.269*T_set_inv)/(237.3+T_set_inv))
        x_int_set = 0.622*(RH_int_set_H*p_intsat/(self.p_atm-(RH_int_set_H*p_intsat)))
        
        self.heatFlow_DDH = 1.2*(self.Htr_op + self.Htr_w)*(T_set_inv - T_e) + Hve_inf*(T_set_inv - T_e) + Hve_vent*(T_set_inv - T_e)

        self.latentFlow_DDH = (G_da_inf*(x_ext-x_int_set) + G_da_vent*(x_ext-x_int_set))*(self.r_0+self.cpv*T_set_inv) + G_v*(self.r_0+self.cpv*T_set_inv)
        self.latentFlow_DDH = -1*self.latentFlow_DDH
        
        
    def solveTZ(self,t,T_e,RH_e,p_extsat,P_DD,tau,mode = '1C'):
        '''
        Solves the thermal zone in the time step t
        
        Parameters
            ----------
            t : int
                timestep of the simulation
            T_e : float
                external temperature [°C]
            RH_e : float
                external relative humidity [0-1]
            p_extsat : float
                atmospheric pressure [Pa]
            P_DD: list of two floats
                Heating and cooling design power [W]
            tau: int
                time constant (3600 s for hourly sim)
            mode : string
                '1C' model or '2C' model
                    
        Returns
        -------
        None.        
        ''' 

        # Check input data type
        
        if mode != '1C' and mode != '2C':
            raise TypeError(f"Ops... thermal zone solveTZ, bd {self.bd_name}, mode input is not a '1C' or '2C': mode {mode}")
        if not isinstance(P_DD, list) or not isinstance(P_DD[0],float) or not isinstance(P_DD[1],float):
            try: 
                P_DD[0]=float(P_DD[0])
                P_DD[1]=float(P_DD[1])
            except ValueError:
                raise TypeError(f"Ops... thermal zone solveTZ, bd {self.bd_name}, P_DD input is not a list of floats: P_DD {P_DD}")
        if not isinstance(t, int):
            try:
                t = int(t)
            except ValueError:
                raise TypeError(f'Ops... thermal zone solveTZ, bd {self.bd_name}, t input is not an int: t {t}')
        if not isinstance(T_e, float):
            try:
                T_e = float(T_e)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone solveTZ, bd {self.bd_name}, T_e input is not a float: T_e {T_e}')
        if not isinstance(RH_e, float):
            try:
                RH_e = float(RH_e)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone solveTZ, bd {self.bd_name}, RH_e input is not a float: RH_e {RH_e}')
        if not isinstance(p_extsat, float):
            try:
                p_extsat = float(p_extsat)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone solveTZ, bd {self.bd_name}, p_extsat input is not a float: p_extsat {p_extsat}')
        if not isinstance(tau, int):
            try:
                tau = int(tau)
            except ValueError:  
                raise TypeError(f'Ops... thermal zone solveTZ, bd {self.bd_name}, tau input is not a float: tau {tau}')
        
        # Check input data quality
        
        if not -50. <= T_e < 50.:
            wrn(f"\n\n Thermal zone solveTZ, bd {self.bd_name}, the T_e is outside limits [-50,50] °C: T_e {T_e}.")
        if not .0 <= RH_e <= 1.:
            wrn(f"\n\n Thermal zone solveTZ, bd {self.bd_name}, the RH_e is outside limits [-50,50] °C: RH_e {RH_e}.")
        
        flag_AHU = True
        P_max = P_DD[0]
        P_min = P_DD[1]
        
        # Inizialization 
        
        T_set_inv = self.schedules.heatingTSP[t]
        T_set_est =self.schedules.coolingTSP[t]
        RH_int_set_H = self.schedules.HeatingRHSP[t]
        RH_int_set_C = self.schedules.CoolingRHSP[t]
        AHUOnOff = self.schedules.AHUOnOff[t]
        AHUHUM = self.schedules.AHUHUM[t]
        self.T_sup[t] = self.schedules.AHUTSupp[t]
        x_sup = self.schedules.AHUxSupp[t]
        Sens_Recovery_eff = self.schedules.sensRec[t]
        Lat_Recovery_eff = self.schedules.latRec[t]
        OutAirRatio = self.schedules.outdoorAirRatio[t]
              
        self.x_ext[t] = 0.622*(RH_e*p_extsat/(self.p_atm-(RH_e*p_extsat)))
        self.G_v[t] = self.schedules.vapour[t]*self.zone_area
        
        # Infiltration 
        
        self.G_da_inf[t] = self.schedules.infFlowRate[t]*self.rho_air*self.V/3600  
        Hve_inf = self.schedules.infFlowRate[t]*self.cp_air*self.rho_air*self.V/3600
        
        # Ventilation
        
        self.G_da_vent[t] =self.schedules.ventFlowRate[t]*self.rho_air*self.zone_area #/3600
        Hve_vent = self.schedules.ventFlowRate[t]*self.cp_air*self.rho_air*self.zone_area #/3600
        
        Hve = [Hve_vent , Hve_inf]
        
        # Set internal gains and supply temperature
        
        if mode == '1C':
            phi_load = [self.phi_ia[t], self.phi_st[t],self.phi_m[t]]
        elif mode == '2C':
            phi_load = [self.Q_il_kon_I[t], self.Q_il_str_aw[t], self.Q_il_str_iw[t]]
            theta_bound = [self.theta_eq_tot[t],T_e]
        
        if AHUOnOff  == 0:
            self.T_sup[t] = T_e
            x_sup = self.x_ext[t]
            
        # Stats the zone system solution.. The logic depend on the heating cooling and latent sensible
        
        while flag_AHU:
            
            # SENSIBLE HEAT LOAD CALCULATION 
        
            # Heating mode 
            
            if self.schedules.plantOnOffSens[t] == 1:
                
                if mode == '1C':
                    pot, Ta, Ts, Tm = self.Sensible1C('Tset',Hve,T_set_inv, self.T_sup[t], T_e , phi_load,tau)
                elif mode == '2C':
                    Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.Sensible2C('Tset',phi_load, Hve, theta_bound, self.T_sup[t], tau, theta_set = T_set_inv)
                
                if pot>0 and pot<P_max*1000:
                    self.heatFlow[t] = pot
                    self.Air_temp[t] = Ta
                    if mode == '1C':
                        self.theta_m = Tm
                    elif mode == '2C':
                        self.theta_m_vdi = [Tm_aw,Tm_iw]
                else:
                    if pot > P_max*1000:
                        phi_set = P_max*1000
                    else:
                        phi_set = 0
                    
                    if mode == '1C':
                        pot, Ta, Ts, Tm = self.Sensible1C('phiset',Hve,T_set_inv, self.T_sup[t], T_e , phi_load,tau, phi_HC_set = phi_set)
                    elif mode == '2C':
                        Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.Sensible2C('phiset',phi_load, Hve, theta_bound, self.T_sup[t],  tau, Q_hk = phi_set)
                    self.heatFlow[t] = pot
                    self.Air_temp[t] = Ta
                    if mode == '1C':
                        self.theta_m = Tm
                    elif mode == '2C':
                        self.theta_m_vdi = [Tm_aw,Tm_iw]
            
            # Cooling mode
            
            if self.schedules.plantOnOffSens[t] == -1:
                
                if mode == '1C':
                    pot, Ta, Ts, Tm = self.Sensible1C('Tset',Hve,T_set_est, self.T_sup[t], T_e , phi_load,tau)
                elif mode == '2C':
                    Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.Sensible2C('Tset',phi_load, Hve, theta_bound,self.T_sup[t],  tau, theta_set = T_set_est)
                
                if pot<0 and pot>P_min*1000:
                    self.heatFlow[t] = pot
                    self.Air_temp[t] = Ta
                    if mode == '1C':
                        self.theta_m = Tm
                    elif mode == '2C':
                        self.theta_m_vdi = [Tm_aw,Tm_iw]    
                else:
                    if pot < P_min*1000:
                        phi_set = P_min*1000
                    else:
                        phi_set = 0
                    if mode == '1C':
                        pot, Ta, Ts, Tm = self.Sensible1C('phiset',Hve,T_set_est, self.T_sup[t], T_e , phi_load,tau, phi_HC_set = phi_set)
                    elif mode == '2C':
                        Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.Sensible2C('phiset',phi_load, Hve, theta_bound,self.T_sup[t],  tau, Q_hk = phi_set)
                    self.heatFlow[t] = pot
                    self.Air_temp[t] = Ta
                    if mode == '1C':
                        self.theta_m = Tm
                    elif mode == '2C':
                        self.theta_m_vdi = [Tm_aw,Tm_iw]
            
            # Plant OFF
            
            if self.schedules.plantOnOffSens[t] == 0:
                
                if mode == '1C':
                    pot, Ta, Ts, Tm = self.Sensible1C('phiset',Hve,T_set_inv, self.T_sup[t], T_e , phi_load,tau)
                elif mode == '2C':
                    Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.Sensible2C('phiset',phi_load, Hve, theta_bound,self.T_sup[t],tau)
                self.heatFlow[t] = pot
                self.Air_temp[t] = Ta
                if mode == '1C':
                    self.theta_m = Tm
                elif mode == '2C':
                    self.theta_m_vdi = [Tm_aw,Tm_iw]          
            
            
            # LATENT HEAT LOAD CALCULATION
                        
            # Dehumidification mode
            
            if self.schedules.plantOnOffLat[t] == -1:
                
                if self.Air_temp[t] < 0:
                    p_intsat = 610.5*np.exp((21.875*self.Air_temp[t])/(265.5+self.Air_temp[t]))
                else:
                    p_intsat = 610.5*np.exp((17.269*self.Air_temp[t])/(237.3+self.Air_temp[t]))
                
                x_int_set = 0.622*(RH_int_set_C*p_intsat/(self.p_atm-(RH_int_set_C*p_intsat)))
                self.latentFlow[t] = (self.G_da_inf[t]*(self.x_ext[t]-x_int_set) + self.G_da_vent[t]*(x_sup-x_int_set)-self.rho_air*self.V*(x_int_set-self.x_m0)/tau)*(self.r_0+self.cpv*self.Air_temp[t]) + self.G_v[t]*(self.r_0+self.cpv*self.Air_temp[t])
                
                if self.latentFlow[t] < 0:
                    self.latentFlow[t] = 0
                    self.x_int[t] =(self.G_da_inf[t]*self.x_ext[t] + self.G_da_vent[t]*x_sup + self.G_v[t] + self.rho_air*self.V*self.x_m0/tau)/(self.G_da_inf[t] + self.G_da_vent[t] + self.rho_air*self.V/tau)
                    p_int = self.p_atm*self.x_int[t]/(0.622+self.x_int[t])
                    RH_int = p_int/p_intsat
                    self.RH_i[t] = RH_int*100
                   
                else:
                    self.RH_i[t] = RH_int_set_C*100
                    self.x_int[t] = x_int_set
            
            # Humidification mode
            
            if self.schedules.plantOnOffLat[t] == 1:
                
                if self.Air_temp[t] < 0:
                    p_intsat = 610.5*np.exp((21.875*self.Air_temp[t])/(265.5+self.Air_temp[t]))
                else:
                    p_intsat = 610.5*np.exp((17.269*self.Air_temp[t])/(237.3+self.Air_temp[t]))
                
                x_int_set = 0.622*(RH_int_set_H*p_intsat/(self.p_atm-(RH_int_set_H*p_intsat)))
                self.latentFlow[t] = (self.G_da_inf[t]*(self.x_ext[t]-x_int_set) + self.G_da_vent[t]*(x_sup-x_int_set)-self.rho_air*self.V*(x_int_set-self.x_m0)/tau)*(self.r_0+self.cpv*self.Air_temp[t]) + self.G_v[t]*(self.r_0+self.cpv*self.Air_temp[t])
                
                if self.latentFlow[t] > 0:
                    self.latentFlow[t] = 0
                    self.x_int[t] =(self.G_da_inf[t]*self.x_ext[t] + self.G_da_vent[t]*x_sup + self.G_v[t] + self.rho_air*self.V*self.x_m0/tau)/(self.G_da_inf[t] + self.G_da_vent[t] + self.rho_air*self.V/tau)
                    p_int = self.p_atm*self.x_int[t]/(0.622+self.x_int[t])
                    RH_int = p_int/p_intsat
                    self.RH_i[t] = RH_int*100
                    
                    if RH_int > 1:
                        RH_int = 0.99
                        self.RH_i[t] = RH_int*100
                        p_int = p_intsat*RH_int
                        self.x_int[t] = 0.622*(p_int/(self.p_atm-p_int))
                    
                else:
                    self.RH_i[t] = RH_int_set_H*100
                    self.x_int[t] = x_int_set
            
            # Plant Off
            
            if self.schedules.plantOnOffLat[t] == 0:
                if self.Air_temp[t] < 0:
                    p_intsat = 610.5*np.exp((21.875*self.Air_temp[t])/(265.5+self.Air_temp[t]))
                else:
                    p_intsat = 610.5*np.exp((17.269*self.Air_temp[t])/(237.3+self.Air_temp[t]))
                    self.x_int[t] = (self.G_da_inf[t]*self.x_ext[t] + self.G_da_vent[t]*x_sup + self.G_v[t] + self.rho_air*self.V*self.x_m0/tau)/(self.G_da_inf[t] + self.G_da_vent[t] + self.rho_air*self.V/tau)
                    p_int = self.p_atm*self.x_int[t]/(0.622+self.x_int[t])
                    RH_int = p_int/p_intsat
                    self.RH_i[t] = RH_int*100
                    
                    if RH_int > 1:
                        RH_int = 0.99
                        self.RH_i[t] = RH_int*100
                        p_int = p_intsat*RH_int
                        self.x_int[t] = 0.622*(p_int/(self.p_atm-p_int))
                
                self.latentFlow[t] = 0
            
            self.latentFlow[t] = -1*self.latentFlow[t]
            
            # AIR HANDLING UNIT HEAT DEMAND 
            
            self.ZoneAHU.AHUCalc(t,self.G_da_vent[t],AHUOnOff,AHUHUM,Sens_Recovery_eff,Lat_Recovery_eff,OutAirRatio,T_e,self.x_ext[t],self.Air_temp[t],self.x_int[t],self.T_sup[t],x_sup)
        
            # Check if supply conditions are changed
            
            err_T_sup = abs(self.T_sup[t] - self.ZoneAHU.T_supAHU[t])
            err_x_sup = abs(x_sup - self.ZoneAHU.x_sup)
            if err_T_sup > 0.1 or err_x_sup > 0.0001:
                self.T_sup[t] = self.ZoneAHU.T_supAHU[t]
                x_sup = self.ZoneAHU.x_sup
            else:
                if mode == '1C':
                    self.theta_m0 = self.theta_m
                    self.T_wall_0 = T_e + (Tm - T_e)*self.Htr_em/self.H_lim_ext_wall_1C
                elif mode == '2C':
                    self.theta_m0_2c[t] = self.theta_m_vdi
                    self.ext_surf_heat_gain[t] = self.C1AW/tau*(self.theta_m_vdi[0]-self.theta_m0_vdi[0])
                    self.theta_m0_vdi = self.theta_m_vdi
                    self.T_wall_0 = self.theta_eq_tot[t] + (Tm_aw - self.theta_eq_tot[t])*self.R_lim_ext_wall_2C/self.RrestAW
                self.Ta0 = Ta
                self.x_m0 = self.x_int[t]
                self.T_wall_0_vector[t] = self.T_wall_0
                flag_AHU = False

                
    def reset_init_values(self):
        '''
        This method allows to reset temperatures starting values after
        the DesignDays calculation and before running simulation
        
        Parameters
            ----------
            None
                    
        Returns
        -------
        None.        
        ''' 
            
        self.Ta0 = 15
        self.theta_m0 = 15
        self.theta_m0_vdi = [15, 15]                                           # aw the first, iw the second
        
        
#%% ---------------------------------------------------------------------------------------------------
#%% Building class     

class Building:
    
    '''
    building class
    generates the single thermal zone or multifloor thermal zone
    sets the internal heat gains 
    
    Methods:
        init
        BDParamsandLoads
        BDdesigndays_Heating
        BDdesigndays_Cooling
        BDplants
        solve
        checkExtWalls
        checkExtWallsMod
        checkExtWallsMod2
        plotBuilding
        printBuildingInfo
    
    '''
    
    # Class variables
    
    Rat_iw = 2.5    # Internal walls area multiplicator [-]. Multiplies the total floor area to calc internal wall area (ISO13790)
    eps_se = 0.9
    h_r = 5*eps_se  # external radiative heat transfer coefficient [W/(m2 K)]
    
    def __init__(self,buildingName,mode,surfList,n_Floors,
                 end_use,envelopes,archId,
                 rh_net,rh_gross,
                 heating_plant,cooling_plant,weather):
        '''
        Initializes the building object
        Builds up teh building geometry with respet to the input data
        creates some np.array which will contain the zone parameters and variables
        
        Parameters
            ----------
            buildingName : string
                string with the building name or id
            mode : str
                This is or 'cityjson' or 'geojson'. Depending on this the number of floors is calculated differently
            surfList : list of lists
                this list includes all surfaces vertices. three layer list [m]
                
                list of surfaces [ surf1 
                                        [v1[x1 y1 z1],  v2[x2 y2 z2],   v3[x3 y3 z3], ..........],
                                   surf2
                                        [v1[x1 y1 z1],  v2[x2 y2 z2],   v3[x3 y3 z3], ..........],
                                   surf3
                                        [v1[x1 y1 z1],  v2[x2 y2 z2],   v3[x3 y3 z3], ..........],
                                   .
                                   .
                                   .
                                   .
                                            
                ]
            n_Floors : number of floors
                used only for geojson
            end_use : str
                the string with the end use archetype name
            envelopes : dictionary
                dictionary of Envelopes objects that returns from loadEnvelopes
            archId : string
                name of the envelope archetype. It is used as key of the envelopes dictionary
            rh_net : float
                multiplicator of the building's volume
            rh_gross : float
                multiplicator of the surface area
            heating_plant : str
                string with the heating plant typology
            cooling_plant : str
                string with the cooling plant typology


        Returns
        -------
        None.        
        '''         
        
        # Check input data type
        
        if not isinstance(buildingName, str):
            try:
                buildingName = str(buildingName)
            except ValueError:
                raise TypeError(f'Ops... Building class, input buildingName is not a string: buildingName {buildingName}')        
        if not isinstance(surfList, list) and not isinstance(surfList, tuple) :
            raise TypeError(f'Ops... Building class, bd {buildingName}, input surfList is not a list: surfList {surfList}')        
        if not isinstance(end_use, str):
            raise TypeError(f'Ops... Building class, bd {buildingName}, input end_use is not a string: end_use {end_use}')    
        if not isinstance(mode, str):
            raise TypeError(f'Ops... Building class, bd {buildingName}, input mode is not a string: mode {mode}')                   
        if not isinstance(n_Floors, int):
            try:
                n_Floors = int(n_Floors)
            except:
                raise TypeError(f'Ops... Building class, bd {buildingName}, n_Floors is not a int: n_Floors {n_Floors}') 
        if not isinstance(envelopes, dict):
            raise TypeError(f'Ops... Building class, bd {buildingName}, envelopes is not a dictionary: envelopes {envelopes}... print help for more info') 
        if not isinstance(archId, str):
            raise TypeError(f'Ops... Building class, bd {buildingName}, input archId is not a string: archId {archId}')  
        if not isinstance(rh_gross, float):
            try:
                rh_gross = float(rh_gross)
            except ValueError:
                raise TypeError(f'Ops... Building class, bd {buildingName}, rh_gross is not a float: rh_gross {rh_gross}') 
        if not isinstance(rh_net, float):
            try:
                rh_net = float(rh_net)
            except ValueError:
                raise TypeError(f'Ops... Building class, bd {buildingName}, rh_net is not a float: rh_net {rh_net}') 
        if not isinstance(heating_plant, str):
            raise TypeError(f'Ops... Building class, bd {buildingName}, input heating_plant is not a string: heating_plant {heating_plant}')  
        if not isinstance(cooling_plant, str):
            raise TypeError(f'Ops... Building class, bd {buildingName}, input cooling_plant is not a string: cooling_plant {cooling_plant}')  
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')

        # check input data quality
        
        if mode != 'cityjson' and mode != 'geojson':
            raise ValueError(f'Ops... Building class, bd {buildingName}, input mode must be "cityjson" or "geojson": mode {mode}')   
        if rh_gross < 0.5 or rh_gross > 1.5:
            wrn(f"\n\nBuilding class, init, bd {buildingName}, are you sure about the external walls multiplication coeff?? rh_gross {rh_gross}\n")
        if rh_net < 0.5 or rh_net > 1.5:
            wrn(f"\n\nBuilding class, init, bd {buildingName}, are you sure about the external walls multiplication coeff?? rn_net {rn_net}\n")
        
        '''
        one thermal zone for building
        ---
        THIS IS FOR LOD1
        ---
        '''
        ts = weather.sim_time[0]
        hours = weather.sim_time[1]
        
        # Output vectors inizialization 
        
        self.l = ts*(hours-1) + 1
        self.Air_tempBD = np.zeros(self.l)
        self.RH_iBD = np.zeros(self.l)
        self.heatFlowBD = np.zeros(self.l)
        self.latentFlowBD = np.zeros(self.l)
        self.AHUDemandBD = np.zeros(self.l)
        self.AHUDemand_latBD = np.zeros(self.l)
        self.AHUDemand_sensBD = np.zeros(self.l)
        
        # set some building's attributes
        self.__flagCheckWalls = False
        self.name = buildingName
        self.archId = archId
        self.end_use = end_use
        self.footprint = 0.0000001
        self.oneFloorHeight = 3.3
        self.zones={}
        self.buildingSurfaces = {}

        i = 0
        
        self.hmax = []
        self.hmin = []
        self.extWallArea = 0
        self.extRoofArea = 0
        self.extWallWinArea = 0
        self.extWallOpaqueArea = 0
        self.T_out_inf = 15.
        self.T_wall_0 = 15.
        self.G_inf_0 = 0.
        self.T_out_AHU = 15.
        self.G_vent_0 = 0.
        self.H_waste = 0.
        
        # Calculation of wwr and surfaces init
        
        try:
            for surf in surfList:
                self.wwr = envelopes[archId].Window.wwr
                i += 1
                surface = Surface('Building Surface '+str(i),weather.azSubdiv,weather.hSubdiv,self.wwr,rh_gross,surf)
                self.buildingSurfaces['Building Surface '+str(i)]=surface
                self.ii = i
        except KeyError:
            raise KeyError(f"\n\n Building class: bd {self.name}, the envelope archetype not found. key {archId}\nList of possible archetypes: \n {envelopes.keys()}")
        
        # check if some surface is coincident
        
        self.checkExtWallsMod2()    
        
        # Calculation of external wall area and other geometrical params
        
        for surface in self.buildingSurfaces.values():
            if surface.type == 'GroundFloor':
                self.footprint += surface.area
            if surface.type == 'ExtWall':
                self.extWallWinArea += surface.glazedArea
                self.extWallOpaqueArea += surface.opaqueArea
                self.hmax.append(surface.maxHeight())
                self.hmin.append(surface.minHeight())
            if surface.type == 'Roof':
                self.extRoofArea += surface.area            
        
        self.hmax = np.mean(self.hmax)
        self.hmin = np.mean(self.hmin)
        self.buildingHeight = self.hmax - self.hmin
        self.Vertsurf = []
        for surf in self.buildingSurfaces.values():
            if surf.type == 'ExtWall':
                self.extWallArea += surf.area
                surf.surfHeight = self.buildingHeight
                self.Vertsurf.append([surf,[]])
        
        # Number of floor calc
        
        if mode == 'geojson':
            self.nFloors = n_Floors
        else:
            self.nFloors = round(self.buildingHeight/self.oneFloorHeight)
        
        # Footprint and volume calc
        
        if self.footprint == 0:
            self.footprint = 1.
        self.total_area = (self.nFloors)*self.footprint
        self.Volume = rh_net*(self.oneFloorHeight*self.footprint*self.nFloors)
        if self.Volume == 0:
            self.Volume = 0.0001
        
        # Set plants data
        
        self.H_plant_type = heating_plant
        self.C_plant_type = cooling_plant
        self.BDPlant = Plants(self.H_plant_type,self.C_plant_type,self.l,ts)
        self.Pnom_H_BD = 1e20
        self.Pnom_C_BD = -1e20

        
    def BDParamsandLoads(self,model,envelopes,sched_db,weather,splitInZone=False):
        '''
        This method allows to initialize thermal zones and to calculate
        Parameters and Internal and Solar loads for each zone of the buildings
        
        Parameters
            ----------
            model  : string
                '1C' model or '2C' model
            envelopes : dictionary
                dictionary of Envelopes objects that returns from loadEnvelopes 
            sched_db : dictionary
                dictionary of Archetypes objects that returns from loadArchetype 
            weather : RC_classes.WeatherData.Weather obj
                object of the class weather WeatherData module
            splitInZone : boolean
                Only True available

        Returns
        -------
        None.        
        '''   
 
        # Check input data type
 
        if not isinstance(model, str):
            raise TypeError(f'Ops... Building class, BDParamsandLoads, bd {self.name}, input model is not a string: model {model}')         
        if not isinstance(envelopes, dict):
            raise TypeError(f'Ops... Building class, BDParamsandLoads, bd {self.name}, envelopes is not a dictionary: envelopes {envelopes}... print help for more info') 
        if not isinstance(sched_db, dict):
            raise TypeError(f'Ops... Building class, BDParamsandLoads, bd {self.name}, sched_db is not a dictionary: sched_db {sched_db}... print help for more info') 
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
        if not isinstance(splitInZone, bool):
            raise TypeError(f'Ops... Building class, BDParamsandLoads, bd {self.name}, splitInZone is not a boolean: splitInZone {splitInZone}')
        
        # Check input data quality
        
        if splitInZone:
            raise TypeError(f'Ops... Building class, BDParamsandLoads, bd {self.name}, splitInZone must be False: splitInZone {splitInZone}. MultiZone building not implemented yet')
        
        # Creates thermal zone and calculates the params
        
        if not splitInZone:
            self.buildingSurfaces['Building Surface '+str(self.ii+1)] = SurfaceInternalMass(('Building Surface '+str(self.ii+1)),(self.total_area-self.footprint),surfType='IntCeiling')
            self.buildingSurfaces['Building Surface '+str(self.ii+2)] = SurfaceInternalMass(('Building Surface '+str(self.ii+1)),(self.total_area-self.footprint),surfType='IntFloor')
            self.buildingSurfaces['Building Surface '+str(self.ii+3)] = SurfaceInternalMass(('Building Surface '+str(self.ii+2)),self.Rat_iw*self.nFloors*self.footprint,surfType='IntWall')
            try:
                self.zones['Zone'] = ThermalZone(1,self.name,
                                                 envelopes[self.archId],sched_db[self.end_use],
                                                 self.buildingSurfaces.values(),
                                                 self.Volume,
                                                 self.footprint*self.nFloors,
                                                 self.l)
            except KeyError:
                raise KeyError(f"\n\n Building class, BDParamsandLoads: bd {self.name}, the end use archetype not found. key {self.end_use}\nList of possible archetypes: \n {sched_db.keys()}")
        
            if model == '1C':
                self.zones['Zone'].zoneParameter13790()
                self.zones['Zone'].calculate_zone_loads_ISO13790(weather,self.h_r)
            
            elif model == '2C':
                self.zones['Zone'].zoneParameter6007()
                self.zones['Zone'].calculate_zone_loads_vdi6007(weather)
            
            else:
                raise ValueError(f'Building class, BDParamsandLoads, bd {self.name}. Give a proper model: "2C" or "1C"')
        
        # PV calc
        
        self.BuildPV = DistrictPVGIS(self.extRoofArea,self.l)
        self.BuildPV.PV_calc(weather)
        
        
    def BDdesigndays_Heating(self,Plant_calc):
        '''
        This method allows to calculate Heating maximum power required by
        buildings during design days
        
        Parameters
            ----------
            Plant_calc  : string
                'YES' or 'NO' 
        
        Returns
        -------
        None.             
        '''
        
        # Check input data type
        
        if not isinstance(Plant_calc, str):
            raise TypeError(f'Ops... Building class, BDdesigndays_Heating, bd {self.name}, input Plant_calc is not a string: Plant_calc {Plant_calc}')         
        
        # Check input data quality
        
        if Plant_calc!='YES' and Plant_calc!='NO':
            raise ValueError(f"Building class, BDdesigndays_Heating, bd {self.name}. Give a proper Plant_calc: 'YES' or 'NO' ")
        
        for z in self.zones.values():
            z.solveDD_Heating()
            self.Pnom_H_BD = (self.zones['Zone'].heatFlow_DDH + self.zones['Zone'].latentFlow_DDH)/1000  # [kW]
            self.PDesignH = self.Pnom_H_BD
            
        if Plant_calc == 'NO':
            self.Pnom_H_BD = 1e20
        
    
    def BDdesigndays_Cooling(self,t,T_e,RH_e,p_extsat,tau,Plant_calc,model = '1C'):

        '''
        This method allows to calculate Cooling maximum power required by
        buildings during design days
        
        Parameters
            ----------
            t : int
                timestep of the simulation
            T_e : float
                external temperature [°C]
            RH_e : float
                external relative humidity [0-1]
            p_extsat : float
                atmospheric pressure [Pa]
            tau: int
                time constant (3600 s for hourly sim)
            Plant_calc  : string
                'YES' or 'NO' 
            model : string
                '1C' model or '2C' model
                    
        Returns
        -------
        None.        
        ''' 

        # Check input data type
        
        if model != '1C' and model != '2C':
            raise TypeError(f"Ops... Building class, BDdesigndays_Cooling, bd {self.name}, model input is not a '1C' or '2C': model {model}")
        if not isinstance(t, int):
            try:
                t = int(t)
            except ValueError:
                raise TypeError(f'Ops... Building class, BDdesigndays_Cooling, bd {self.name}, t input is not an int: t {t}')
        if not isinstance(T_e, float):
            try:
                T_e = float(T_e)
            except ValueError:  
                raise TypeError(f'Ops... Building class, BDdesigndays_Cooling, bd {self.name}, T_e input is not a float: T_e {T_e}')
        if not isinstance(RH_e, float):
            try:
                RH_e = float(RH_e)
            except ValueError:  
                raise TypeError(f'Ops... Building class, BDdesigndays_Cooling, bd {self.name}, RH_e input is not a float: RH_e {RH_e}')
        if not isinstance(p_extsat, float):
            try:
                p_extsat = float(p_extsat)
            except ValueError:  
                raise TypeError(f'Ops... Building class, BDdesigndays_Cooling, bd {self.name}, p_extsat input is not a float: p_extsat {p_extsat}')
        if not isinstance(tau, int):
            try:
                tau = int(tau)
            except ValueError:  
                raise TypeError(f'Ops... Building class, BDdesigndays_Cooling, bd {self.name}, tau input is not a float: tau {tau}')
        if not isinstance(Plant_calc, str):
            raise TypeError(f'Ops... Building class, BDdesigndays_Cooling, bd {self.name}, input Plant_calc is not a string: Plant_calc {Plant_calc}')         
        
        # Check input data quality
        
        if Plant_calc!='YES' and Plant_calc!='NO':
            raise ValueError(f"Building class, BDdesigndays_Cooling, bd {self.name}. Give a proper Plant_calc: 'YES' or 'NO' ")
        if not -50. <= T_e < 50.:
            wrn(f"\n\n Building class, BDdesigndays_Cooling, bd {self.name}, the T_e is outside limits [-50,50] °C: T_e {T_e}.")
        if not .0 <= RH_e <= 1.:
            wrn(f"\n\n Building class, BDdesigndays_Cooling, bd {self.name}, RH_e is outside limits [-50,50] °C: RH_e {RH_e}.")
        
        # Simulation for the summer design period
        
        for z in self.zones.values():
            z.solveTZ(t,T_e,RH_e,p_extsat,1e20*np.array([1,-1]),tau,model)
            self.Pnom_C_BD = min((self.zones['Zone'].heatFlow + self.zones['Zone'].latentFlow)/1000 + self.zones['Zone'].ZoneAHU.AHUDemand)   # [kW]
            self.PDesignC = self.Pnom_C_BD
            
        if Plant_calc == 'NO':
            self.Pnom_C_BD = -1e20


    def BDplants(self,Plants_list,T_ext_H_avg):
        
        '''
        This method allows to set the plant of each building and to
        check minimum plant efficiency
        
        Parameters
            ----------
            Plants_list: dictionary
                dictionary with plant_key/Plants data
            T_ext_H_avg : float
                external average temperature during the heating seasn
                    
        Returns
        -------
        None.        
        '''
        
        # Check input data type
        
        if not isinstance(Plants_list, dict):
            raise TypeError(f' Building class, BDplants, bd {self.name}, Plants_list must be a dictionary: Plants_list {Plants_list}')
        if not isinstance(T_ext_H_avg, float):
            raise TypeError(f' Building class, BDplants, bd {self.name}, T_ext_H_avg must be a float: T_ext_H_avg {T_ext_H_avg}')
                
        # Set building plants
        
        self.BDPlant.setPlant(Plants_list,self.Pnom_H_BD,self.Pnom_C_BD,T_ext_H_avg)
        
        
    def solve(self,t,Plants_list,T_e,RH_e,p_extsat,tau,Plant_calc,model = '1C'):

        '''
        This method allows to calculate Cooling maximum power required by
        buildings during design days
        
        Parameters
            ----------
            t : int
                timestep of the simulation
            Plants_list: dictionary
                dictionary with plant_key/Plants data
            T_e : float
                external temperature [°C]
            RH_e : float
                external relative humidity [0-1]
            p_extsat : float
                atmospheric pressure [Pa]
            tau: int
                time constant (3600 s for hourly sim)
            Plant_calc  : string
                'YES' or 'NO' 
            model : string
                '1C' model or '2C' model
                    
        Returns
        -------
        None.        
        ''' 

        # Check input data type
        
        if model != '1C' and model != '2C':
            raise TypeError(f"Ops... Building class, solve, bd {self.name}, model input is not a '1C' or '2C': model {model}")
        if not isinstance(Plants_list, dict):
            raise TypeError(f'Building class, solve, bd {self.name}, Plants_list must be a dictionary: Plants_list {Plants_list}')        
        if not isinstance(t, int):
            try:
                t = int(t)
            except ValueError:
                raise TypeError(f'Ops... Building class, solve, bd {self.name}, t input is not an int: t {t}')
        if not isinstance(T_e, float):
            try:
                T_e = float(T_e)
            except ValueError:  
                raise TypeError(f'Ops... Building class, solve, bd {self.name}, T_e input is not a float: T_e {T_e}')
        if not isinstance(RH_e, float):
            try:
                RH_e = float(RH_e)
            except ValueError:  
                raise TypeError(f'Ops... Building class, solve, bd {self.name}, RH_e input is not a float: RH_e {RH_e}')
        if not isinstance(p_extsat, float):
            try:
                p_extsat = float(p_extsat)
            except ValueError:  
                raise TypeError(f'Ops... Building class, solve, bd {self.name}, p_extsat input is not a float: p_extsat {p_extsat}')
        if not isinstance(tau, int):
            try:
                tau = int(tau)
            except ValueError:  
                raise TypeError(f'Ops... Building class, solve, bd {self.name}, tau input is not a float: tau {tau}')
        if not isinstance(Plant_calc, str):
            raise TypeError(f'Ops... Building class, solve, bd {self.name}, input Plant_calc is not a string: Plant_calc {Plant_calc}')         
        
        # Check input data quality
        
        if Plant_calc!='YES' and Plant_calc!='NO':
            raise ValueError(f"Building class, solve, bd {self.name}. Give a proper Plant_calc: 'YES' or 'NO' ")
        if not -50. <= T_e < 50.:
            wrn(f"\n\n Building class, solve, bd {self.name}, the T_e is outside limits [-50,50] °C: T_e {T_e}.")
        if not .0 <= RH_e <= 1.:
            wrn(f"\n\n Building class, solve, bd {self.name}, RH_e is outside limits [-50,50] °C: RH_e {RH_e}.")
        
        # Thermal zone simulation
        
        for z in self.zones.values():
            z.solveTZ(t,T_e,RH_e,p_extsat,[self.Pnom_H_BD,self.Pnom_C_BD],tau,model)
            self.Air_tempBD[t] = self.zones['Zone'].Air_temp[t]
            self.RH_iBD[t] = self.zones['Zone'].RH_i[t]
            self.heatFlowBD[t] = self.zones['Zone'].heatFlow[t]/1000           # [kW]
            self.latentFlowBD[t] = self.zones['Zone'].latentFlow[t]/1000       # [kW]
            self.AHUDemandBD[t] = self.zones['Zone'].ZoneAHU.AHUDemand[t]      # [kW]
            self.AHUDemand_latBD[t] = self.zones['Zone'].ZoneAHU.AHUDemand_lat[t]
            self.AHUDemand_sensBD[t] = self.zones['Zone'].ZoneAHU.AHUDemand_sens[t]
        
        # Calculates some data for urban climate
        
        self.T_out_inf = np.mean(self.Air_tempBD[t])                           # for Urban Canyon
        self.T_wall_0 = np.mean(self.zones['Zone'].T_wall_0)                   # for Urban Canyon
        self.G_inf_0 = self.zones['Zone'].G_da_inf[t] / self.zones['Zone'].rho_air
        self.T_out_AHU = self.zones['Zone'].ZoneAHU.T_out
        self.G_vent_0 = self.zones['Zone'].G_da_vent[t] / self.zones['Zone'].rho_air
        
        if Plant_calc == 'YES':
            self.BDPlant.solvePlant(Plants_list,self.heatFlowBD[t] + self.latentFlowBD[t] + self.AHUDemandBD[t],t,T_e,self.Air_tempBD[t],self.RH_iBD[t])
            self.H_waste = self.BDPlant.H_waste[t]
        elif Plant_calc == 'NO':
            pass
            
    def checkExtWalls(self):
        
        '''
        DEPRECATED
        This method checks if there are coincident external walls
        and, in case, corrects their area
        
        Parameters
            ----------
            None
                    
        Returns
        -------
        None.    
        
        '''
        
        # Lists surface and checks their coincidance        

        if not self.__flagCheckWalls:
            self.__flagCheckWalls = True
            for surf in self.buildingSurfaces.values():
                if surf.type=='GroundFloor':
                    for i in range(len(surf.vertList)):
                        v1=surf.vertList[i-1]
                        v2=surf.vertList[i]
                        coincidentWalls=[]
                        for wall in self.buildingSurfaces.values():
                            if wall.type=='ExtWall':
                                if (v1 in wall.vertList and v2 in wall.vertList):
                                    coincidentWalls.append(wall)
                        if len(coincidentWalls)==2:
                            #n += 1
                            #surf.printInfo()
                            #coincidentWalls[0].printInfo()
                            #coincidentWalls[1].printInfo()
                            if coincidentWalls[0].area >= coincidentWalls[1].area:
                                coincidentWalls[0].area = coincidentWalls[0].area - coincidentWalls[1].area
                                coincidentWalls[1].area = 0
                            else:
                                coincidentWalls[1].area = coincidentWalls[1].area - coincidentWalls[0].area
                                coincidentWalls[0].area = 0
                                
                        
                        if len(coincidentWalls)>2:
                            print('WARNING: building "'+buildingName+'" could have 3 or more coincident walls')
                            [i.printInfo() for i in coincidentWalls]
    
    
    def checkExtWallsMod(self):
        '''
        DEPRECATED
        This method checks if there are coincident external walls
        and, in case, corrects their area
        
        Parameters
            ----------
            None
                    
        Returns
        -------
        None.    
        
        '''
        
        # Lists surface and checks their coincidance  
        
        if not self.__flagCheckWalls:
            self.__flagCheckWalls = True
            for wall1 in self.buildingSurfaces.values():
                if wall1.type=='ExtWall':
                    for wall2 in self.buildingSurfaces.values():
                        if (wall2.type=='ExtWall' and wall2 != wall1 and wall1.checkSurfaceCoincidence(wall2)):
                            
                            if wall1.area >= wall2.area:
                                wall1.area = wall1.area - wall2.area
                                wall2.area = 0
                            else :
                                wall2.area = wall2.area - wall1.area
                                wall1.area = 0 
    
    
    def checkExtWallsMod2(self):
        '''
        This method checks if there are coincident external walls
        and, in case, corrects their area
        
        Parameters
            ----------
            None
                    
        Returns
        -------
        None.    
        
        '''
        
        # Lists surface and checks their coincidance  
        
        if not self.__flagCheckWalls:
            self.__flagCheckWalls = True
            
            for wall1 in self.buildingSurfaces.values():
                if wall1.type=='ExtWall':
                    for wall2 in self.buildingSurfaces.values():
                        if (wall2.type=='ExtWall' and wall2 != wall1 and wall1.checkSurfaceCoincidence(wall2)):
                            
                            #print(wall1.name,wall2.name)
                            intersectionArea = wall1.calculateIntersectionArea(wall2)
                            wall1.reduceArea(intersectionArea)
                            '''
                            wall2.correctArea(wall2.Area - intersectionArea)
                            this line is not necessary, because the first for cycle will consider wall1 firstly, then after a while wall2
                            and will reduce wall1 and wall2 area in two distinct moments
                            '''
            
            
    def plotBuilding(self,addSurf):
        '''
        DEPRECATED
        '''
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.view_init(elev=90,azim=0)
        for surface in self.buildingSurfaces.values():
            if surface.type=='ExtWall':
                x_v=np.zeros(0)
                y_v=np.zeros(0)
                z_v=np.zeros(0)
                for vert in surface.vertList:
                    x_v=np.append(x_v,vert[0])
                    y_v=np.append(y_v,vert[1])
                    z_v=np.append(z_v,vert[2])
    
                ax.plot(x_v,y_v,z_v, 'b-.x')
                
        x_v=np.zeros(0)
        y_v=np.zeros(0)
        z_v=np.zeros(0)       
        for vert in addSurf.vertList:
            x_v=np.append(x_v,vert[0])
            y_v=np.append(y_v,vert[1])
            z_v=np.append(z_v,vert[2])
        ax.plot(x_v,y_v,z_v, 'r-.x')
    
    
    def printBuildingInfo(self):
        '''
        Print some info of the building
        
        Parameters
            ----------
            None
                    
        Returns
        -------
        None.    
        
        '''
        print('nome '+str(self.name)+
              '\nfootprint: '+str(self.footprint)+
              '\nN_floors: '+str(self.nFloors)+
              '\nheight: '+str(self.buildingHeight)+
              '\nexternal wall area: '+str(self.extWallArea)
              +'\n')
    
#%% ---------------------------------------------------------------------------------------------------
#%% Complex class

class Complex:
    
    ''' 
    Complex class in order to merge the output vectors of the buildings
    belonging to the same complex
    
    Methods:
        init
    '''
    
    def __init__(self, Complexname):
        '''
        Initializes the complex object
        
        Parameters
            ----------
            Complexname : string
                string with the complex name or id
        '''
    
        # Check input data type
        
        if not isinstance(Complexname, str):
            try:
                Complexname = str(Complexname)
            except ValueError:
                raise TypeError(f'Ops... Complex class, input Complexname is not a string: Complexname {Complexname}')    
        
        # Vectors inizialization 
        
        year = 8760
        self.ComplexName = Complexname
        self.heatFlowC = np.zeros(year)
        self.latentFlowC = np.zeros(year)
        self.AHUDemandC = np.zeros(year)
        self.AHUDemand_latC = np.zeros(year)
        self.AHUDemand_sensC = np.zeros(year)
        

#%% ---------------------------------------------------------------------------------------------------
        
'''
TEST METHOD
'''
import os

if __name__ == '__main__':
    buildingName = 'edificio Prova'
    archId = 1
    
    env_path = os.path.join('..','Input','buildings_envelope_V02_EP.xlsx')
    envelopes = loadEvelopes(env_path)
    
    # sched_path = os.path.join('..','Input','Schedule_ICA-RC.xlsx')
    # Archetype = pd.read_excel(sched_path,sheet_name="Archetype",header=[1],index_col=[0])
    # Arch = Archetype.loc['Office']
    # schedule = pd.read_excel(sched_path,sheet_name="Schedule",usecols="C:end",header=[2])
    
    iopath = os.path.join('..','Input', 'ITA_Venezia-Tessera.161050_IGDG.epw')
    year=2007
    epw_0=pvlib.iotools.read_epw(iopath,coerce_year=year)
    epw = epw_0[0].reset_index(drop=True)
    #time = epw[0].index
    time = np.arange(8760)
    tz='Europe/Rome'
    T_ext = epw['temp_air']
    ghi_infrared = epw['ghi_infrared']
    Sigma = 5.6697e-08
    T_as = (ghi_infrared/Sigma)**0.25-(273+T_ext)
    dT_er = statistics.mean(T_ext-T_as)
    lat, lon = epw_0[1]['latitude'], epw_0[1]['longitude']
    site = pvlib.location.Location(lat, lon, tz=tz)
    SolarGains = pd.read_csv(os.path.join('..','Input','PlanesIrradiances.csv'),header=[0,1,2],index_col=[0]).set_index(time)
    rho_air = 1.2   #kg/m3
    cp_air = 1000   #[J/(kg K)]
    
    schedpath2 = os.path.join('..','Input','Schedule_ICA-RC.xlsx')
    sched2 = loadArchetype(schedpath2,time)
    sched_zone = sched2['Office']
    
    #zone1 = ThermalZone(1, building)
    #zone2 = ThermalZone(2, building)
    
    #wall1=Surface([[0,0,0],[1,0,0],[1,0,9.5],[0,0,9.5]])
    #wall2=Surface([[0,0,0],[0,0,9.5],[0,1,9.5],[0,1,0]])
    #wall3=Surface([[0,1,0],[0,1,9.5],[1,1,9.5],[1,1,0]])
    #wall4=Surface([[1,0,0],[1,1,0],[1,1,9.5],[1,0,9.5]])
    
    #floorPT = Surface([[0,0,0],[0,1,0],[1,1,0],[1,0,0]])
    #roofP1 = Surface([[0,0,9.5],[1,0,9.5],[1,1,9.5],[0,1,9.5]])
    
    #intWall = SurfaceInternalMass(2)
    #celinig1 = SurfaceInternalAdjacent(2,adjacentZone = zone2)
    #celinig2 = SurfaceInternalAdjacent(2,adjacentZone = zone1)
    
    #surfaces1=[wall1,wall2,wall3,wall4,floorPT,intWall,celinig1]
    #surfaces2=[wall1,wall2,wall3,wall4,roofP1,intWall,celinig2]
    
    #zone1.zoneParameter(surfaces1,envelopes[archId])
    #zone2.zoneParameter(surfaces2,envelopes[archId])
    
    wall1=[[0,0,0],[1,0,0],[1,0,9.5],[0,0,9.5]]
    wall2=[[0,0,0],[0,0,9.5],[0,1,9.5],[0,1,0]]
    wall3=[[0,1,0],[0,1,9.5],[1,1,9.5],[1,1,0]]
    wall4=[[1,0,0],[1,1,0],[1,1,9.5],[1,0,9.5]]
    
    floorPT = [[0,0,0],[0,1,0],[1,1,0],[1,0,0]]
    roofP1 = [[0,0,9.5],[1,0,9.5],[1,1,9.5],[0,1,9.5]]
    
    surfaces=[wall1,wall2,wall3,wall4,floorPT,roofP1]
    
    
    edifici = {}
    
    for i in range(1):
        edifici[str(i)] = Building(buildingName, surfaces,envelopes,archId,sched2,'Office',SolarGains,T_ext,dT_er)
        #edifici[str(i)].load_calculation(SolarGains,dT_er)
  
    # for t in time:

    #     T_e = T_ext.iloc[t]
    #     sched = sched_zone.sched_df.iloc[t]
    #     Solar_gain = SolarGains.iloc[t]
    #     print(t)
        
    #     for i in range(1):
    #         #edifici['0'].load_calculation(t,Solar_gain,dT_er)
    #         edifici['0'].zones['Zone'].solve(t,T_e)
