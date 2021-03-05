'''IMPORTING MODULES'''

import numpy as np
import pandas as pd
from RC_classes.auxiliary_functions import wrn
from RC_classes.WeatherData import  Weather

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions

def PVGIS_production(tilt,GR,BR,DR,RR,Ta,W,P_nom,i):
    
    '''
    This function allows to calculate the photovoltaic electricity production
    
    Parameters
    ----------
    tilt : float
        inclination angle of the photovoltaic panel [°]
    GR : float
        global radiation on the oriented surface [W/m2]
    BR : float
        beam radiation on the oriented surface [W/m2]
    DR : float
        diffure radiation on the oriented surface [W/m2]
    RR : float
        reflected radiation on the oriented surfaces [W/m2]
    Ta : float
        external air temperature [°C]
    W : float
        wind speed [m/s]
    P_nom : float
        nominal power of the photovoltaic plant in Standard conditions
    i : float
        angle of incidence [°]
        
    Returns
    -------
    P : float
        electrical power production [W]
    etaRel : float
        plant relative efficiency [-]
    FB : float
        angular losses factor for beam radiation [-]
    
    '''
    
    # Check input data type
    
    if not isinstance(tilt, float):
        raise TypeError(f'Ops... PVGIS_production function, tilt is not a float: tilt {tilt}')
    if not isinstance(GR, float):
        raise TypeError(f'Ops... PVGIS_production function, GR is not a float: GR {GR}')
    if not isinstance(BR, float):
        raise TypeError(f'Ops... PVGIS_production function, BR is not a float: BR {BR}')
    if not isinstance(DR, float):
        raise TypeError(f'Ops... PVGIS_production function, DR is not a float: DR {DR}')
    if not isinstance(RR, float):
        raise TypeError(f'Ops... PVGIS_production function, RR is not a float: RR {RR}')
    if not isinstance(Ta, float):
        raise TypeError(f'Ops... PVGIS_production function, Ta is not a float: Ta {Ta}')
    if not isinstance(W, float):
        raise TypeError(f'Ops... PVGIS_production function, W is not a float: W {W}')
    if not isinstance(P_nom, float):
        raise TypeError(f'Ops... PVGIS_production function, P_nom is not a float: P_nom {P_nom}')
    if not isinstance(i, float):
        raise TypeError(f'Ops... PVGIS_production function, i is not a float: i {i}')
    
    # Check input data quality
    
    if not 0 <= GR <= 3000:
        wrn(f"\n\nPVGIS_production function, input GR is out of plausible range: GR {GR}\n")
    if not 0 <= BR <= 3000:
        wrn(f"\n\nPVGIS_production function, input BR is out of plausible range: BR {BR}\n")
    if not 0 <= DR <= 3000:
        wrn(f"\n\nPVGIS_production function, input DR is out of plausible range: DR {DR}\n")
    if not 0 <= RR <= 3000:
        wrn(f"\n\nPVGIS_production function, input RR is out of plausible range: RR {RR}\n")
    if not -50 <= Ta <= 60:
        wrn(f"\n\nPVGIS_production function, input Ta is out of plausible range: Ta {Ta}\n")
    if not 0 <= W <= 25:
        wrn(f"\n\nPVGIS_production function, input W is out of plausible range: W {W}\n")

    # Degree to radians conversion
    i = np.radians(i)
    tiltrad = np.radians(tilt)
    
    # Angular losses calculation (N. Martin and JM. Ruiz)
    ar = 0.157                                                                 # Reference value for Air/Glass/Silicon
    c1 = 4/(3*np.pi)
    c2 = -0.074                                                                # Reference value for Air/Glass/Silicon
    
    FB = (np.exp(-np.cos(i)/ar)-np.exp(-1/ar))/(1-np.exp(-1/ar))
    FD = np.exp(-1/ar*(c1*(np.sin(tiltrad)+(np.pi-tilt*np.pi/180-np.sin(tiltrad))/(1+np.cos(tiltrad)))+c2*(np.sin(tiltrad)+(np.pi-tilt*np.pi/180-np.sin(tiltrad))/(1+np.cos(tiltrad)))**2))
    FA = np.exp(-1/ar*(c1*(np.sin(tiltrad)+(tilt*np.pi/180-np.sin(tiltrad))/(1-np.cos(tiltrad)))+c2*(np.sin(tiltrad)+(tilt*np.pi/180-np.sin(tiltrad))/(1-np.cos(tiltrad)))**2))
    
    # Setting curve coefficients (Reference values for c-Si modules)
    k1 = -0.017237
    k2 = -0.040465
    k3 = -0.004702
    k4 = 0.000149
    k5 = 0.000170
    k6 = 0.000005
    
    # Module temperature evaluation (c-Si modules)
    U0 = 26.92
    U1 = 6.24
    Tm = Ta + GR/(U0+U1*W)                                                     # Module temperature [°C]
    tm = Tm - 25
    
    # Relative efficiency and power production calculation
    g = GR/1000
    if g == 0:
        etaRel = 0
    else:
        etaRel = 1+k1*np.log(g)+k2*np.log(g)**2+k3*tm+k4*tm*np.log(g)+k5*tm*np.log(g)**2+k6*tm**2
    P = 1/1000*P_nom*etaRel*(BR*(1-FB)+DR*(1-FD)+RR*(1-FA))                    # [W]
    if P < 0:
        P = 0
    return P, etaRel, FB                                                       # Output variables
    

#%% ---------------------------------------------------------------------------------------------------
#%% DistrictPVGIS class

class DistrictPVGIS:
    '''
    DistrictPVGIS class is used to determine the photovoltaic electricity
    production of a reference tilted photovoltaic plant on buildings roof
    
    Methods
        init
        PV_calc
    '''
    
    # class variables
    Tilt = 30.                                                                 # [°]
    Az = 0.                                                                    # [°]
    etaNom = 0.15                                                              # Standard condition efficiency [-]
    area = 0.4                                                                 # Percentage of roof area
    rho_g = 0.2                                                                # Albedo coefficient
    
    def __init__(self,ExtRoofArea,l):
        '''
        Initilization of the photovoltaic plant
        
        Parameters
        ----------
        ExtRoofArea : float
            building's external roof area [m2]
        l : int
            number of simulation time steps [-]
        
        Returns
        -------
        None
        
        '''
        
        # Check input data 
        
        if not isinstance(ExtRoofArea, float):
            raise TypeError(f'Ops... DistrictPVGIS class, ExtRoofArea is not a float: ExtRoofArea {ExtRoofArea}')
        if not isinstance(l, int):
            raise TypeError(f'Ops... DistrictPVGIS class, l is not a int: l {l}')
        
        # Check input data quality
        
        if ExtRoofArea < 0:
            wrn(f"\n\nDistrictPVGIS class, input ExtRoofArea must be positive: ExtRoofArea {ExtRoofArea}\n")
        if l < 0:
            wrn(f"\n\nDistrictPVGIS class, input l must be positive: l {l}\n")
        
        # Total area and nominal power evaluation
        self.A_all = ExtRoofArea*self.area                                     # Total plant area [m2]
        self.P_nom = 1000*self.A_all*self.etaNom                               # Nominal power in Standard condition [W]
        
        # Output vectors initialization
        self.P = np.zeros(l)
        self.etaRel = np.zeros(l)
        self.FB = np.zeros(l)  
        self.P_el = np.zeros(l)
        self.l = l
        

    def PV_calc(self,weather):
        '''
        PV_calc calculate the electricity production
        
        Parameters
        ----------
        weather : RC_classes.WeatherData.Weather obj
            object of the class weather WeatherData module
        
        Returns
        -------
        None
        '''

        # Check input data 
        
        if not isinstance(weather, Weather):
            raise TypeError(f'Ops... JsonCity class, weather is not a RC_classes.WeatherData.Weather: weather {weather}')
        
        # Irradiance on the tilted surface
        Tilt_rad = np.radians(self.Tilt)                                             # Degree to radians conversion
        irradiance = weather.SolarGains[str(self.Az)][str(self.Tilt)]
        irradiance0 = weather.SolarGains[str(0.0)][str(0.0)]
        self.AOI = irradiance['AOI'].to_numpy()                                      # Angle of incidence [°]
        self.BRTS = irradiance['direct'].to_numpy()                                  # Direct radiation [W/m2]
        self.DRTS = irradiance['global'].to_numpy()-irradiance['direct'].to_numpy()  # Diffuse radiation [W/m2]
        TRH = irradiance0['global'].to_numpy()                                       # Global radiation on horizontal surface [W/m2]
        self.RRTS = TRH*0.5*(1-np.cos(Tilt_rad))*self.rho_g                          # Reflected radiation [W/m2]
        TRTS = self.BRTS + self.RRTS + self.DRTS                                     # Global radiation [W/m2]
        time = np.arange(self.l)
        
        # Electrical power calculation
        for t in time:
            self.P[t], self.etaRel[t], self.FB[t] = PVGIS_production(self.Tilt,TRTS[t],self.BRTS[t],self.DRTS[t],self.RRTS[t],weather.Text[t],weather.w[t],self.P_nom,self.AOI[t])
            self.P_el[t] = self.P[t]/1000                                      # [kW]
        self.Tot_prod = sum(self.P_el)/1000                                    # [MWh]