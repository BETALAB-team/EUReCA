'''IMPORTING MODULES'''

import numpy as np
import statistics
import pandas as pd
from warnings import warn as wrn

#%% ---------------------------------------------------------------------------------------------------

def TskyCalc(T_ext,T_dp,P_,n_opaque):
    '''
    Apparent sky temperature calculation procedure
    Martin Berdhal model used by TRNSYS 
    
    Parameters
    ----------
    T_ext : dataframe column
        External Temperature [°C]
    T_dp : dataframe column
        External dew point temperature [°C]
    P_ : dataframe column
        External pressure [Pa]
    n_opaque : dataframe column
        Opaque sky covering [-]

    Returns
    -------
    dTer: int, Average temperature difference between External air temperature and Apparent sky temperature [°C]

    '''
    
    # Check input data type
    
    if not isinstance(T_ext, pd.core.series.Series):
        raise TypeError(f'Ops... input T_ext is not a pandas series: T_ext {T_ext}')    
    if not isinstance(T_dp, pd.core.series.Series):
        raise TypeError(f'Ops... input T_dp is not a pandas series: T_dp {T_dp}')    
    if not isinstance(P_, pd.core.series.Series):
        raise TypeError(f'Ops... input P_ is not a pandas series: P_ {P_}')    
    if not isinstance(n_opaque, pd.core.series.Series):
        raise TypeError(f'Ops... input n_opaque is not a pandas series: n_opaque {n_opaque}')    
    
    # Calculation Martin Berdhal model used by TRNSYS 
        
    day = np.arange(24)                                                        # Inizialization day vector
    T_sky = np.zeros(24)
    Tsky = []
    T_sky_anno = []
    for i in range(365):
        for x in day:
            t = i*24 + x
            Tdp = T_dp.iloc[t]
            P = P_.iloc[t]/100                                                 # [mbar]
            nopaque = n_opaque.iloc[t]*0.1                                     # [0-1]
            eps_m = 0.711+0.56*Tdp/100 + 0.73*(Tdp/100)**2
            eps_h = 0.013*np.cos(2*np.pi*(x+1)/24)
            eps_e = 0.00012*(P-1000)
            eps_clear = eps_m + eps_h + eps_e                                  # Emissivity under clear sky condition
            C = nopaque*0.9                                                    # Infrared cloud amount
            eps_sky = eps_clear + (1 - eps_clear)*C                            # Sky emissivity
            T_sky[x] = ((T_ext.iloc[t]+273)*(eps_sky**0.25))-273               # Apparent sky temperature [°C]
        Tsky = np.append(Tsky,T_sky)                                           # Annual Tsky created day by day
    T_sky_anno.append(Tsky)
                                            
    # Average temperature difference between External air temperature and Apparent sky temperature
    dT_er = statistics.mean(T_ext-Tsky)
    
    return dT_er

def rescale_weather(ts,w,T_ext,RH_ext,T_dp,P_,n_opaque):
    '''
    rescales weather variables arrays
    returns numpy arrays!!!!
    
    Parameters
    ----------
    ts : int
        Number of time steps per hour.
    w : dataframe column
        External wind speed [m/s]
    T_ext : dataframe column
        External dry bulb temperature [°C]
    RH_ext : dataframe column
        External relative humidity [0-1]    
    T_dp : dataframe column
        External dew point temperature [°C]
    P_ : dataframe column
        External pressure [Pa]
    n_opaque : dataframe column
        Opaque sky covering [-]

    Returns
    -------
    w: np array, wind speed [m/s]
    t: np array, temperature [°C]
    rh: np array, relative humidity [0-1] 
    dp: np array, dew point temperature [°C]
    p: np array, pressure [Pa]
    nop: np array, opaque sky covering [-]   
    
    '''
    
    # Check input data type
    
    if not isinstance(ts, int):
        raise TypeError(f'Ops... input ts is not an integer: ts {ts}') 
    if not isinstance(w, pd.core.series.Series):
        raise TypeError(f'Ops... input w is not a pandas series: w {w}')  
    if not isinstance(T_ext, pd.core.series.Series):
        raise TypeError(f'Ops... input T_ext is not a pandas series: T_ext {T_ext}')  
    if not isinstance(RH_ext, pd.core.series.Series):
        raise TypeError(f'Ops... input RH_ext is not a pandas series: RH_ext {RH_ext}')     
    if not isinstance(T_dp, pd.core.series.Series):
        raise TypeError(f'Ops... input T_dp is not a pandas series: T_dp {T_dp}')    
    if not isinstance(P_, pd.core.series.Series):
        raise TypeError(f'Ops... input P_ is not a pandas series: P_ {P_}')    
    if not isinstance(n_opaque, pd.core.series.Series):
        raise TypeError(f'Ops... input n_opaque is not a pandas series: n_opaque {n_opaque}')     
    
    # Check input data quality
    
    if not np.all(np.greater(w,-0.001)) or not np.all(np.less(w,25.001)):
        wrn(f"\n\nrescale_weather function, input w is out of plausible range: w {w}\n")
    
    # Rescale the data with pandas rescale
    
    data = pd.DataFrame(zip(w,T_ext,RH_ext,T_dp,P_,n_opaque), columns = ['w','t','rh','dp','p','nop'])
    m = str(60/float(ts)) + 'min'
    time = pd.date_range('2020-01-01', periods=8760, freq='1h')
    data.set_index(time, inplace=True) 
    data = data.resample(m).interpolate(method='linear')                               # Steps interpolation Resampling
    return data['w'].to_numpy() , data['t'].to_numpy() , data['rh'].to_numpy() , data['dp'].to_numpy() , data['p'].to_numpy() , data['nop'].to_numpy()

def rescale_sol_gain(ts,df):
    '''
    rescales solar gains dataframe
    
    Parameters
    ----------
    ts : int
        Number of time steps per hour.
    df : pandas dataframe 
        the irradiances dataframe
    
    Returns
    -------
    : pd dataframe of irradiances
    '''
    
    # Check input data type
    
    if not isinstance(ts, int):
        raise TypeError(f'Ops... input ts is not an integer: ts {ts}') 
    if not isinstance(df, pd.core.frame.DataFrame):
        raise TypeError(f'Ops... input df is not a pandas dataframe: df {df}') 

    # Rescale the data with pandas rescale
        
    m = str(60/float(ts)) + 'min'
    time = pd.date_range('2020-01-01', periods=8760, freq='1h')
    df.set_index(time, inplace=True) 
    return df.resample(m).interpolate(method='linear').reset_index(drop=True)                               # Steps interpolation Resampling

#%% ---------------------------------------------------------------------------------------------------
# SolarPosition class
class SolarPosition():
    '''
    This class only contains the solar position during the simulation run period 
    
    Methods:
        init
    '''

    def __init__(self,ts,sp):
        '''
        initialize SolarPosition object
        
        Parameters
        ----------
        ts : int
            Number of time steps per hour.
        sp : pandas dataframe 
            the solar position dataframe from pvlib
        
        Returns
        -------
        None
        '''
        
        # Check input data type
        
        if not isinstance(ts, int):
            raise TypeError(f'Ops... input ts is not an integer: ts {ts}') 
        if not isinstance(sp, pd.core.frame.DataFrame):
            raise TypeError(f'Ops... input sp is not a pandas dataframe: sp {sp}') 
            
        # Rescale the data with pandas rescale    
    
        m = str(60/float(ts)) + 'min'
        time = pd.date_range('2020-01-01', periods=8760, freq='1h')
        sp.set_index(time, inplace=True) 
        sp = sp.resample(m).interpolate(method='linear')                               # Steps interpolation Resampling
        
        # Set some attributes of the class
        
        self.apparent_zenith = sp['apparent_zenith'].to_numpy()
        self.zenith = sp['zenith'].to_numpy()
        self.apparent_elevation = sp['apparent_elevation'].to_numpy()
        self.elevation = sp['elevation'].to_numpy()
        self.azimuth = sp['azimuth'].to_numpy() - 180
        self.EOfTime = sp['equation_of_time'].to_numpy()