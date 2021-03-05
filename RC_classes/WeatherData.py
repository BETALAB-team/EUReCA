'''IMPORTING MODULES'''

import numpy as np
import statistics
import pandas as pd
from warnings import warn as wrn
import os
import pvlib

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


def get_irradiance(site,time,solar_position, surf_tilt, surf_az,irradiance,year):
    '''
    function from pvlib to calculate irradiance on a specific surface
    https://pvlib-python.readthedocs.io/en/stable/auto_examples/plot_ghi_transposition.html#sphx-glr-auto-examples-plot-ghi-transposition-py    
    
    Parameters
    ----------
    site : pvlib.location.Location
        Location object from pvlib.
    time: pandas series
        time index of the epw object from pvlib
    solar_position : pandas dataframe 
            the solar position dataframe from pvlib
    surf_tilt: float
        tilt of the surface
    surf_az: float
        azimuth of the surface
    irradiance:   pandas dataframe  
        dataframe from the epw with ghi, dni, dhi
    year : int
        year of the simulation, just to set the dataframe index
        
    Returns
    -------
    pandas DataFrame with the irradiances on the surface 
    '''
    
    # Check input data type
        
    if not isinstance(site, pvlib.location.Location):
        raise TypeError(f'Ops... input site is not a pvlib location object: site {site}') 
    if not isinstance(time, pd.core.indexes.datetimes.DatetimeIndex):
        raise TypeError(f'Ops... input time is not a pandas series: time {time}') 
    if not isinstance(solar_position, pd.core.frame.DataFrame):
        raise TypeError(f'Ops... input solar_position is not a pandas dataframe: solar_position {solar_position}') 
    if not isinstance(surf_tilt, float):
        try:
            surf_tilt = float(surf_tilt)
        except:
            raise TypeError(f'Ops... input surf_tilt is not a float: surf_tilt {surf_tilt}') 
    if not isinstance(surf_az, float):
        try:
            surf_az = float(surf_az)
        except:
            raise TypeError(f'Ops... input surf_az is not a float: surf_az {surf_az}') 
    if not isinstance(irradiance, pd.core.frame.DataFrame):
        raise TypeError(f'Ops... input irradiance is not a pandas dataframe: irradiance {irradiance}') 
    if not isinstance(year, int):
        try:
            year = int(year)
        except:
            raise TypeError(f'Ops... input year is not an integer: year {year}') 
    
    # Control input data quality
    
    if surf_tilt < 0 or surf_tilt > 90 or  surf_az < -180 or surf_az > 180:
        wrn(f"\n\nget_irradiance funtion, are you sure that the surface orientation is correct?? surf_tilt {surf_tilt}, surf_az {surf_az}\n")

    # Use pvlib function to calculate the irradiance on the surface
    
    surf_az = surf_az+180
    POA_irradiance = pvlib.irradiance.get_total_irradiance(
        surface_tilt=surf_tilt,
        surface_azimuth=surf_az,
        dni_extra=pvlib.irradiance.get_extra_radiation(time,solar_constant=1366.1, method='spencer', epoch_year=year),
        dni=irradiance['dni'],
        ghi=irradiance['ghi'],
        dhi=irradiance['dhi'],
        solar_zenith=solar_position['apparent_zenith'],
        solar_azimuth=solar_position['azimuth'],
        model='isotropic',
        model_perez='allsitescomposite1990',
        airmass=site.get_airmass(solar_position=solar_position))
    AOI = pvlib.irradiance.aoi(
        surface_tilt=surf_tilt,
        surface_azimuth=surf_az,
        solar_zenith=solar_position['apparent_zenith'],
        solar_azimuth=solar_position['azimuth'])
    
    # Cleaning AOI vector
    
    for i in range(len(AOI)): 
        if AOI[i] > 90 or solar_position['apparent_zenith'][i] > 90:
            AOI[i] = 90
    
    return pd.DataFrame({'GHI': irradiance['ghi'],
                         'POA': POA_irradiance['poa_global'],
                         'POA_B': POA_irradiance['poa_direct'],
                         'POA_D': POA_irradiance['poa_global']-POA_irradiance['poa_direct'],
                         'AOI':AOI,
                         'solar zenith':solar_position['apparent_zenith']})

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
        
#%%--------------------------------------------------------------------------------------------------- 
#%% PlanesIrradiances class

class PlanesIrradiances:
    
    '''
    This class is a preprocessor of solar radiation
    Estimates the annual huorly solar radiation for a specific set of decided surfaces
    
    When initialized, the number of azimuth orientations and the number of height orientation must be given
    Solar irradiance will be calculated for the entire set of these orientations
    
    example: if number of azimuth subdivision is 4 solar irradiance will be calculate for N E S W

    Methods:
        init
    '''
        
    def __init__(self,site,epw,year = 2020,azSubdiv = 8,hSubdiv = 3):
    
        '''
        initialize PlaneIrradiances object
        
        Parameters
        ----------
        site : pvlib.location.Location
            Location object from pvlib.
        epw : tuple from pvlib read_epw
            contains 2 dataframe with the epw characteristics
        year : int
            year of the simulation, just to set the dataframe index
        azSubdiv: int
            number of the different direction (azimuth) solar radiation will be calculated
        hSubdiv: int
            number of the different direction (solar height) solar radiation will be calculated      
            
        Returns
        -------
        None
        '''
        
        # Check input data type
        
        if not isinstance(site, pvlib.location.Location):
            raise TypeError(f'Ops... input site is not a pvlib location object: site {site}') 
        if not isinstance(epw, tuple):
            raise TypeError(f'Ops... input epw is not a tuple: epw {epw}') 
        if not isinstance(epw[0], pd.core.frame.DataFrame):
            raise TypeError(f'Ops... input epw[0] is not a pandas dataframe: epw[0] {epw[0]}') 
        if not isinstance(epw[1], dict):
            raise TypeError(f'Ops... input epw[1] is not a dictionary: epw[1] {epw[1]}') 
        if not isinstance(year, int):
            raise TypeError(f'Ops... input year is not an integer: year {year}') 
        if not isinstance(azSubdiv, int):
            raise TypeError(f'Ops... input azSubdiv is not an integer: azSubdiv {azSubdiv}') 
        if not isinstance(hSubdiv, int):
            raise TypeError(f'Ops... input hSubdiv is not an integer: hSubdiv {hSubdiv}') 
        
        # Control input data quality
    
        if azSubdiv > 10 or hSubdiv > 5:
            wrn(f"\n\nPlanesIrradiances class, init, solar calculation could be long..... azSubdiv {azSubdiv}, hSubdiv {hSubdiv}\n")
     
        # Creates the dataframe with global beam and Angle of Incidence for many directions
        
        time = epw[0].index
        irradiance = epw[0][['ghi','dni','dhi']]
        self.solar_position = site.get_solarposition(times=time)
        
        azimuth = np.linspace(-180,180,azSubdiv+1)[:-1]
        height = np.linspace(90,0,hSubdiv+1)[:-1]
        components = ['global','direct','AOI']
        
        col = pd.MultiIndex.from_product([azimuth,height,components], names=['Azimuth', 'Height','Component'])
        col = col.union(pd.MultiIndex.from_product([[0],[0],components], names=['Azimuth', 'Height','Component']))
        Irradiances = pd.DataFrame(0., index=time, columns=col)
        
        for az in azimuth:
            for h in height:
                POA = get_irradiance(site,time,self.solar_position, h, az,irradiance,year)
                Irradiances[az,h,'global']= POA['POA']
                Irradiances[az,h,'direct']= POA['POA_B']
                Irradiances[az,h,'AOI']= POA['AOI']
                
        POA = get_irradiance(site,time,self.solar_position, 0, 0,irradiance,year)
        
        Irradiances[0,0,'global']= POA['POA']
        Irradiances[0,0,'direct']= POA['POA_B']
        Irradiances[0,0,'AOI']= POA['AOI']
        
        self.Irradiances=Irradiances

#%% ---------------------------------------------------------------------------------------------------
# Weather class
class Weather():
    '''
    This class is a container for all weather data.
    It processes the epw file to extract arrays of temperature, wind, humidity etc......
    
    Methods:
        init
    '''

    def __init__(self,epw_name, input_path = '.', tz = 'Europe/Rome', 
                 year = 2020, ts = 2, hours = 8760, n_years = 1, 
                 irradiances_calc = True, 
                 azSubdiv = 8, hSubdiv = 3, 
                 shad_tol = [80.,
                             100.,
                             80.]):
        '''
        initialize weather obj
        It processes the epw file to extract arrays of temperature, wind, humidity etc......
    
        Parameters
        ----------
        epw_name : str
            path of the epw file.
        input_path : str
            path of the Input files.
        tz : str
            this is the time zone. For a list of the avilable see: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
        year : int 
            the year of simulation. it is used only to create a pd.DataFrame.
        ts : int
            number of time steps in a hour.
        hours : int 
            number of hours
        n_years : int
            number of years of sim
        azSubdiv: int
            number of the different direction (azimuth) solar radiation will be calculated
        hSubdiv: int
            number of the different direction (solar height) solar radiation will be calculated      
        shad_tol: list of floats
            list of the tollerances for urban shading calc (azimuth, disance, theta)
            
        Returns
        -------
        None
        '''
        
        # Check input data type
        
        if not isinstance(epw_name, str):
            raise TypeError(f'Weather object, init, input epw_file is not a string: epw_file {epw_file}') 
        if not isinstance(input_path, str):
            raise TypeError(f'Weather object, init, input input_path is not a string: input_path {input_path}') 
        if not isinstance(tz, str):
            raise TypeError(f'Weather object, init, input tz is not a string: tz {tz}') 
        if not isinstance(year, int):
            try:
                year = int(year)
            except:
                raise TypeError(f'Weather object, init, input year is not an int: year {year}') 
        if not isinstance(ts, int):
            try:
                ts = int(ts)
            except:
                raise TypeError(f'Weather object, init, input ts is not an int: ts {ts}') 
        if not isinstance(hours, int):
            try:
                hours = int(hours)
            except:
                raise TypeError(f'Weather object, init, input hours is not an int: hours {hours}') 
        if not isinstance(n_years, int):
            try:
                n_years = int(n_years)
            except:
                raise TypeError(f'Weather object, init, input n_years is not an int: n_years {n_years}') 
        if not isinstance(azSubdiv, int):
            raise TypeError(f'Weather object, init, input azSubdiv is not an integer: azSubdiv {azSubdiv}') 
        if not isinstance(hSubdiv, int):
            raise TypeError(f'Weather object, init, input hSubdiv is not an integer: hSubdiv {hSubdiv}') 
        if not isinstance(shad_tol,list):
            raise TypeError(f'Weather object, init, input shad_tol is not list: shad_tol {shad_tol}')
        if not isinstance(shad_tol[0],float) or not isinstance(shad_tol[1],float) or not isinstance(shad_tol[2],float):
            try:
                shad_tol[0] = float(shad_tol[0])
                shad_tol[1] = float(shad_tol[1])
                shad_tol[2] = float(shad_tol[2])
            except:
                raise TypeError(f'Weather object, init, input shad_tol is not a list of floats: shad_tol {shad_tol}') 
            
        # Check input data quality
        
        if ts > 4:
            wrn(f"\n\nWeather object, init,  input ts is higher than 4, this means more than 4 time steps per hours were set: ts {ts}\n")
        if azSubdiv > 10 or hSubdiv > 5:
            wrn(f"\n\nWeather object, init,  solar calculation could be long..... azSubdiv {azSubdiv}, hSubdiv {hSubdiv}\n")
        if shad_tol[0] < 0.0 or shad_tol[0] > 90.0:
            wrn(f"\n\Weather object, init,  toll_az is out of range [0-90]..... toll_az {shad_tol[0]}\n")
        if shad_tol[1] < 0.0 or shad_tol[1] > 200.0:
            wrn(f"\n\Weather object, init,  toll_dist is out of range [0-200]..... toll_dist {shad_tol[1]}\n")
        if shad_tol[2] < 0.0 or shad_tol[2] > 90.0:
            wrn(f"\n\Weather object, init,  toll_theta is out of range [0-90]..... toll_theta {shad_tol[2]}\n")
                    
        epw_file = os.path.join(input_path,epw_name)
     
        # Importing and processing weather data from .epw 
        try:
            epw = pvlib.iotools.read_epw(epw_file, coerce_year = year)                       # Reading the epw via pvlib 
        except FileNotFoundError:
            raise FileNotFoundError(f"Weather epw file not found in the Input folder: epw name {epw_name}, input folder {input_path}")
        
        epw_res = epw[0].reset_index(drop=True)                                        # Exporting the hourly values
        lat, lon = epw[1]['latitude'], epw[1]['longitude']                             # Extracting latitude and longitude from the epw
        site = pvlib.location.Location(lat, lon, tz = tz)                              # Creating a location variable
        time = np.arange(8760)                                                         # Time vector inizialization
        
        # Weather Data and Average temperature difference between Text and Tsky 
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
        
        # Inizialization of Solar Gain DataFrame 
        # Creates a dataframe with the hourly solar radiation ond angle of incidence for several directions 
        
        if irradiances_calc:
            Irradiances=PlanesIrradiances(site,epw,year,azSubdiv,hSubdiv)
            Irradiances.Irradiances.to_csv(os.path.join(input_path,'PlanesIrradiances.csv'))
        
        # Solar Gain DataFrame 
        Solar_Gains = pd.read_csv(os.path.join(input_path,'PlanesIrradiances.csv'),header=[0,1,2],index_col=[0])
        Solar_Gains = rescale_sol_gain(ts,Solar_Gains)
        
        # Check some weather data values
        if not np.all(np.greater(T_ext,-50.)) or not np.all(np.less(T_ext,60.)):
            wrn(f"\n\nWeather class, input T_ext is out of plausible range: T_ext {T_ext}\n")
        if not np.all(np.greater(w,-0.001)) or not np.all(np.less(w,25.001)):
            wrn(f"\n\nWeather class, input w is out of plausible range: w {w}\n")       
        if not np.all(np.greater(RH_ext,-0.0001)) or not np.all(np.less(RH_ext,1.)):
            wrn(f"\n\nWeather class, input w is out of plausible range: w {w}\n")       
        
        # Memorizing the attributes needed for the sim        
        self.ts = ts                    # number of timesteps in an hour
        self.hours = hours              # number of hours of a sim
        self.sim_time = [ts,hours]      #
        self.tau = 3600/ts              # number of seconds of a time step (for the solution of the networks)
        self.Text = T_ext               # External air temperature [°C]
        self.RHext = RH_ext             # External relative humidity [0-1]
        self.azSubdiv = azSubdiv        # Number of azimuth subdivision [-]
        self.hSubdiv = hSubdiv          # Number of height subdivision [-]
        self.w = w                      # wind velocity [m/s]
        self.dT_er = dT_er              # Average temperature diff between external air and sky vault[°C]
        self.THextAv = T_ext_H_avg      # Average heating season air temperature [°C]
        
        self.SolarPosition = Solar_position
        self.SolarGains = Solar_Gains
        
        # tollerances for urban shading calculation
        self.Az_toll = shad_tol[0]      # [°]
        self.Dist_toll = shad_tol[1]    # [m]
        self.Theta_toll = shad_tol[2]   # [°]