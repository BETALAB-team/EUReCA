'''IMPORTING MODULES'''

import pvlib
import os
import pandas as pd
import numpy as np
from warnings import warn as wrn

#%% ---------------------------------------------------------------------------------------------------

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
    if not isinstance(time, pd.core.series.Series):
        raise TypeError(f'Ops... input time is not a pandas series: time {time}') 
    if not isinstance(solar_position, pd.core.frame.DataFrame):
        raise TypeError(f'Ops... input solar_position is not a pandas dataframe: solar_position {solar_position}') 
    if not isinstance(surf_tilt, float):
        raise TypeError(f'Ops... input surf_tilt is not a float: surf_tilt {surf_tilt}') 
    if not isinstance(surf_az, float):
        raise TypeError(f'Ops... input surf_az is not a float: surf_az {surf_az}') 
    if not isinstance(irradiance, pd.core.frame.DataFrame):
        raise TypeError(f'Ops... input irradiance is not a pandas dataframe: irradiance {irradiance}') 
    if not isinstance(year, int):
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

#%%--------------------------------------------------------------------------------------------------- 

'''
TEST METHOD
'''
if __name__ == '__main__':
    iopath = os.path.join('..','Input', 'ITA_Venezia-Tessera.161050_IGDG.epw')
    year=2007
    epw=pvlib.iotools.read_epw(iopath,coerce_year=year)
    tz='Europe/Rome'
    lat, lon = epw[1]['latitude'], epw[1]['longitude']
    site = pvlib.location.Location(lat, lon, tz=tz)
    Irradiances=PlanesIrradiances(site,epw,year,8,3)
    Irradiances.Irradiances.to_csv(os.path.join('..','Input','PlanesIrradiances.csv')) 
