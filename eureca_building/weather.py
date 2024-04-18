"""This module includes classes and functions to manage weather file
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging
import io

import pvlib
import requests
import numpy as np
import pandas as pd

from eureca_building.config import CONFIG
from eureca_building._auxiliary_function_for_monthly_calc import get_monthly_value_from_annual_vector

# %% ---------------------------------------------------------------------------------------------------
# Weather class
class WeatherFile():
    '''This class is a container for all weather data.
    It processes the epw file to extract arrays of temperature, wind, humidity etc......
    '''

    def __init__(self,
                 epw: str,
                 year=None,
                 time_steps: int = 1,
                 irradiances_calculation: bool = True,
                 azimuth_subdivisions: int = 8,
                 height_subdivisions: int = 3,
                 urban_shading_tol=[80., 100., 80.]
                 ):
        '''Initialize weather obj. It processes the epw file to extract arrays of temperature, wind, humidity etc......

        Parameters
        ----------
        epw_name : str
            path of the epw file.
        year : int, default None
            the year of simulation. it is used only to create a pd.DataFrame.
        time_steps : int, default 1
            number of time steps in a hour.
        irradiances_calculation : bool, default False
            Whether to do or not the irradiances calculation
        azimuth_subdivisions : int, default 8
            number of the different direction (azimuth) solar radiation will be calculated
        height_subdivisions : int, default 3
            number of the different direction (solar height) solar radiation will be calculated
        urban_shading_tol: list, default [80.,100.,80.]
            list of three floats with the tolerances for urban shading calc (azimuth, distance, theta)
        '''
        self.weatherepw=epw
        # Importing and processing weather data from .epw
        try:
            epw = pvlib.iotools.read_epw(epw, coerce_year=year)  # Reading the epw via pvlib
        except OSError:
            epw = pvlib.iotools.parse_epw(epw)
        except FileNotFoundError:
            raise FileNotFoundError \
                (f"ERROR Weather epw file not found in the Input folder.")

        # Check some integer inputs
        if not isinstance(time_steps, int):
            raise TypeError(f"WeatherFile: time_steps must be a integer: time_steps = {time_steps}")
        if not isinstance(azimuth_subdivisions, int):
            raise TypeError(
                f"WeatherFile: azimuth_subdivisions must be a integer: azimuth_subdivisions = {azimuth_subdivisions}")
        if not isinstance(height_subdivisions, int):
            raise TypeError(f"WeatherFile: time_steps must be a integer: height_subdivisions = {height_subdivisions}")
        if time_steps > 6:
            logging.warning(f"WeatherFile: time_steps higher than 6")
        if azimuth_subdivisions > 10:
            logging.warning(f"WeatherFile: azimuth_subdivisions higher than 10")
        if height_subdivisions > 4:
            logging.warning(f"WeatherFile: height_subdivisions higher than 4")

        self._epw_hourly_data = epw[0]  # Hourly values
        self._epw_general_data = epw[1]  # Extracting general data
        self._site = pvlib.location.Location(self._epw_general_data['latitude'],
                                             self._epw_general_data['longitude'],
                                             tz=self._epw_general_data['TZ'])  # Creating a location variable

        if time_steps > 1:
            m = str(60 / float(time_steps)) + 'min'
            self._epw_hourly_data = epw[0].resample(m).interpolate(method='linear')

        # Weather Data and Average temperature difference between Text and Tsky
        self.hourly_data = {}
        self.general_data = {}
        self.hourly_data["time_index"] = self._epw_hourly_data.index
        self.hourly_data["wind_speed"] = self._epw_hourly_data['wind_speed'].values  # [m/s]
        self.hourly_data["wind_direction"] = self._epw_hourly_data['wind_direction'].values  # [°]
        self.hourly_data["out_air_db_temperature"] = self._epw_hourly_data['temp_air'].values  # [°C]
        self.hourly_data["out_air_dp_temperature"] = self._epw_hourly_data['temp_dew'].values  # [°C]
        self.hourly_data["out_air_relative_humidity"] = self._epw_hourly_data['relative_humidity'].values / 100  # [0-1]
        self.hourly_data["out_air_pressure"] = self._epw_hourly_data['atmospheric_pressure'].values  # Pa
        self.hourly_data["opaque_sky_coverage"] = self._epw_hourly_data['opaque_sky_cover'].values  # [0-10]
        # Average temperature difference between Text and Tsky
        self.general_data['average_dt_air_sky'] = _TskyCalc(self.hourly_data["out_air_db_temperature"],
                                                            self.hourly_data["out_air_dp_temperature"],
                                                            self.hourly_data["out_air_pressure"],
                                                            self.hourly_data["opaque_sky_coverage"],
                                                            time_steps)
        self.general_data['number_of_time_steps'] = len(self._epw_hourly_data.index)
        self.general_data['time_steps_per_hour'] = time_steps
        self.general_data['azimuth_subdivisions'] = azimuth_subdivisions
        self.general_data['height_subdivisions'] = height_subdivisions

        # Check some weather data values
        if not np.all(np.greater(self.hourly_data["out_air_db_temperature"], -50.)) or not np.all(
                np.less(self.hourly_data["out_air_db_temperature"], 60.)):
            logging.warning(f"WeatherFile class, input drybulb outdoor temperature is out of range [-50:60] °C")
        if not np.all(np.greater(self.hourly_data["wind_speed"], -0.001)) or not np.all(
                np.less(self.hourly_data["wind_speed"], 25.001)):
            logging.warning(f"WeatherFile, input wind speed is out of range [-0.001, 25.] m/s")
        if not np.all(np.greater(self.hourly_data["out_air_relative_humidity"], -0.0001)) or not np.all(
                np.less(self.hourly_data["out_air_relative_humidity"], 1.)):
            logging.warning(f"WeatherFile, input relative humidity is out of range [-0.001, 1] [-]")

        # Humidity calculation
        self.hourly_data["out_air_saturation_pressure"] = np.zeros(len(self.hourly_data["out_air_db_temperature"]))
        higher_filt = self.hourly_data["out_air_db_temperature"] > 0
        lower_filt = np.logical_not(higher_filt)
        self.hourly_data["out_air_saturation_pressure"][lower_filt] = 610.5 * np.exp(
            (21.875 * self.hourly_data["out_air_db_temperature"][lower_filt]) / (
                    265.5 + self.hourly_data["out_air_db_temperature"][lower_filt]))
        self.hourly_data["out_air_saturation_pressure"][higher_filt] = 610.5 * np.exp(
            (17.269 * self.hourly_data["out_air_db_temperature"][higher_filt]) / (
                    237.3 + self.hourly_data["out_air_db_temperature"][higher_filt]))

        self.hourly_data["out_air_specific_humidity"] = 0.622 * (
                self.hourly_data["out_air_relative_humidity"] * self.hourly_data["out_air_saturation_pressure"] / (
                self.hourly_data["out_air_pressure"] - (self.hourly_data["out_air_relative_humidity"] *
                                                        self.hourly_data[
                                                            "out_air_saturation_pressure"])))

        if irradiances_calculation:
            self.irradiances_calculation()

        self.general_data['urban_shading_tol'] = urban_shading_tol

        # Definition of average T_ext during heating season
        av_t_daily = np.array([np.mean(self.hourly_data["out_air_db_temperature"][i:i+24*time_steps]) for i in range(0,8760 * time_steps, 24*time_steps)])
        self.general_data['heating_degree_days'] =  np.sum(20-av_t_daily[av_t_daily < 12])
        self.general_data['average_out_air_db_temperature_heating_season'] = np.mean(np.hstack([
            self.hourly_data["out_air_db_temperature"][CONFIG.heating_season_start_time_step:],
            self.hourly_data["out_air_db_temperature"][:CONFIG.heating_season_end_time_step]
        ]))
        self.general_data['average_out_air_db_temperature_cooling_season'] = np.mean(self.hourly_data["out_air_db_temperature"][CONFIG.cooling_season_start_time_step:CONFIG.cooling_season_end_time_step])
        self.general_data['average_out_air_db_temperature'] = np.mean(self.hourly_data["out_air_db_temperature"])

        self.monthly_data = {}
        self.monthly_data["out_air_specific_humidity"] = get_monthly_value_from_annual_vector(self.hourly_data["out_air_specific_humidity"], method='mean')
        self.monthly_data["out_air_relative_humidity"] = get_monthly_value_from_annual_vector(self.hourly_data["out_air_relative_humidity"], method='mean')
        self.monthly_data["out_air_db_temperature"] = get_monthly_value_from_annual_vector(self.hourly_data["out_air_db_temperature"], method='mean')

    def irradiances_calculation(self):
        """Internal method to tun the irrandiances calculation
        """
        # Get and store solar position arrays
        self._solar_position = self._site.get_solarposition(times=self._epw_hourly_data.index)
        if self.general_data['time_steps_per_hour'] > 1:
            m = str(60 / float(self.general_data['time_steps_per_hour'])) + 'min'
            self._solar_position = self._solar_position.resample(m).interpolate(
                method='ffill')  # Bfill for azimuth that is not continuous
        self.hourly_data["solar_position_apparent_zenith"] = self._solar_position['apparent_zenith'].values
        self.hourly_data["solar_position_zenith"] = self._solar_position['zenith'].values
        self.hourly_data["solar_position_apparent_elevation"] = self._solar_position['apparent_elevation'].values
        self.hourly_data["solar_position_elevation"] = self._solar_position['elevation'].values
        self.hourly_data["solar_position_azimuth"] = self._solar_position['azimuth'].values - 180.
        self.hourly_data["solar_position_equation_of_time"] = self._solar_position['equation_of_time'].values

        # Dataframe with hourly solar radiation per each direction
        azimuth_array = np.linspace(-180, 180, self.general_data['azimuth_subdivisions'] + 1, dtype=int)[:-1]
        height_array = np.linspace(90, 0, self.general_data['height_subdivisions'] + 1, dtype=int)[:-1]
        self.hourly_data_irradiances = {}
        for az in azimuth_array:
            self.hourly_data_irradiances[az] = {}
            for h in height_array:
                self.hourly_data_irradiances[az][h] = {}
                POA = _get_irradiance(self, h, az)
                self.hourly_data_irradiances[az][h]['global'] = POA['POA'].values
                self.hourly_data_irradiances[az][h]['direct'] = POA['POA_B'].values
                self.hourly_data_irradiances[az][h]['AOI'] = POA['AOI'].values
        # Horizontal
        POA = _get_irradiance(self, 0., 0.)
        self.hourly_data_irradiances[0][0] = {}
        self.hourly_data_irradiances[0][0]['global'] = POA['POA'].values
        self.hourly_data_irradiances[0][0]['direct'] = POA['POA_B'].values
        self.hourly_data_irradiances[0][0]['AOI'] = POA['AOI'].values

    @classmethod
    def from_pvgis(cls,
                   lat: float,
                   long: float,
                   country = None,
                   city = None,
                   year=None,
                   time_steps: int = 1,
                   irradiances_calculation: bool = True,
                   azimuth_subdivisions: int = 8,
                   height_subdivisions: int = 3,
                   urban_shading_tol=[80., 100., 80.]
        ):

        api_url = f'https://re.jrc.ec.europa.eu/api/v5_2/tmy?lat={lat:.4f}&lon={long:.4f}&outputformat=epw'

        response = requests.get(api_url).content.decode()
        loc = "unknown" if city == None else city
        ctr = "unknown" if country == None else country
        response = response.replace("LOCATION,unknown,-,unknown", f"LOCATION,{loc},-,{ctr}")
        # with open("test_epw.epw", 'w') as f:
        #     f.write("\n".join(response.splitlines()))

        w = cls(
            io.StringIO(response),
            year=None,
            time_steps= time_steps,
            irradiances_calculation = irradiances_calculation,
            azimuth_subdivisions = azimuth_subdivisions,
            height_subdivisions = height_subdivisions,
            urban_shading_tol = urban_shading_tol
        )

        return w




def _TskyCalc(T_ext, T_dp, P_, n_opaque, time_steps):
    '''Apparent sky temperature calculation procedure
    Martin Berdhal model used by TRNSYS

    Parameters
    ----------
    T_ext : numpy.array
        External Temperature [°C]
    T_dp : pandas.DataFrame
        External dew point temperature [°C]. It must be a pandas.DataFrame column
    P_ : pandas.DataFrame
        External pressure [Pa]. It must be a pandas.DataFrame column
    n_opaque : pandas.DataFrame
        Opaque sky covering [-]. It must be a pandas.DataFrame column

    Returns
    -------
    float
        Average temperature difference between External air temperature and Apparent sky temperature [°C]

    '''

    # Check input data type
    # Calculation Martin Berdhal model used by TRNSYS

    day = np.arange(24 * time_steps)  # Inizialization day vector
    T_sky = np.zeros(24 * time_steps)
    Tsky = []
    T_sky_year = []
    for i in range(365):
        for x in day:
            t = i * 24 + x
            Tdp = T_dp[t]
            P = P_[t] / 100  # [mbar]
            nopaque = n_opaque[t] * 0.1  # [0-1]
            eps_m = 0.711 + 0.56 * Tdp / 100 + 0.73 * (Tdp / 100) ** 2
            eps_h = 0.013 * np.cos(2 * np.pi * (x + 1) / 24)
            eps_e = 0.00012 * (P - 1000)
            eps_clear = eps_m + eps_h + eps_e  # Emissivity under clear sky condition
            C = nopaque * 0.9  # Infrared cloud amount
            eps_sky = eps_clear + (1 - eps_clear) * C  # Sky emissivity
            T_sky[x] = ((T_ext[t] + 273) * (eps_sky ** 0.25)) - 273  # Apparent sky temperature [°C]
        Tsky = np.append(Tsky, T_sky)  # Annual Tsky created day by day

    # Average temperature difference between External air temperature and Apparent sky temperature
    if time_steps > 1:
        dT_er = np.mean(T_ext - Tsky[:-time_steps + 1])
    else:
        dT_er = np.mean(T_ext - Tsky)

    return dT_er


def _get_irradiance(weather_obj, surf_tilt, surf_az):
    '''function from pvlib to calculate irradiance on a specific surface
    https://pvlib-python.readthedocs.io/en/stable/auto_examples/plot_ghi_transposition.html#sphx-glr-auto-examples-plot-ghi-transposition-py

    Parameters
    ----------
    weather_obj : eureca_building.weather.WeatherFile
        Weather file object
    surf_tilt: float
        tilt of the surface
    surf_az: float
        azimuth of the surface

    Returns
    -------
    pandas.DataFrame
        pandas DataFrame with the irradiances on the surface
    '''

    if surf_tilt < 0 or surf_tilt > 90 or surf_az < -180 or surf_az > 180:
        logging.warning(
            f"WARNING get_irradiance function, are you sure that the surface orientation is correct?? surf_tilt {surf_tilt}, surf_az {surf_az}")

    # Use pvlib function to calculate the irradiance on the surface

    surf_az = surf_az
    POA_irradiance = pvlib.irradiance.get_total_irradiance(
        surface_tilt=surf_tilt,
        surface_azimuth=surf_az,
        dni_extra=pvlib.irradiance.get_extra_radiation(weather_obj.hourly_data["time_index"], solar_constant=1366.1,
                                                       method='spencer'),
        dni=weather_obj._epw_hourly_data['dni'],
        ghi=weather_obj._epw_hourly_data['ghi'],
        dhi=weather_obj._epw_hourly_data['dhi'],
        solar_zenith=weather_obj.hourly_data["solar_position_apparent_zenith"],
        solar_azimuth=weather_obj.hourly_data["solar_position_azimuth"],
        model='isotropic',
        model_perez='allsitescomposite1990',
        airmass=weather_obj._site.get_airmass(solar_position=weather_obj._solar_position))
    AOI = pvlib.irradiance.aoi(
        surface_tilt=surf_tilt,
        surface_azimuth=surf_az,
        solar_zenith=weather_obj.hourly_data["solar_position_apparent_zenith"],
        solar_azimuth=weather_obj.hourly_data["solar_position_azimuth"])

    # Cleaning AOI vector

    for i in range(len(AOI)):
        if AOI[i] > 90 or weather_obj.hourly_data["solar_position_apparent_zenith"][i] > 90:
            AOI[i] = 90

    return pd.DataFrame({'GHI': weather_obj._epw_hourly_data['ghi'],
                         'POA': POA_irradiance['poa_global'],
                         'POA_B': POA_irradiance['poa_direct'],
                         'POA_D': POA_irradiance['poa_global'] - POA_irradiance['poa_direct'],
                         'AOI': AOI,
                         'solar zenith': weather_obj.hourly_data["solar_position_apparent_zenith"]})

if __name__ == "__main__":
    lat = 45.234
    long = 11.154
    country = "Italia"
    city = "Venezia"
    api_url = f'https://re.jrc.ec.europa.eu/api/v5_2/tmy?lat={lat:.4f}&lon={long:.4f}&outputformat=epw'

    response = requests.get(api_url).content.decode()
    loc = "unknown" if city == None else city
    ctr = "unknown" if country == None else country
    response = response.replace("LOCATION,unknown,-,unknown", f"LOCATION,{loc},-,{ctr}")
    # with open("test_epw.epw", 'w')as f:
    #     f.write("\n".join(response.splitlines()))

    epw = pvlib.iotools.parse_epw(io.StringIO(response))  # Reading the epw via pvlib

    w = WeatherFile.from_pvgis(
                   lat= 45.124,
                   long= 11.124,
                   country = "Italy",
                   city = "Venice",
        year=None,
        time_steps=1,
        irradiances_calculation=True,
        azimuth_subdivisions=10,
        height_subdivisions=4,
    )