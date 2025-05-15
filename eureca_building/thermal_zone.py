"""
This module includes functions to model the thermal zone
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import logging
import warnings

#import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
from scipy import interpolate

from eureca_building.config import CONFIG
from eureca_building.surface import Surface, SurfaceInternalMass
from eureca_building.fluids_properties import air_properties, vapour_properties
from eureca_building._VDI6007_auxiliary_functions import impedence_parallel, tri2star, long_wave_radiation, loadHK
from eureca_building.internal_load import Lights, ElectricLoad, People, InternalLoad
from eureca_building.ventilation import NaturalVentilation, Infiltration
from eureca_building.air_handling_unit import AirHandlingUnit, HeatRecoveryUnit
from eureca_building.domestic_hot_water import DomesticHotWater
from eureca_building.weather import WeatherFile
from eureca_building._auxiliary_function_for_monthly_calc import get_monthly_value_from_annual_vector
from eureca_building.setpoints import Setpoint
from eureca_building.exceptions import (
    Non3ComponentsVertex,
    SurfaceWrongNumberOfVertices,
    WindowToWallRatioOutsideBoundaries,
    InvalidSurfaceType,
    NonPlanarSurface,
    NegativeSurfaceArea,
)


# %% ThermalZone class


class ThermalZone(object):
    """Thermal zone class. Manages all the models and time step solution of the sensible and latent systems
    """

    def __init__(self, name: str, surface_list: list, net_floor_area=None, volume=None, number_of_units: int=1):
        """Int method. Creates the thermal zone object from a list of eureca_building.surface.Surface.
        Checks the inputs through properties

        Parameters
        ----------
        name : str
            Name of the zone
        surface_list : list
            list of eureca_building.surface.Surface or eureca_building.surface.SurfaceInternalMass objects
        net_floor_area : float, default None
            footprint area of the zone in m2. If None searches for a GroundFloor surface
        volume : float, default None
            volume of the zone in m3. If None sets 0 m3.
        number_of_units : int, default 1
            number of units/dwellings in the thermal zone (for DHW)
        NV_ACH_limit : float, default None
            maximum accepted ACH for the thermal zone in Vol/h. If None, leave None.
        """
        self.name = name
        self._surface_list = surface_list
        if volume == None:
            logging.warning(f"Thermal zone {self.name}, the volume is not set. Initialized with 0 m3")
            self._volume = 0.
        else:
            self._volume = volume
        if net_floor_area == None:
            floors_area = [surf._area for surf in self._surface_list if
                           surf._surface_type == 'GroundFloor']
            if len(floors_area) == 0:
                logging.warning(f"Thermal zone {self.name}, the footprint area is not set. Initialized with 0 m2")
                self._net_floor_area = 0.
            else:
                self._net_floor_area = np.array(floors_area).sum()
        else:
            self._net_floor_area = net_floor_area

        # Number of unit/dwelling in the thermal zone (for DHW)
        self.number_of_units = number_of_units

        self.internal_loads_list = []
        self.infiltration_list = []
        self.natural_ventilation = None
        self.domestic_hot_water_list = []
        self.design_heating_system_power = 1e20  # W
        self.design_cooling_system_power = -1e20  # W
        self.design_sensible_cooling_system_power = -1e20
        self.heating_sigma = {
            '1C': [0., 1.],
            '2C': [0., 0., 1.],
        } # Default convective load
        self.cooling_sigma = {
            '1C': [0., 1.],
            '2C': [0., 0., 1.],
        }  # Default convective load
        self.reset_init_values()

    @property
    def _surface_list(self) -> float:
        return self.__surface_list

    @_surface_list.setter
    def _surface_list(self, value: list):
        try:
            value = list(value)
        except ValueError:
            raise TypeError(f"Thermal zone {self.name}, the surface_list must be a list or a tuple: {type(value)}")
        for surface in value:
            if not isinstance(surface, Surface) and not isinstance(surface, SurfaceInternalMass):
                raise TypeError(f"Thermal zone {self.name}, non SUrface object in surface_list. ")
        self.__surface_list = value

    @property
    def _net_floor_area(self) -> float:
        return self.__net_floor_area

    @_net_floor_area.setter
    def _net_floor_area(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Thermal zone {self.name}, footprint area is not an float: {value}")
        if value < 0.0:
            logging.error(
                f"Thermal zone {self.name}, negative footprint area: {value}. Simulation will continue with absolute value"
            )
        if float(abs(value)) < 1e-5:
            self.__net_floor_area = 1e-5
        else:
            self.__net_floor_area = abs(value)

    @property
    def _volume(self) -> float:
        return self.__volume

    @_volume.setter
    def _volume(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Thermal zone {self.name}, volume is not an float: {value}")
        if value < 0.0:
            logging.error(
                f"Thermal zone {self.name}, negative volume: {value}. Simulation will continue with the absolute value"
            )
        if float(abs(value)) < 1e-5:
            self.__volume = 1e-5
        else:
            self.__volume = abs(value)
        self._air_thermal_capacity = self.__volume * air_properties["density"] * air_properties["specific_heat"]

    @property
    def number_of_units(self) -> int:
        return self._number_of_units

    @number_of_units.setter
    def number_of_units(self, value: int):
        try:
            value = int(value)
        except ValueError:
            raise TypeError(f"Thermal zone {self.name}, number of units must be an int: {value}")
        self._number_of_units = value

    @property
    def _air_thermal_capacity(self) -> float:
        return self.__air_thermal_capacity

    @_air_thermal_capacity.setter
    def _air_thermal_capacity(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Thermal zone {self.name}, air thermal capacity is not an float: {value}")
        if value < 0.0:
            logging.error(
                f"Thermal zone {self.name}, negative air thermal capacity: {value}. Simulation will continue with the absolute value"
            )
        if float(abs(value)) < 1e-5:
            self.__air_thermal_capacity = 1e-5
        else:
            self.__air_thermal_capacity = abs(value)

    def add_temperature_setpoint(self, setpoint, mode='air'):
        """Function to associate a setpoint objects to the thermal zone

        Parameters
        ----------
        setpoint : eureca_building.setpoint.Setpoint
            object of the class Setpoint
        mode : str
            setpoint mode: ['air', 'operative', 'radiant'].
            FOR NOW ONLY AIR IMPLEMENTED

        """
        self.temperature_setpoint = setpoint
        self.temperature_setpoint_mode = mode

    @property
    def temperature_setpoint(self):
        return self._temperature_setpoint

    @temperature_setpoint.setter
    def temperature_setpoint(self, value: Setpoint):
        if not isinstance(value, Setpoint):
            raise TypeError(
                f"Thermal zone {self.name}, setpoint must be a Setpoint object"
            )
        self._temperature_setpoint = value

    @property
    def temperature_setpoint_mode(self):
        return self._temperature_setpoint_mode

    @temperature_setpoint_mode.setter
    def temperature_setpoint_mode(self, value: str):
        if not isinstance(value, str):
            raise TypeError(
                f"Thermal zone {self.name}, setpoint mode must be a str"
            )
        if value not in ['air', 'operative', 'radiant']:
            raise ValueError(
                f"Thermal zone {self.name}, setpoint mode must be chosen from 'air', 'operative', 'radiant'"
            )
        self._temperature_setpoint_mode = value

    def add_humidity_setpoint(self, setpoint):
        """Function to associate a humidity setpoint to the thermal zone

        Parameters
        ----------
        setpoint : eureca_building.setpoint.Setpoint
            object of the class Setpoint

        """
        self.humidity_setpoint = setpoint
        # the mode is in the SP object
        self.humidity_setpoint_mode = setpoint.setpoint_type

    @property
    def humidity_setpoint(self):
        return self._humidity_setpoint

    @humidity_setpoint.setter
    def humidity_setpoint(self, value: Setpoint):
        if not isinstance(value, Setpoint):
            raise TypeError(
                f"Thermal zone {self.name}, setpoint must be a Setpoint object"
            )
        self._humidity_setpoint = value

    @property
    def humidity_setpoint_mode(self):
        return self._humidity_setpoint_mode

    @humidity_setpoint_mode.setter
    def humidity_setpoint_mode(self, value: str):
        if not isinstance(value, str):
            raise TypeError(
                f"Thermal zone {self.name}, setpoint mode must be a str"
            )
        if value not in ['relative_humidity']:
            raise ValueError(
                f"Thermal zone {self.name}, setpoint mode must be chosen from 'relative_humidity'"
            )
        self._humidity_setpoint_mode = value

    @staticmethod
    def get_specific_humidity(air_t, air_rh, p_atm):
        """Static method. Give the specific humidity from temperature, relative humidity, and atmospheric pressure
        
        Parameters
        ----------
        air_t : float
            air temperature [°C]
        air_rh : float
            air relative humidity [-]
        p_atm : float
            atmospheric pressure [pa]

        Returns
        ----------
        float
            specific humidity [kg_vap/kg_da]

        """
        if air_rh < 0.:
            logging.warning("A relative humidity is below zero. Simulation will continue with a minimum of zero")
            air_rh = 0.
        if air_rh > 1.:
            logging.warning("A relative humidity is above one. Simulation will continue with a maximum of one")
            air_rh = 1.
        if air_t < 0:
            p_intsat = 610.5 * np.exp((21.875 * air_t) / (265.5 + air_t))
        else:
            p_intsat = 610.5 * np.exp((17.269 * air_t) / (237.3 + air_t))

        x_int = 0.622 * (air_rh * p_intsat / (p_atm - (air_rh * p_intsat)))
        return x_int, p_intsat

    def add_internal_load(self, *internal_load):
        """Function to associate a load to the thermal zone

        Parameters
        ----------
        internal_load : eureca_building.internal_load.InternalLoad
            As many eureca_building.internal_load.InternalLoad can be provided

        """
        for int_load in internal_load:
            if not isinstance(int_load, InternalLoad):
                raise TypeError(
                    f"ThermalZone {self.name}, add_internal_load() method: internal_load not of InternalLoad type: {type(int_load)}"
                )
            self.internal_loads_list.append(int_load)

    def extract_convective_radiative_latent_electric_load(self):
        """
        Extracts internal load components from the thermal zone.
        
        Returns
        -------
        dict
            Dictionary with keys: 'convective [W]', 'radiative [W]', 'latent [kg_vap/s]', 'electric [W]'
            Each mapped to a NumPy array of length equal to CONFIG.number_of_time_steps_year.
        """
        convective = np.zeros(CONFIG.number_of_time_steps_year)
        radiant = np.zeros(CONFIG.number_of_time_steps_year)
        latent = np.zeros(CONFIG.number_of_time_steps_year)
        electric = np.zeros(CONFIG.number_of_time_steps_year)
        for load in self.internal_loads_list:
            load_conv, load_rad, load_lat, load_el = load.get_loads(area=self._net_floor_area)
            convective += load_conv
            radiant += load_rad
            latent += load_lat
            electric += load_el
        self.convective_load = convective
        self.radiative_load = radiant
        self.latent_load = latent
        self.electric_load = electric
        return {'convective [W]': convective,
                'radiative [W]': radiant,
                'latent [kg_vap/s]': latent,
                'electric [W]': electric}

    def add_infiltration(self, *infiltration):
        """Function to associate a natural ventilation object to the thermal zone

        Parameters
        ----------
        natural_ventilation : eureca_building.ventilation.Infiltration
            As many eureca_building.ventilation.Infiltration can be provided

        """
        for inf in infiltration:
            if not isinstance(inf, Infiltration):
                raise TypeError(
                    f"ThermalZone {self.name}, add_infiltration() method: infiltrion not of Infiltration type: {type(inf)}"
                )
            self.infiltration_list.append(inf)

    def calc_infiltration(self, weather):
        """From the infiltration_list calculates 2 arrays (len equal to 8769 * number of time steps per hour):
        {
        air mass flow rate [kg/s] : np.array
        vapour mass flow rate [kg/s] : np.array
        }

        Parameters
        ----------
        weather : eureca_building.weather.WeatherFile
            WeatherFile object

        Returns
        ----------
        dict
        """
        infiltration_air_flow_rate = np.zeros(CONFIG.number_of_time_steps_year)
        infiltration_vapour_flow_rate = np.zeros(CONFIG.number_of_time_steps_year)
        for inf in self.infiltration_list:
            air_rate, vapour_rate = inf.get_flow_rate(weather, area=self._net_floor_area, volume=self._volume)
            infiltration_air_flow_rate += air_rate
            infiltration_vapour_flow_rate += vapour_rate
        self.infiltration_air_flow_rate = infiltration_air_flow_rate
        self.infiltration_vapour_flow_rate = infiltration_vapour_flow_rate
        self.nat_vent_air_flow_rate = np.zeros(CONFIG.number_of_time_steps_year)
        self.nat_vent_info = {
            'airflow_rate' : {'kg/s' : np.zeros([CONFIG.number_of_time_steps]),
                              'm3/s' : np.zeros([CONFIG.number_of_time_steps]),
                              'L/s' : np.zeros([CONFIG.number_of_time_steps]),
                              'm3/h' : np.zeros([CONFIG.number_of_time_steps]),
                              'vol/h' : np.zeros([CONFIG.number_of_time_steps]),
                              },
            'windows_opening' : {'open_fraction' : np.zeros([CONFIG.number_of_time_steps]),
                                 },
            }
        return {'infiltration_air_flow_rate [kg/s]': infiltration_air_flow_rate,
                'infiltration_vapour_flow_rate [kg/s]': infiltration_vapour_flow_rate, }
    
    def add_natural_ventilation(self, natural_ventilation, weather):
        """Function to associate a natural ventilation object to the thermal zone

        Parameters
        ----------
        natural_ventilation : eureca_building.ventilation.NaturalVentilation
            As many eureca_building.ventilation.NaturalVentilation can be provided

        """
        if not isinstance(natural_ventilation, NaturalVentilation):
            raise TypeError(
                f"ThermalZone {self.name}, add_infiltration() method: infiltrion not of Infiltration type: {type(natural_ventilation)}"
            )
        self.natural_ventilation = natural_ventilation
        ext_wall_list = [s for s in self._surface_list if s.surface_type == "ExtWall"]
        self.natural_ventilation.define_pressure_coef(
            weather = weather,
            surfaces_with_opening = ext_wall_list,
        )
        

    def add_air_handling_unit(self, ahu, weather):
        """Function to associate an air_handling_unit object to the thermal zone

        Parameters
        ----------
        ahu : eureca_building.air_handling_unit.AirHandlingUnit
            AirHandlingUnit object
        weather : eureca_building.weather.WeatherFile
            WeatherFile object
        """
        self.air_handling_unit = ahu

    @property
    def air_handling_unit(self) -> AirHandlingUnit:
        return self._air_handling_unit

    @air_handling_unit.setter
    def air_handling_unit(self, value: AirHandlingUnit):
        if not isinstance(value, (AirHandlingUnit, HeatRecoveryUnit)):
            raise TypeError(
                f"Thermal zone {self.name}, air_handling_unit must be a AirHandlingUnit object"
            )
        self._air_handling_unit = value

        # mechanical_ventilation_air_flow_rate = np.zeros(CONFIG.number_of_time_steps_year)
        # mechanical_ventilation_vapour_flow_rate = np.zeros(CONFIG.number_of_time_steps_year)
        # air_rate, vapour_rate = self._air_handling_unit.mechanical_ventilation.get_flow_rate(weather, area=self._net_floor_area, volume=self._volume)
        # mechanical_ventilation_air_flow_rate = air_rate
        # mechanical_ventilation_vapour_flow_rate = vapour_rate
        # self.mechanical_ventilation_air_flow_rate = mechanical_ventilation_air_flow_rate
        # self.mechanical_ventilation_vapour_flow_rate = mechanical_ventilation_vapour_flow_rate

    # def _plot_Zone_Natural_Ventilation(self, weather_file):
    #     fig, [ax1, ax2] = plt.subplots(nrows=2)
    #     ax1_ = ax1.twinx()
    #     pd.DataFrame({
    #         'infiltration_air_flow_rate [kg/s]': self.infiltration_air_flow_rate,
    #     }).plot(ax=ax1)
    #     pd.DataFrame({
    #         'infiltration_vapour_flow_rate [kg/s]': self.infiltration_vapour_flow_rate,
    #     }).plot(ax=ax1_, color='r')
    #     pd.DataFrame({
    #         'oa_specific_humidity': weather_file.hourly_data['out_air_specific_humidity'],
    #         'theta_ext': weather_file.hourly_data['out_air_db_temperature'],
    #         'oa_relative_humidity': weather_file.hourly_data['out_air_relative_humidity']
    #     }).plot(ax=ax2)

    def add_domestic_hot_water(self, weather_obj, *dhw_obj):
        """
        Function to associate a domestic hot water object to the thermal zone

        Parameters
        ----------
        dhw_obj : eureca_building.domestic_hot_water.DomesticHotWater
            DomesticHotWater object
        weather : eureca_building.weather.WeatherFile
            WeatherFile object
        """
        self.domestic_hot_water_volume_flow_rate = 0 # m3/s
        self.domestic_hot_water_demand = 0 # W
        for dhw in dhw_obj:
            if not isinstance(dhw, DomesticHotWater):
                raise TypeError(
                    f"Thermal zone {self.name}: add_domestic_hot_water input is not a DomesticHotWater object"
                )
            self.domestic_hot_water_list.append(dhw)
            volume, demand = dhw.get_dhw_yearly_mass_flow_rate(self._net_floor_area, self.number_of_units, weather_obj)
            self.domestic_hot_water_volume_flow_rate += volume  # m3/s
            self.domestic_hot_water_demand += demand  # W

    def _ISO13790_params(self):
        '''Calculates the thermal zone parameters of the ISO 13790
        it does not require input

        Htr_is
        Htr_w
        Htr_ms
        Htr_emCm
        DenAm
        Atot
        Htr_op
        UA_tot
        '''

        self.Htr_is = 0.
        self.Htr_w = 0.
        self.Htr_ms = 0.
        self.Htr_em = 0.
        self.Cm = 0.
        self.DenAm = 0.
        self.Atot = 0.
        self.Htr_op = 0.
        self.ext_wall_opaque_area = 0.
        self.ext_wall_glazed_area = 0.
        self.ext_roof_area = 0.

        # list all surface to extract the window and opeque area and other thermo physical prop

        for surface in self._surface_list:
            try:
                if surface.surface_type == 'IntFloor':
                    self.Cm += surface._opaque_area * surface.construction.k_est
                    self.DenAm += surface._opaque_area * surface.construction.k_est ** 2
                else:
                    self.Cm += surface._opaque_area * surface.construction.k_int
                    self.DenAm += surface._opaque_area * surface.construction.k_int ** 2
                self.Atot += surface._area

                if surface.surface_type in ["ExtWall", "GroundFloor", "Roof"]:
                    self.Htr_op += surface._opaque_area * surface.construction._u_value
                    if surface._glazed_area > 0.:
                        self.Htr_w += surface._glazed_area * surface.window._u_value
                        
                if surface.surface_type == "Roof":
                    self.ext_roof_area += surface._area

                if surface.surface_type == "ExtWall":
                    self.ext_wall_opaque_area += surface._opaque_area
                    self.ext_wall_glazed_area += surface._glazed_area

            except AttributeError:
                raise AttributeError(
                    f"Thermal zone {self.name}, surface {surface.name} construction or window not specified"
                )

        # Final calculation
        self.htr_ms = 9.1  # heat tranfer coeff. ISO 13790 [W/(m2 K)]
        self.h_is = 3.45  # heat tranfer coeff. ISO 13790 [W/(m2 K)]

        self.Am = self.Cm ** 2 / self.DenAm
        self.Htr_ms = self.Am * self.htr_ms
        self.Htr_em = 1 / (1 / self.Htr_op - 1 / self.Htr_ms)
        self.Htr_is = self.h_is * self.Atot
        self.UA_tot = self.Htr_op + self.Htr_w

    def print_ISO13790_params(self):
        """Just a beuty print of ISO 13790 parameters

        Returns
        -------
        str

        """
        return f"""
Thermal zone {self.name} 1C params:        
        Cm: {self.Cm:.5f}
        Htr_is: {self.Htr_is:.5f}
        Htr_w: {self.Htr_w:.5f}
        Htr_ms: {self.Htr_ms:.5f}
        Htr_em: {self.Htr_em:.5f}
        Htr_op: {self.Htr_op:.5f}
        """

    def _VDI6007_params(self):
        '''Calculates the thermal zone parameters of the VDI 6007
        it does not require input

        Araum_tot
        Aaw_tot
        Araum_opaque
        Aaw_opaque
        R1AW, R1IW, C1AW, C1IW
        RgesAW
        RrestAW
        RalphaStarIL, RalphaStarAW, RalphaStarIW
        UA_tot
        Htr_op
        Htr_w
        '''

        # Creation of some arrys of variables and parameters

        alphaStr = 5  # vdi Value
        alphaKonA = 20  # vdi Value
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
        RalphaStrAW = np.array([])
        RalphaStrIW = np.array([])
        RalphaStrAF = np.array([])
        AreaAW = np.array([])
        AreaAF = np.array([])
        AreaIW = np.array([])
        self.Araum_tot = 0
        self.Aaw_tot = 0
        self.Araum_opaque = 0
        self.Aaw_opaque = 0
        self.ext_wall_opaque_area = 0.
        self.ext_wall_glazed_area = 0.
        self.ext_roof_area = 0.

        # Cycling surface to calculates the Resistance and capacitance of the vdi 6007

        for surface in self._surface_list:
            self.Araum_tot += surface._area
            self.Araum_opaque += surface._opaque_area
            if surface.surface_type == "ExtWall":
                self.ext_wall_opaque_area += surface._opaque_area
                self.ext_wall_glazed_area += surface._glazed_area
            if surface.surface_type == "Roof":
                self.ext_roof_area += surface._area
            if surface.surface_type in ["ExtWall", "GroundFloor", "Roof"]:
                self.Aaw_tot += surface._area
                self.Aaw_opaque += surface._opaque_area
                surface_R1, surface_C1 = surface.get_VDI6007_surface_params(asim=True)
                C1AW_v = np.append(C1AW_v, [surface_C1], axis=0)
                # Opaque params
                HAW_v = np.append(HAW_v, surface.construction._u_value * surface._opaque_area)
                alphaKonAW = np.append(alphaKonAW,
                                       [surface._opaque_area * (surface.construction._conv_heat_trans_coef_int)],
                                       axis=0)
                RalphaStrAW = np.append(RalphaStrAW, [1 / (surface._opaque_area * surface.construction.rad_heat_trans_coef)])

                try:
                    # Glazed params
                    R_AF_v = (surface.window.Rl_w / surface._glazed_area) if abs(surface._glazed_area) > 1e-5 else 1e15
                    # Eq 26
                    HAF_v = np.append(HAF_v, surface.window._u_value * surface._glazed_area)
                    alphaKonAF = np.append(alphaKonAF,
                                           [surface._glazed_area * (1 / surface.window.Ri_w - surface.construction.rad_heat_trans_coef)], axis=0)
                    RalphaStrAF = np.append(RalphaStrAF, [1 / (surface._glazed_area * surface.construction.rad_heat_trans_coef)])
                except (AttributeError, ZeroDivisionError):
                    # case of no glazed area
                    R_AF_v = 1e15
                    HAF_v = np.append(HAF_v, 0.)
                    alphaKonAF = np.append(alphaKonAF,
                                           [0.], axis=0)
                    RalphaStrAF = np.append(RalphaStrAF, [1e15])
                # this part is a little different in Jacopo model,
                # However this part calculates opaque R, glazed R and insert the parallel as wall R
                # R1AW_v = np.append(R1AW_v, [1 / (1 / surface_R1 + 1 / R_AF_v)], axis=0)
                R1AW_v = np.append(R1AW_v,[1/(1/surface_R1+6/R_AF_v)], axis=0) # ALTERNATIVA NORMA

                AreaAW = np.append(AreaAW, surface._opaque_area)
                AreaAF = np.append(AreaAF, surface._glazed_area)

            elif surface.surface_type in ["IntWall", "IntCeiling", "IntFloor"]:
                surface_R1, surface_C1 = surface.get_VDI6007_surface_params(asim=False)
                R1IW_m = np.append(R1IW_m, [surface_R1], axis=0)
                C1IW_m = np.append(C1IW_m, [surface_C1], axis=0)
                R_IW = np.append(R_IW, [sum(surface.construction.thermal_resistances)], axis=0)
                alphaKonIW = np.append(alphaKonIW,
                                       [surface._opaque_area * (surface.construction._conv_heat_trans_coef_int)],
                                       axis=0)

                # if surface.opaqueArea*alphaStr == 0:
                #    print(surface.name)
                RalphaStrIW = np.append(RalphaStrIW, [1 / (surface._opaque_area * surface.construction.rad_heat_trans_coef)])

                AreaIW = np.append(AreaIW, surface._area)
            else:
                raise TypeError(f'Surface {surface.name}: surface type not found: {surface.surface_type}.')

        # Doing the parallel of the impedances

        self.R1AW, self.C1AW = impedence_parallel(R1AW_v, C1AW_v)  # eq 22
        self.R1IW, self.C1IW = impedence_parallel(R1IW_m, C1IW_m)

        # Final params

        self.RgesAW = 1 / (sum(HAW_v) + sum(HAF_v))  # eq 27

        RalphaKonAW = 1 / (sum(alphaKonAW) + sum(alphaKonAF))  # scalar
        RalphaKonIW = 1 / sum(alphaKonIW)  # scalar

        if sum(AreaAW) <= sum(AreaIW):
            RalphaStrAWIW = 1 / (sum(1 / RalphaStrAW) + sum(1 / RalphaStrAF))  # eq 29
        else:
            RalphaStrAWIW = 1 / sum(1 / RalphaStrIW)  # eq 31

        self.RrestAW = self.RgesAW - self.R1AW - 1 / (1 / RalphaKonAW + 1 / RalphaStrAWIW)  # eq 28

        RalphaGesAW_A = 1 / (alphaKonA * (sum(AreaAF) + sum(AreaAW)))

        if self.RgesAW < RalphaGesAW_A:  # this is different from Jacopo's model but equal to the standard
            self.RrestAW = RalphaGesAW_A  # eq 28a
            self.R1AW = self.RgesAW - self.RrestAW - 1 / (1 / RalphaKonAW + 1 / RalphaStrAWIW)  # eq 28b

            if self.R1AW < 10 ** (-10):
                self.R1AW = 10 ** (-10)  # Thresold (only numerical to avoid division by zero)  #eq 28c

        self.RalphaStarIL, self.RalphaStarAW, self.RalphaStarIW = tri2star(RalphaStrAWIW, RalphaKonIW, RalphaKonAW)
        self.UA_tot = sum(HAW_v) + sum(HAF_v)
        self.Htr_op = sum(HAW_v)
        self.Htr_w = sum(HAF_v)

    def print_VDI6007_params(self):
        """Just a beuty print of VDI6007 parameters

        Returns
        -------
        str

        """
        return f"""
Thermal zone {self.name} 2C params:        
        R1AW: {self.R1AW:.10f}
        R1IW: {self.R1IW:.10f}   
        C1AW: {self.C1AW:.1f}
        C1IW: {self.C1IW:.1f}
        RrestAW: {self.RrestAW:.10f}
        RalphaStarIL: {self.RalphaStarIL:.10f}
        RalphaStarAW: {self.RalphaStarAW:.10f}
        RalphaStarIW: {self.RalphaStarIW:.10f}
        """

    def calculate_zone_loads_ISO13790(self, weather):
        '''Calculates the heat gains (internal and solar) on the three nodes of the ISO 13790 network
        Vectorial calculation

        Parameters
        ----------
        weather : eureca_building.weather.WeatherFile
            WeatherFile object

        '''

        # Check input data type

        if not isinstance(weather, WeatherFile):
            raise TypeError(f'ThermalZone {self.name}, weather type is not a WeatherFile: weather type {type(weather)}')
        # Check input data quality
        # First calculation of internal heat gains

        phi_int = self.extract_convective_radiative_latent_electric_load()
        phi_sol_gl_tot = 0
        phi_sol_op_tot = 0

        # Solar radiation
        irradiances = weather.hourly_data_irradiances

        for surface in self._surface_list:

            phi_sol_op = 0
            phi_sol_gl = 0

            if surface.surface_type in ['ExtWall', 'Roof']:
                F_sh_urban_shading = surface.shading_coefficient if hasattr(surface, "shading_coefficient") else 1.
                if hasattr(surface, 'window'):
                    h_r = surface.get_surface_external_radiative_coefficient()

                    irradiance = irradiances[float(surface._azimuth_round)][float(surface._height_round)]
                    
                    
                    BRV = irradiance['direct']
                    TRV = irradiance['global']
                    DRV = TRV - BRV
                    A_ww = surface._glazed_area
                    A_op = surface._opaque_area

                    # Glazed surfaces
                    F_sh_w = surface.window._shading_coef_ext
                    F_w = surface.window._shading_coef_int
                    F_f = surface.window._frame_factor
                    AOI = irradiance['AOI']
                    shgc = interpolate.splev(AOI, surface.window.solar_heat_gain_coef_profile, der=0)
                    shgc_diffuse = interpolate.splev(70, surface.window.solar_heat_gain_coef_profile, der=0)
                    # if i.OnOff_shading == 'On':
                    #     phi_sol_gl = F_so * (BRV * F_sh * F_w * (
                    #             1 - F_f) * shgc * A_ww * i.shading_effect) + F_so * DRV * F_sh * F_w * (
                    #                          1 - F_f) * shgc_diffuse * A_ww
                    # else:
                    phi_sol_gl = F_sh_urban_shading * BRV * F_sh_w * F_w * (1 - F_f) * shgc * A_ww \
                                 + DRV * F_sh_w * F_w * (1 - F_f) * shgc_diffuse * A_ww

                # Opaque surfaces
                F_r = surface._sky_view_factor
                alpha = surface._construction.ext_absorptance
                sr = surface._construction._R_se
                U_net = surface._construction._u_value_net
                # if i.OnOff_shading == 'On':
                #     phi_sol_op = F_sh * (
                #             BRV * i.shading_effect + DRV) * self.alpha * self.sr_ew * self.U_ew_net * self.A_ew - self.F_r * self.sr_ew * self.U_ew_net * self.A_ew * h_r * weather.dT_er
                # else:
                phi_sol_op = F_sh_urban_shading * TRV * alpha * sr * U_net * A_op - \
                             F_r * sr * U_net * A_op * h_r * \
                             weather.general_data['average_dt_air_sky']

            # Total solar gain
            phi_sol_gl_tot += phi_sol_gl
            phi_sol_op_tot += phi_sol_op

        phi_sol = phi_sol_gl_tot + phi_sol_op_tot

        # Distribute heat gains to temperature nodes
        self.phi_ia = phi_int['convective [W]']
        self.phi_st = (1 - self.Am / self.Atot - self.Htr_w / (9.1 * self.Atot)) * (phi_int['radiative [W]'] + phi_sol)
        self.phi_m = self.Am / self.Atot * (phi_int['radiative [W]'] + phi_sol)

        # For quasi-steadystate
        int_gain = phi_int['convective [W]'] + phi_int['radiative [W]']
        self.monthly_int_gain = get_monthly_value_from_annual_vector(int_gain, method = 'integral') # [J]
        self.monthly_sol_gain = get_monthly_value_from_annual_vector(phi_sol, method = 'integral') # [J]

    def calculate_zone_loads_VDI6007(self, weather):
        '''Calculates zone loads for the vdi 6007 standard
        Also the external equivalent temperature

        Parameters
        ----------
        weather : eureca_building.weather.WeatherFile
            WeatherFile obj

        '''

        # Check input data type

        if not isinstance(weather, WeatherFile):
            raise TypeError(f'ThermalZone {self.name}, weather type is not a WeatherFile: weather type {type(weather)}')


        # Eerd = Solar_gain['0.0']['0.0']['Global']*rho_ground                          #
        # Eatm = Solar_gain['0.0']['0.0']['Global']-Solar_gain['0.0']['0.0']['Direct']
        #
        # T_ext vettore


        Eatm, Eerd, theta_erd, theta_atm = long_wave_radiation(weather.hourly_data['out_air_db_temperature'])

        alpha_str_lw = (Eatm + Eerd)/(theta_atm - theta_erd)
        alpha_str_lw[(alpha_str_lw == 0.) & ((theta_atm - theta_erd) == 0)] = 5.

        # Creates some vectors and set some parameters

        T_ext = weather.hourly_data['out_air_db_temperature']
        theta_eq = np.zeros([len(T_ext), len(self._surface_list)])
        delta_theta_eq_lw = np.zeros([len(T_ext), len(self._surface_list)])
        delta_theta_eq_kw = np.zeros([len(T_ext), len(self._surface_list)])
        theta_eq_w = np.zeros([len(T_ext), len(self._surface_list)])
        Q_il_str_A_iw = 0
        Q_il_str_A_aw = 0
        Q_il_kon_A = 0

        # Solar radiation
        irradiances = weather.hourly_data_irradiances

        i = -1

        # Lists all surfaces to calculate the irradiance on each one and creates the solar gains

        for surface in self._surface_list:
            i += 1
            if surface.surface_type in ['ExtWall', 'Roof']:
                # Some value loaded from surface and weather object
                h_r = surface.get_surface_external_radiative_coefficient()
                alpha_str_a = surface.construction.rad_heat_trans_coef
                alpha_a = surface.construction._conv_heat_trans_coef_ext + surface.construction.rad_heat_trans_coef
                phi = surface._sky_view_factor
                F_sh_urban_shading = surface.shading_coefficient if hasattr(surface, "shading_coefficient") else 1.
                irradiance = irradiances[float(surface._azimuth_round)][float(surface._height_round)]
                AOI = irradiance['AOI']
                BRV = irradiance['direct']
                TRV = irradiance['global']
                # Delta T long wave
                delta_theta_eq_lw[:, i] = ((theta_erd - T_ext) * (1 - phi) +
                                           (theta_atm - T_ext) * phi) * alpha_str_lw * 0.9 / (0.93 * alpha_a)
                delta_theta_eq_kw[:, i] = (BRV * F_sh_urban_shading + (TRV - BRV)) * surface.construction.ext_absorptance / alpha_a
                theta_eq[:, i] = (T_ext + delta_theta_eq_lw[:, i] + delta_theta_eq_kw[:, i]) * \
                                 surface.construction._u_value * surface._opaque_area / self.UA_tot
                if hasattr(surface, 'window'):
                    theta_eq_w[:, i] = (T_ext + delta_theta_eq_lw[:, i]) * \
                                       surface.window.u_value * surface._glazed_area / self.UA_tot
                    frame_factor = 1 - surface.window._frame_factor
                    F_sh = surface.window._shading_coef_ext
                    F_w = surface.window._shading_coef_int
                    F_sh_w = F_sh * F_w * F_sh_urban_shading

                    shgc_rad = interpolate.splev(AOI, surface.window.solar_heat_gain_coef_profile_rad, der=0)
                    shgc_diffuse_rad = interpolate.splev(70, surface.window.solar_heat_gain_coef_profile_rad, der=0)
                    shgc_conv = interpolate.splev(AOI, surface.window.solar_heat_gain_coef_profile_conv, der=0)
                    shgc_diffuse_conv = interpolate.splev(70, surface.window.solar_heat_gain_coef_profile_conv, der=0)

                    radiative_inward_solar_radiation = frame_factor * surface._glazed_area * ( shgc_rad * BRV * F_sh_w +
                                                                                        shgc_diffuse_rad * (TRV - BRV))
                    convective_inward_solar_radiation = frame_factor * surface._glazed_area * (shgc_conv * BRV * F_sh_w +
                                                                                              shgc_diffuse_conv * (
                                                                                                          TRV - BRV))

                    # Jacopo quì usa come A_v l'area finestrata, mentre la norma parla di area finestrata + opaca per la direzione
                    Q_il_str_A_iw += radiative_inward_solar_radiation * (
                                             (self.Araum_tot - self.Aaw_tot) / (self.Araum_tot - surface._area))
                    Q_il_str_A_aw += radiative_inward_solar_radiation * (
                                             (self.Aaw_tot - surface._area) / (self.Araum_tot - surface._area))
                    Q_il_kon_A += convective_inward_solar_radiation

            if surface.surface_type == 'GroundFloor':
                theta_eq[:, i] = T_ext * surface.construction._u_value * surface._opaque_area / self.UA_tot

        # self.Q_il_str_A = self.Q_il_str_A.to_numpy()
        # self.carichi_sol = (Q_il_str_A_iw+ Q_il_str_A_aw).to_numpy()

        self.theta_eq_tot = theta_eq.sum(axis=1) + theta_eq_w.sum(axis=1)

        # Calculates internal heat gains
        phi_int = self.extract_convective_radiative_latent_electric_load()

        Q_il_str_I = phi_int['radiative [W]']
        self.Q_il_kon_I = phi_int['convective [W]'] + Q_il_kon_A

        Q_il_str_I_iw = Q_il_str_I * (self.Araum_tot - self.Aaw_tot) / self.Araum_tot
        Q_il_str_I_aw = Q_il_str_I * self.Aaw_tot / self.Araum_tot

        self.Q_il_str_iw = Q_il_str_A_iw + Q_il_str_I_iw
        self.Q_il_str_aw = Q_il_str_A_aw + Q_il_str_I_aw


        # sigma_fhk: fraction of radiant heating/cooling surfaces on total heating/cooling load
        # sigma_fhk_aw: fraction of radiant heating/cooling inside external walls on the total radiant heating/cooling load
        # sigma_hk_str: this is the radiant fraction the heating/cooling systems that are not embedded in the surfaces
        #
        # Let's say that a roof has 3 systems:
        #     1) a radiant internal ceiling (20% of the total load)
        #     2) a radiant floor (30% of the total load) in contact to the ground
        #     3) a radiator working 80% convective and 20% radiant (50% of the total load)
        #
        # sigma_fhk: 0.5
        #     sum of 20% for the radiant ceiling and 30% for the radiant floor
        # sigma_fhk_aw: 0.6
        #     because 0.3/(0.2+0.3) is equal to 0.6, i.e. the part of the radiant surface load
        #     associated to external surfaces
        # sigma_hk_str: 0.2
        #     because the radiator is radiative for the 20%

        # self.sigma = loadHK(sigma_fhk, sigma_fhk_aw, sigma_hk_str, self.Aaw, self.Araum)

    # def _plot_ISO13790_IHG(self):
    #     fig, ax = plt.subplots()
    #     pd.DataFrame({
    #         'phi_ia': self.phi_ia,
    #         'phi_st': self.phi_st,
    #         'phi_m': self.phi_m
    #     }).plot(ax=ax)

    # def _plot_VDI6007_IHG(self, weather_file):
    #     fig, [ax1, ax2] = plt.subplots(nrows=2)
    #     pd.DataFrame({
    #         'Q_il_kon_I': self.Q_il_kon_I,
    #         'Q_il_str_iw': self.Q_il_str_iw,
    #         'Q_il_str_aw': self.Q_il_str_aw
    #     }).plot(ax=ax1)
    #     pd.DataFrame({
    #         'theta_eq_tot': self.theta_eq_tot,
    #         'theta_ext': weather_file.hourly_data['out_air_db_temperature']
    #     }).plot(ax=ax2)

    def sensible_balance_1C(self, flag, Hve, T_e, T_sup_AHU, phi_load, sigma=[0., 1.], T_set=20., phi_HC_set=0.):
        '''Solves ISO 13790 network for a specific time step

        Parameters
        ----------
        flag : string
            flag to set the calculation method string 'Tset' or 'phiset'
        Hve : list
            list of two ositive floats. Ventilation and infiltration heat tranfer coeff [W/K]
            [Hve_vent, Hve_inf]
        T_e : float
            timestep external temperature [°C]
        T_sup_AHU : float
            timestep ventilation supply temperature [°C]
        phi_load : list
            list of three floats. The load on the three nodes network (ia, sm, m respectively) [W]
            [phi_ia, phi_sm, phi_m]
        sigma : list, default [0., 1.]
            list of 2 floats
            portion of the heating/cooling load to:
            sigma[0]: radiant to surface
            sigma[1]: convective to air node
            sum must be 1
        T_set : float, default 20.
            time step setpoint temperature [°C]
        phi_HC_set : float, default 0.
            time step thermal load to the ambient [W]

        Returns
        -------
        numpy.array
            array with system load, air temperature, surfaces equivalent temperature and mass equivalent temperature
            e.g.
            [
            demand [W],
            T_air [°C],
            T_s [°C],
            T_m [°C]
            ]
        '''

        # Check input data type

        if flag != 'Tset' and flag != 'phiset':
            raise TypeError(
                f"ERROR Thermal zone {self.name}, Sensible1C flag input is not a 'Tset' or 'phiset': flag {flag}")
        # if np.abs((np.array(sigma).sum() - 1)) > 1e-3:
        #     raise ValueError(
        #         f"Thermal Zone {self.name}, Sensible1C: sigma total must be 1. Sigma = {sigma}"
        #     )
        # Set some data and build up the system

        phi_ia = phi_load[0]
        phi_st = phi_load[1]
        phi_m = phi_load[2]
        Hve_vent = Hve[0]
        Hve_inf = Hve[1]
        tau = CONFIG.time_step
        rad_factor = sigma[0]
        conv_factor = sigma[1]
        Y = np.zeros((3, 3))
        q = np.zeros(3)

        if flag == 'Tset':
            Y[0, 0] = conv_factor
            Y[0, 1] = self.Htr_is
            Y[1, 0] = rad_factor
            Y[1, 1] = -(self.Htr_is + self.Htr_w + self.Htr_ms)
            Y[1, 2] = self.Htr_ms
            Y[2, 1] = self.Htr_ms
            Y[2, 2] = -self.Cm / tau - self.Htr_em - self.Htr_ms

            q[0] = Hve_inf * (T_set - T_e) + Hve_vent * (T_set - T_sup_AHU) - phi_ia + self.Htr_is * T_set \
                   + self._air_thermal_capacity * (T_set - self.Ta0) / tau
            q[1] = -self.Htr_is * T_set - phi_st - self.Htr_w * T_e
            q[2] = -self.Htr_em * T_e - phi_m - self.Cm * self.Tm0[0] / tau
            x = np.linalg.inv(Y).dot(q)
            # Seems to be more computationally efficient then np.insert
            return np.array([num for num in x[:1]] + [T_set] + [num for num in x[1:]]) # np.insert(x, 1, T_set)

        elif flag == 'phiset':
            Y[0, 0] = -(self.Htr_is + Hve_inf + Hve_vent) - self._air_thermal_capacity / tau
            Y[0, 1] = self.Htr_is
            Y[1, 0] = self.Htr_is
            Y[1, 1] = -(self.Htr_is + self.Htr_w + self.Htr_ms)
            Y[1, 2] = self.Htr_ms
            Y[2, 1] = self.Htr_ms
            Y[2, 2] = -self.Cm / tau - self.Htr_em - self.Htr_ms

            q[0] = -phi_HC_set * conv_factor - Hve_inf * T_e - Hve_vent * T_sup_AHU - phi_ia \
                   - self._air_thermal_capacity * self.Ta0 / tau
            q[1] = -phi_st * rad_factor - self.Htr_w * T_e
            q[2] = -self.Htr_em * T_e - phi_m - self.Cm * self.Tm0[0] / tau
            y = np.linalg.inv(Y).dot(q)
            # Seems to be more computationally efficient then np.insert
            return np.array([phi_HC_set] + [num for num in y]) # np.insert(y, 0, phi_HC_set)

        else:
            raise ValueError(f'Energy system zone solution: flag must be "phiset" or "Tset", flag: {flag}')

    def sensible_balance_2C(self, flag, Hve, T_e, T_e_eq, T_sup_AHU, phi_load, sigma=[0., 0., 1.], T_set=20.,
                            phi_HC_set=0.):
        """Sensible2C solves the linear system (Y*x = q) of VDI6007 at each
        iteration with a given setpoint temperature or
        (*)Note: if  phi_HC > 0 HEATING LOAD; phi_HC < 0 COOLING LOAD.

        Parameters
        ----------
        flag : string
            string 'Tset' or 'phiset'
        Hve : list
            list of two ositive floats. Ventilation and infiltration heat tranfer coeff [W/K]
            [Hve_vent, Hve_inf]
        T_e : float
            time step external temperature [°C]
        T_e_eq : float
            time step external equivalent temperature [°C]
        T_sup_AHU : float
            time step ventilation supply temperature [°C]
        phi_load : list
            list of three floats: internal and solar gains, three components (convective, aw and iw)[W]
            [phi_conv, phi_aw, phi_iw]
        sigma: list, default [0., 0., 1.]
            list of 3 floats
            portion of the heating cooling load to:
            sigma[0]: radiant to non adiabatic
            sigma[1]: radiant to adiabatic
            sigma[2]: convective to air node
            The sum must be 1
        T_set : float, default 20.
            time step setpoint of considered thermal zone [°C]
        phi_HC_set : float, default  0.
            time step Thermal load entering in the system [W]

        Returns
        -------
        numpy.array
            temperature nodes of the RC model:
            theta_m_aw thermal mass of AW building components [°C]
            theta_s_aw surface of AW building components [°C]
            theta_lu_star No physical meaning (node obtained from the delta-->star transformation) [°C]
            theta_I_lu internal air temperature [°C]  <--- WHAT WE WILL EXTRACT AS OUTPUT outside the function
            Q_hk_ges heating/cooling load for maintaining the given setpoint temperature [W]  <--- WHAT WE WILL EXTRACT AS OUTPUT outside the function
            theta_s_iw surface of IW building components [°C]
            theta_m_iw thermal mass of IW building components [°C]
            e.g.
            [
            theta_m_aw [°C],
            theta_s_aw [°C],
            theta_lu_star [°C],
            theta_I_lu [°C],
            Q_hk_ges [W]
            theta_s_iw [°C],
            theta_m_iw [°C],
            ]
        """

        # Check input data type

        if flag != 'Tset' and flag != 'phiset':
            raise TypeError(
                f"ERROR Thermal zone {self.name}, Sensible2C flag input is not a 'Tset' or 'phiset': flag {flag}")
        # if np.abs((np.array(sigma).sum() - 1)) > 1e-3:
        #     raise ValueError(
        #         f"Thermal Zone {self.name}, Sensible1C: sigma total must be 1. Sigma = {sigma}"
        #     )
        # % Resistances and capacitances of the 7R2C model
        R_lue_ve = 1e20 if Hve[0] == 0 else 1 / Hve[0]
        R_lue_inf = 1e20 if Hve[1] == 0 else 1 / Hve[1]

        Q_il_kon = phi_load[0]  # convective heat gains (internal)
        Q_il_str_aw = phi_load[1]  # radiant heat gains on surface node AW (internal + solar)
        Q_il_str_iw = phi_load[2]  # radiant heat gains on surface node IW (internal + solar)

        tau = CONFIG.time_step

        theta_A_eq = T_e_eq  # equivalent outdoor temperature (sol-air temperature)
        theta_lue = T_e  # outdoor air temperature
        theta_sup = T_sup_AHU

        if flag == 'Tset':
            theta_I_lu = T_set  # internal air setpoint

            # MATRIX OF THERMAL TRANSMITTANCES

            Y = np.zeros([6, 6])

            Y[0, 0] = -1 / self.RrestAW - 1 / self.R1AW - self.C1AW / tau
            Y[0, 1] = 1 / self.R1AW

            Y[1, 0] = 1 / self.R1AW
            Y[1, 1] = -1 / self.R1AW - 1 / self.RalphaStarAW
            Y[1, 2] = 1 / self.RalphaStarAW
            Y[1, 3] = sigma[1]

            Y[2, 1] = 1 / self.RalphaStarAW
            Y[2, 2] = -1 / self.RalphaStarAW - 1 / self.RalphaStarIL - 1 / self.RalphaStarIW
            Y[2, 4] = 1 / self.RalphaStarIW

            Y[3, 2] = 1 / self.RalphaStarIL
            Y[3, 3] = sigma[2]

            Y[4, 2] = 1 / self.RalphaStarIW
            Y[4, 3] = sigma[0]
            Y[4, 4] = -1 / self.RalphaStarIW - 1 / self.R1IW
            Y[4, 5] = 1 / self.R1IW

            Y[5, 4] = 1 / self.R1IW
            Y[5, 5] = -1 / self.R1IW - self.C1IW / tau

            # VECTOR OF KNOWN VALUES

            q = np.zeros(6)

            q[0] = -theta_A_eq / self.RrestAW - self.C1AW * self.Tm0[0] / tau
            q[1] = -Q_il_str_aw
            q[2] = -theta_I_lu / self.RalphaStarIL
            q[3] = theta_I_lu / self.RalphaStarIL - Q_il_kon \
                   - (theta_lue - theta_I_lu) / R_lue_inf \
                   - (theta_sup - theta_I_lu) / R_lue_ve \
                   + self._air_thermal_capacity * (theta_I_lu - self.Ta0) / tau
            q[4] = -Q_il_str_iw
            q[5] = -self.C1IW * self.Tm0[1] / tau

            # OUTPUT (UNKNOWN) VARIABLES OF THE LINEAR SYSTEM

            y = np.linalg.inv(Y).dot(q)
            # Seems to be more computationally efficient then np.insert
            return np.array([a for a in y[:3]] + [T_set] + [a for a in y[3:]]) # np.insert(y, 3, T_set)

        elif flag == 'phiset':
            # Note: the heat load in input is already distributed on the 3 nodes

            Q_hk_iw = phi_HC_set * sigma[0]  # % radiant heat flow from HVAC system (on surface node IW)
            Q_hk_aw = phi_HC_set * sigma[1]  # % radiant heat flow from HVAC system (on surface node AW)
            Q_hk_kon = phi_HC_set * sigma[2]  # % convective heat flow from HVAC system (on air node)

            # MATRIX OF THERMAL TRANSMITTANCES

            Y = np.zeros([6, 6])

            Y[0, 0] = -1 / self.RrestAW - 1 / self.R1AW - self.C1AW / tau
            Y[0, 1] = 1 / self.R1AW

            Y[1, 0] = 1 / self.R1AW
            Y[1, 1] = -1 / self.R1AW - 1 / self.RalphaStarAW
            Y[1, 2] = 1 / self.RalphaStarAW

            Y[2, 1] = 1 / self.RalphaStarAW
            Y[2, 2] = -1 / self.RalphaStarAW - 1 / self.RalphaStarIL - 1 / self.RalphaStarIW
            Y[2, 3] = 1 / self.RalphaStarIL
            Y[2, 4] = 1 / self.RalphaStarIW

            Y[3, 2] = 1 / self.RalphaStarIL
            Y[3, 3] = -1 / self.RalphaStarIL - 1 / R_lue_inf - 1 / R_lue_ve - self._air_thermal_capacity / tau

            Y[4, 2] = 1 / self.RalphaStarIW

            Y[4, 4] = -1 / self.RalphaStarIW - 1 / self.R1IW
            Y[4, 5] = 1 / self.R1IW

            Y[5, 4] = 1 / self.R1IW
            Y[5, 5] = -1 / self.R1IW - self.C1IW / tau

            # VECTOR OF KNOWN TERMS

            q = np.zeros(6)
            q[0] = -theta_A_eq / self.RrestAW - self.C1AW * self.Tm0[0] / tau
            q[1] = -Q_hk_aw - Q_il_str_aw
            q[2] = 0
            q[3] = -Q_hk_kon - Q_il_kon - theta_lue / R_lue_inf - theta_sup / R_lue_ve - self._air_thermal_capacity * self.Ta0 / tau
            q[4] = -Q_hk_iw - Q_il_str_iw
            q[5] = -self.C1IW * self.Tm0[1] / tau

            # OUTPUT LINEAR SYSTEM

            y = np.linalg.inv(Y).dot(q)
            # Seems to be more computationally efficient then np.insert
            return np.array([num for num in y[:4]] + [phi_HC_set] + [num for num in y[4:]]) # np.insert(y, 4, phi_HC_set)

        else:
            raise ValueError(f'Energy system zone solution: flag must be "phiset" or "Tset", flag: {flag}')

    def latent_balance(self, flag, G_ve, x_ext, x_sup, vapour_int_load, t_air_int, p_atm, rh_int_set=0.5,
                       phi_HC_set=0.):
        """Solves latent balance (vapour balance) for a specific time step

        Parameters
        ----------
        flag : str
             calculation mode"phiset" or "rhset"
        G_ve : list
            List of two floats: ventialtion and infiltration mass flow rate [kg/s]
            [G_ve_vent, G_ve_inf]
        x_ext : float
            time step external specific humidity [kg_vap/kg_da]
        x_sup : float
            time step supply specific humidity [kg_vap/kg_da]
        vapour_int_load : float
            time step vapour mass flow rate of internal loads [kg/s]
        t_air_int : float
            time step air temperature [°C]
        p_atm : float
            time step air temperature [°C]
        rh_int_set : float, default 0.5
            relative humidity set point [-]
        phi_HC_set : float, default 0.
            hvac latent load [W]

        Returns
        -------
        tuple
            array with zone specific humidity [kg_v, kg_as], zone relative humidity [-], latent demand [W]
            e.g.
            [
            zone specific humidity [kg_v, kg_as],
            zone relative humidity [-],
            latent demand [W],
            ]

        """

        x_int_set, p_intsat = self.get_specific_humidity(t_air_int, rh_int_set, p_atm)

        G_da_vent = G_ve[0]
        G_da_inf = G_ve[1]

        rho_air = air_properties['density']
        vapour_lat_heat = vapour_properties['latent_heat']
        vapour_spec_heat = vapour_properties['specific_heat']

        tau = CONFIG.time_step

        if flag == 'rhset':
            phi_lat = (G_da_inf * (x_ext - x_int_set) + G_da_vent * (
                    x_sup - x_int_set) - rho_air * self._volume * (x_int_set - self.xm0) / tau) * (
                              vapour_lat_heat + vapour_spec_heat * t_air_int) + vapour_int_load * (
                              vapour_lat_heat + vapour_spec_heat * t_air_int)
            x_int = x_int_set
            phi_lat = -phi_lat
        elif flag == 'phiset':
            x_int = (G_da_inf * x_ext + G_da_vent * x_sup + phi_HC_set / (
                    vapour_lat_heat + vapour_spec_heat * t_air_int) + vapour_int_load + rho_air * self._volume * self.xm0 / tau) / (
                            G_da_inf + G_da_vent + rho_air * self._volume / tau)
            phi_lat = phi_HC_set
        else:
            raise ValueError(f'Humidity system zone solution: flag must be "phiset" or "rhset", flag: {flag}')

        p_int = p_atm * x_int / (0.622 + x_int)
        RH_int = p_int / p_intsat
        return x_int, RH_int, phi_lat

    def reset_init_values(self, T:float = 15., X:float = 0.0105):
        '''This method allows to reset temperatures starting values, for the first calculation

        Parameters
        ----------
        T : float, default 15.
            Temperature to set as temperature mass [°C]
        X : float, default 0.0105
            Specific humidity [kg_v/kg_as]

        '''

        self.Ta0 = T
        self.Tm0 = np.array([T, T])
        self.xm0 = X
        self.zone_air_temperature = T
        self.zone_air_spec_humidity = X

    def reset_init_values_VDI(self):
        '''This method allows to reset temperatures starting values, according to VDI tests
        Starting T 22 °C
        Starting x 0.0105 kg_v/kg_as
        '''

        self.Ta0 = 22.
        self.Tm0 = np.array([22., 22.])
        self.xm0 = 0.0105
        self.zone_air_temperature = 22.
        self.zone_air_spec_humidity = 0.0105

    def solve_timestep(self, t, weather, model='2C'):
        '''Solves the thermal zone t - time step

        Parameters
        ----------
        t : int
            timestep of the simulation
        weather : eureca_building.weather.WeatherFile
            WeatherFile object
        model : str, default "2C"
            '1C' model or '2C' model
        '''

        # Check input data type

        if model != '1C' and model != '2C':
            raise TypeError(
                f"Thermal Zone {self.name}, time step {t}, model input is not a '1C' or '2C': model {model}")

        flag_AHU = True

        # Weather data
        T_ext = weather.hourly_data['out_air_db_temperature'][t]
        x_ext = weather.hourly_data['out_air_specific_humidity'][t]
        p_atm = weather.hourly_data['out_air_pressure'][t]


        # Internal Loads
        G_IHG_vapour = self.latent_load[t]  # kg_vap/s
        if model == "1C":
            phi_load = [self.phi_ia[t], self.phi_st[t], self.phi_m[t]]
        else:
            T_ext_eq = self.theta_eq_tot[t]
            phi_load = [self.Q_il_kon_I[t], self.Q_il_str_aw[t], self.Q_il_str_iw[t]]

        # Natural Ventilation
        # nat_vent_mass_flow = self.natural_ventilation.get_timestep_ventilation_mass_flow(t, self.zone_air_temperature, weather)
        # nv_outcomes = self.natural_ventilation.get_timestep_ventilation_mass_flow(t, self.zone_air_temperature, weather)
        # nat_vent_vol_flow = nv_outcomes[2]  #m3/s
        nat_vent_vol_flow = 0 if self.natural_ventilation is None else self.natural_ventilation.get_timestep_ventilation_mass_flow(t, self.zone_air_temperature, weather)
        nat_vent_mass_flow = nat_vent_vol_flow * air_properties['density']  # [kg/s]
        self.nat_vent_air_flow_rate[t] = nat_vent_mass_flow  # [kg/s]
        # self.nat_vent_info = {
        #     "airflow_rate" : {'kg/s' : np.zeros([CONFIG.number_of_time_steps]),
        #                       'm3/s' : np.zeros([CONFIG.number_of_time_steps]),
        #                       'L/s' : np.zeros([CONFIG.number_of_time_steps]),
        #                       'vol/h' : np.zeros([CONFIG.number_of_time_steps]),
        #                       },
        #     "windows_opening" : {'open_fraction' : np.zeros([CONFIG.number_of_time_steps]),
        #                          'open_area' : np.zeros([CONFIG.number_of_time_steps]),
        #                          },
        #     }
        self.nat_vent_info['airflow_rate']['kg/s'][t] = nat_vent_mass_flow
        self.nat_vent_info['airflow_rate']['m3/s'][t] = nat_vent_vol_flow
        self.nat_vent_info['airflow_rate']['L/s'][t] = nat_vent_vol_flow/1000
        self.nat_vent_info['airflow_rate']['m3/h'][t] = nat_vent_vol_flow*3600
        self.nat_vent_info['airflow_rate']['vol/h'][t] = nat_vent_vol_flow/self._volume*3600

        G_OA_nat_vent = self.infiltration_air_flow_rate[t] + nat_vent_mass_flow # kg/s outdoor air
        H_ve_nat_vent = G_OA_nat_vent * air_properties['specific_heat']  # W/K

        # Ventilation
        # TODO: Air handler to more zone?
        # Air Handling Unit
        self.air_handling_unit.air_handling_unit_calc(t, weather, self.zone_air_temperature, self.zone_air_spec_humidity)
        T_sup = self.air_handling_unit.T_sup
        x_sup = self.air_handling_unit.x_sup
        # for now the whole ahu flow rate to zone
        G_OA_mec_vent = self.air_handling_unit.air_flow_rate_kg_S[t] # kg/s supply air
        H_ve_mec_vent = G_OA_mec_vent * air_properties['specific_heat']  # W/K

        H_ve = [H_ve_mec_vent, H_ve_nat_vent]

        # Set points
        T_set_heat = self._temperature_setpoint.schedule_lower.schedule[t]
        T_set_cool = self._temperature_setpoint.schedule_upper.schedule[t]
        RH_set_int_H = self._humidity_setpoint.schedule_lower.schedule[t]
        RH_set_int_C = self._humidity_setpoint.schedule_upper.schedule[t]

        # Definition if the zone terminal is on or off for heating and cooling
        zone_equipment_heating_mode = True
        zone_equipment_cooling_mode = True
        zone_humidity_equipment_heating_mode = True
        zone_humidity_equipment_cooling_mode = True
        P_max = self.design_heating_system_power
        P_min = self.design_sensible_cooling_system_power

        # Start the zone system solution.. The logic depends on the heating cooling and latent sensible

        # while flag_AHU:

        # SENSIBLE HEAT LOAD CALCULATION

        # First calc without plant (phi_HC_set = 0)
        if model == '1C':
            pot, Ta, Ts, Tm = self.sensible_balance_1C(
                'phiset',
                H_ve,
                T_ext,
                T_sup,
                phi_load,
                sigma=[0.5,0.5],
                phi_HC_set=0.)
        elif model == '2C':
            Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.sensible_balance_2C(
                'phiset',
                H_ve,
                T_ext,
                T_ext_eq,
                T_sup,
                phi_load,
                sigma=[0.25,0.25,0.5],
                phi_HC_set=0.)

        # Check if heating mode and heating calc
        if Ta < T_set_heat and zone_equipment_heating_mode:
            # This means that we are under setpoint and hetaing is active
            if model == '1C':
                pot, Ta, Ts, Tm = self.sensible_balance_1C(
                    'Tset',
                    H_ve,
                    T_ext,
                    T_sup,
                    phi_load,
                    sigma=self.heating_sigma[model],
                    T_set=T_set_heat)
            elif model == '2C':
                Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.sensible_balance_2C(
                    'Tset',
                    H_ve,
                    T_ext,
                    T_ext_eq,
                    T_sup,
                    phi_load,
                    sigma=self.heating_sigma[model],
                    T_set=T_set_heat)
            if pot > P_max:
                # If the system reaches the maximum power
                if model == '1C':
                    pot, Ta, Ts, Tm = self.sensible_balance_1C(
                        'phiset',
                        H_ve,
                        T_ext,
                        T_sup,
                        phi_load,
                        sigma=self.heating_sigma[model],
                        phi_HC_set=P_max)
                elif model == '2C':
                    Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.sensible_balance_2C(
                        'phiset',
                        H_ve,
                        T_ext,
                        T_ext_eq,
                        T_sup,
                        phi_load,
                        sigma=self.heating_sigma[model],
                        phi_HC_set=P_max)

        # Check if cooling mode and cooling calc
        if Ta > T_set_cool and zone_equipment_cooling_mode:
            # This means that we are over cooling setpoint and cooling is active
            if model == '1C':
                pot, Ta, Ts, Tm = self.sensible_balance_1C(
                    'Tset',
                    H_ve,
                    T_ext,
                    T_sup,
                    phi_load,
                    sigma=self.cooling_sigma[model],
                    T_set=T_set_cool)
            elif model == '2C':
                Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.sensible_balance_2C(
                    'Tset',
                    H_ve,
                    T_ext,
                    T_ext_eq,
                    T_sup,
                    phi_load,
                    sigma=self.cooling_sigma[model],
                    T_set=T_set_cool)
            if pot < P_min:
                # If the system reaches the maximum cooling power
                if model == '1C':
                    pot, Ta, Ts, Tm = self.sensible_balance_1C(
                        'phiset',
                        H_ve,
                        T_ext,
                        T_sup,
                        phi_load,
                        sigma=self.cooling_sigma[model],
                        phi_HC_set=P_min)
                elif model == '2C':
                    Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.sensible_balance_2C(
                        'phiset',
                        H_ve,
                        T_ext,
                        T_ext_eq,
                        T_sup,
                        phi_load,
                        sigma=self.cooling_sigma[model],
                        phi_HC_set=P_min)

        # LATENT HEAT LOAD CALCULATION

        # First calc without plant (phi_HC_set = 0)
        x_int, rh_int, lat_pot = self.latent_balance('phiset',
                                                     [G_OA_mec_vent, G_OA_nat_vent],
                                                     x_ext,
                                                     x_sup,
                                                     G_IHG_vapour,
                                                     Ta,
                                                     p_atm,
                                                     phi_HC_set=0.)

        # Humidifying
        if rh_int < RH_set_int_H and zone_humidity_equipment_heating_mode:
            x_int, rh_int, lat_pot = self.latent_balance('rhset',
                                                         [G_OA_mec_vent, G_OA_nat_vent],
                                                         x_ext,
                                                         x_sup,
                                                         G_IHG_vapour,
                                                         Ta,
                                                         p_atm,
                                                         rh_int_set=RH_set_int_H,
                                                         )
            # TODO: Check if a maximum value is reached
        # Dehumidifying
        elif rh_int > RH_set_int_C and zone_humidity_equipment_cooling_mode:
            x_int, rh_int, lat_pot = self.latent_balance('rhset',
                                                         [G_OA_mec_vent, G_OA_nat_vent],
                                                         x_ext,
                                                         x_sup,
                                                         G_IHG_vapour,
                                                         Ta,
                                                         p_atm,
                                                         rh_int_set=RH_set_int_C,
                                                         )
            # TODO: Check if a maximum value is reached

        # err_T_sup = abs(T_sup - self.air_handling_unit.T_sup)
        # err_x_sup = abs(x_sup - self.air_handling_unit.x_sup)
        #
        # if err_T_sup > 0.1 or err_x_sup > 0.0001:
        #     #case of t_sup or x_sup changed
        #     T_sup = self.air_handling_unit.T_sup
        #     x_sup = self.air_handling_unit.x_sup
        # else:
        if isinstance(self.air_handling_unit, AirHandlingUnit):
            # UPDATE of dynamic variables
            self.air_handling_unit.supply_temperature.schedule[t] = T_sup
            self.air_handling_unit.supply_specific_humidity.schedule[t] = x_sup

        heat_flow = pot
        air_temp = Ta
        if model == '1C':
            self.Tm0 = [Tm, Tm]  # For 1C only 1 Tm
            operative_temp = (Ts + Ta)/2
            mean_radiant_temp = Ts
        elif model == '2C':
            self.Tm0 = [Tm_aw, Tm_iw]
            operative_temp = ((Ts_aw * self.Aaw_tot/self.Araum_tot + Ts_iw * ( 1- self.Aaw_tot/self.Araum_tot)) + Ta)/2
            mean_radiant_temp = (Ts_aw * self.Aaw_tot/self.Araum_tot + Ts_iw * ( 1- self.Aaw_tot/self.Araum_tot))
        self.Ta0 = Ta
        self.xm0 = x_int

        lat_heat_flow = lat_pot
        air_spec_humidity = x_int
        air_rel_humidity = rh_int
        flag_AHU = False

        self.sensible_zone_load = pot
        self.latent_zone_load = lat_heat_flow

        self.sensible_AHU_load = self.air_handling_unit.AHU_demand_sens
        self.latent_AHU_load = self.air_handling_unit.AHU_demand_lat
        self.AHU_electric_consumption = self.air_handling_unit.electric_consumption_W[t]

        self.zone_air_temperature = air_temp
        self.zone_operative_temperature = operative_temp
        self.zone_mean_radiant_temperature = mean_radiant_temp
        self.zone_air_rel_humidity = air_rel_humidity
        self.zone_air_spec_humidity = x_int

        return [air_temp,operative_temp,mean_radiant_temp], air_rel_humidity, self.Tm0, pot, lat_heat_flow

    def solve_quasisteadystate_method(self, weather):
        # TODO: Docstring
        # This runs the ISO13790 methods to calculates envelope params and monthly int and sol heat gains
        self._ISO13790_params()
        self.calculate_zone_loads_ISO13790(weather)

        # Calculation of the hours of heating cooling for each month
        heating_timesteps = self._temperature_setpoint.schedule_lower.schedule > 17.
        monthly_heating_timesteps = get_monthly_value_from_annual_vector(heating_timesteps, method = 'sum')
        cooling_timesteps = self._temperature_setpoint.schedule_upper.schedule < 27.
        monthly_cooling_timesteps = get_monthly_value_from_annual_vector(cooling_timesteps, method = 'sum')

        # Calculation of average setpoints

        with warnings.catch_warnings():
            # I expect to see RuntimeWarnings in this block
            warnings.simplefilter("ignore", category=RuntimeWarning)
            h_av_sp = self._temperature_setpoint.schedule_lower.schedule[heating_timesteps > 0.5].mean()
            c_av_sp = self._temperature_setpoint.schedule_upper.schedule[cooling_timesteps > 0.5].mean()

        h_av_sp = h_av_sp if monthly_heating_timesteps.sum() > 0 else 20
        c_av_sp = c_av_sp if monthly_cooling_timesteps.sum() > 0 else 26
        # [J] = [W/K] * ([°C] - [°C]) * [-] * [s]
        Q_h_tr_monthly = self.UA_tot * (h_av_sp - weather.monthly_data["out_air_db_temperature"]) * monthly_heating_timesteps * CONFIG.time_step
        # Q_h_tr_monthly[Q_h_tr_monthly < 0.] = 0.
        Q_c_tr_monthly = self.UA_tot * (c_av_sp - weather.monthly_data["out_air_db_temperature"]) * monthly_cooling_timesteps * CONFIG.time_step
        # Q_c_tr_monthly[Q_c_tr_monthly > 0.] = 0. # For cooling negative values

        # Infiltration
        # [W/K] = [kg/s] * [J/kgK]
        H_ve_inf = self.infiltration_air_flow_rate * air_properties['specific_heat']
        Q_h_inf = H_ve_inf * (h_av_sp - weather.hourly_data["out_air_db_temperature"]) * heating_timesteps # [W]
        Q_c_inf = H_ve_inf * (c_av_sp - weather.hourly_data["out_air_db_temperature"]) * cooling_timesteps # [W]
        Q_h_inf_monthly = get_monthly_value_from_annual_vector(Q_h_inf, method = 'integral') # [J]
        Q_c_inf_monthly = get_monthly_value_from_annual_vector(Q_c_inf, method = 'integral') # [J]

        # Effect of Mech Ventilation in zone
        # [W/K] = [kg/s] * [J/kgK]
        H_ve_mec = self.air_handling_unit.air_flow_rate_kg_S * air_properties['specific_heat']
        Q_h_ve = H_ve_mec * (h_av_sp - self.air_handling_unit.supply_temperature.schedule) * heating_timesteps # [W]
        Q_c_ve = H_ve_mec * (c_av_sp - self.air_handling_unit.supply_temperature.schedule) * cooling_timesteps # [W]
        Q_h_ve_monthly = get_monthly_value_from_annual_vector(Q_h_ve, method = 'integral') # [J]
        Q_c_ve_monthly = get_monthly_value_from_annual_vector(Q_c_ve, method = 'integral') # [J]

        Q_h_ht = Q_h_tr_monthly + Q_h_inf_monthly + Q_h_ve_monthly
        Q_c_ht = Q_c_tr_monthly + Q_c_inf_monthly + Q_c_ve_monthly

        Q_h_gn = self.monthly_int_gain + self.monthly_sol_gain
        Q_c_gn = Q_h_gn

        # Eta calculation
        # Heating
        tau = self.Cm/3600/(self.UA_tot + np.max(H_ve_inf))
        a_0 = 1.
        tau_0 = 15.
        a_h = a_0 + tau/tau_0
        Q_h_ht[Q_h_ht<1e-5]=1e-5
        gamma_h = Q_h_gn/Q_h_ht
        gamma_h[gamma_h>1e5] = 1e5
        eta_h_gn = np.zeros(12)
        eta_h_gn[gamma_h > 0.] = (1-gamma_h[gamma_h > 0.] ** a_h)/(1-gamma_h[gamma_h > 0.] ** (a_h+1))
        eta_h_gn[gamma_h == 1.] = (a_h)/(a_h+1)
        eta_h_gn[gamma_h <= 0.] = 1/gamma_h[gamma_h <= 0.]
        # Cooling
        a_c = a_0 + tau/tau_0
        Q_c_ht[Q_c_ht<1e-5]=1e-5
        gamma_c = Q_c_gn / Q_c_ht
        gamma_c[gamma_c>1e15] = 1e15
        eta_c_is = np.zeros(12)
        eta_c_is[gamma_c > 0.] = (1 - gamma_c[gamma_c > 0.] ** (-1*a_c)) / (1 - gamma_c[gamma_c > 0.] ** (-1*a_c - 1))
        eta_c_is[gamma_c == 1.] = (a_c) / (a_c + 1)
        eta_c_is[gamma_c <= 0.] = 1

        # Final monthly calc
        Q_h_sens = (Q_h_ht - eta_h_gn * Q_h_gn) * (monthly_heating_timesteps > 0)
        Q_c_sens = -1*(Q_c_gn - eta_c_is * Q_c_ht) * (monthly_cooling_timesteps > 0)
        Q_h_sens[Q_h_sens<0.] = 0.
        Q_c_sens[Q_c_sens>0.] = 0.


        # Latent
        humidifying_timesteps = self._humidity_setpoint.schedule_lower.schedule > 0.00725
        lat_int = self.extract_convective_radiative_latent_electric_load()['latent [kg_vap/s]']
        net_vapuor_in_zone_h = (self.infiltration_vapour_flow_rate + lat_int + self.air_handling_unit.air_flow_rate_kg_S * self.air_handling_unit.supply_specific_humidity.schedule) \
                             - (self.air_handling_unit.air_flow_rate_kg_S + self.infiltration_vapour_flow_rate) * 0.00725
        Q_hum_demand = -1 * net_vapuor_in_zone_h * vapour_properties['latent_heat']
        Q_hum_demand[~humidifying_timesteps] = 0.
        Q_hum_demand[Q_hum_demand < 0.] = 0.

        dehumidifying_timesteps = self._humidity_setpoint.schedule_upper.schedule < 0.0105
        lat_int = self.extract_convective_radiative_latent_electric_load()['latent [kg_vap/s]']
        net_vapuor_in_zone_c = (
                                           self.infiltration_vapour_flow_rate + lat_int + self.air_handling_unit.air_flow_rate_kg_S * self.air_handling_unit.supply_specific_humidity.schedule) \
                               - (
                                           self.air_handling_unit.air_flow_rate_kg_S + self.infiltration_vapour_flow_rate) * 0.0105
        Q_dehum_demand = -1 * net_vapuor_in_zone_c * vapour_properties['latent_heat']
        Q_dehum_demand[~dehumidifying_timesteps] = 0.
        Q_dehum_demand[Q_dehum_demand > 0.] = 0.

        Q_h_lat = get_monthly_value_from_annual_vector(Q_hum_demand, method = 'integral') # [J]
        Q_c_lat = get_monthly_value_from_annual_vector(Q_dehum_demand, method = 'integral') # [J]

        # Mechanical ventilation
        H_ve_mec = self.air_handling_unit.air_flow_rate_kg_S * air_properties['specific_heat']
        Q_h_mec_sens = H_ve_mec * (self.air_handling_unit.supply_temperature.schedule - weather.hourly_data['out_air_db_temperature']) # [W]
        Q_h_mec_sens = get_monthly_value_from_annual_vector(Q_h_mec_sens, method = 'integral') # [J]
        H_ve_mec = self.air_handling_unit.air_flow_rate_kg_S * vapour_properties['latent_heat']
        Q_h_mec_lat = H_ve_mec * (self.air_handling_unit.supply_specific_humidity.schedule - weather.hourly_data['out_air_specific_humidity']) # [W]
        Q_h_mec_lat = get_monthly_value_from_annual_vector(Q_h_mec_lat, method = 'integral') # [J]

        self.heat_sensible_zone_demand_qss_method = (Q_h_sens)/3600000 # [kWh]
        self.heat_latent_zone_demand_qss_method = (Q_h_lat)/3600000 # [kWh]
        self.cool_sensible_zone_demand_qss_method = (Q_c_sens)/3600000 # [kWh]
        self.cool_latent_zone_demand_qss_method = (Q_c_lat)/3600000 # [kWh]
        self.sensible_AHU_demand_qss_method = (Q_h_mec_sens)/3600000 # [kWh]
        self.latent_AHU_demand_qss_method = (Q_h_mec_lat)/3600000 # [kWh]

    def design_heating_load(self, weather):
        """Preliminary calculation to calculate the heating design temperature.
        Static calculation considering the product UA of all surfaces, the maximum infiltration flow rate, and the outdoor design temperature

        Parameters
        ----------
        t_ext_design : float
            outdoor design temperature [°C]
        """
        t_ext_design = weather.general_data['heating_oa_design_temperature']

        m_ve = self.infiltration_air_flow_rate.max()
        self.design_heating_system_power = 1.3 * (self.Htr_op + self.Htr_w) * (20. - t_ext_design) + m_ve * air_properties['specific_heat'] * (20. - t_ext_design)
        return self.design_heating_system_power

    def design_sensible_cooling_load(self, weather, model='2C'):
        """This method calculates the peak cooling according to the thermal load and weather of the thermal zone
        In order to do this, checks the hottest day (outdoor air temperature for 1C model, equivalent temperature for 2C, including also solar loads)
        And runs a simulation of the two days in between this maximum (similar to Carrier method, but using the 1C or 2C method)

        Parameters
        ----------
        weather : eureca_building.weather.WeatherFile
            WeatherFile object
        model : str, default "2C"
            '1C' model or '2C' model
        """
        self.reset_init_values(T = 27.)
        # Define point of max load
        if model == "1C":
            max = weather.hourly_data['out_air_db_temperature'].argmax()
        else:
            max = self.theta_eq_tot.argmax()
        lim_air = [max - 24*CONFIG.ts_per_hour, max + 24*CONFIG.ts_per_hour]
        sensible_max_cooling_load = 0.
        if model == "1C":
            for t in range(lim_air[0],lim_air[1]):
                pot, Ta, Ts, Tm = self.sensible_balance_1C("Tset",
                                         [0., self.infiltration_air_flow_rate[t] * air_properties['specific_heat']],
                                         weather.hourly_data['out_air_db_temperature'][t],
                                         26.,
                                         [self.phi_ia[t], self.phi_st[t], self.phi_m[t]],
                                         sigma = self.cooling_sigma["1C"],
                                         T_set = self._temperature_setpoint.schedule_upper.schedule.min(),
                                         )
                if pot < sensible_max_cooling_load:
                    sensible_max_cooling_load = pot

        if model == "2C":
            for t in range(lim_air[0],lim_air[1]):
                Tm_aw, Ts_aw, T_lu_star, Ta, pot, Ts_iw, Tm_iw = self.sensible_balance_2C(
                    'Tset',
                    [0., self.infiltration_air_flow_rate[t] * air_properties['specific_heat']],
                    weather.hourly_data['out_air_db_temperature'][t],
                    self.theta_eq_tot[t],
                    26.,
                    [self.Q_il_kon_I[t], self.Q_il_str_aw[t], self.Q_il_str_iw[t]],
                    sigma=self.cooling_sigma["2C"],
                    T_set=self._temperature_setpoint.schedule_upper.schedule.min()
                )
                if pot < sensible_max_cooling_load:
                    sensible_max_cooling_load = pot

        self.design_sensible_cooling_system_power = sensible_max_cooling_load * 1.2
        self.reset_init_values()
        return self.design_sensible_cooling_system_power

    def get_zone_info(self):
        """Just to print a beuty string of some info of the zone
        """
        return {
            "Zone volume [m3]": self._volume,
            "Zone net floor area [m2]": self._net_floor_area,
            'Number of units [-]': self.number_of_units,
            "Zone external walls opaque area [m2]": self.ext_wall_opaque_area,
            "Zone external walls glazed area [m2]": self.ext_wall_glazed_area,
            "Zone opaque heat transfer loss [W/K]": self.Htr_op,
            "Zone glazed heat transfer loss [W/K]": self.Htr_w,
            "Zone design heating load [kW]": self.design_heating_system_power/1000,
            "Zone design cooling load [kW]": self.design_sensible_cooling_system_power/1000,
        }


