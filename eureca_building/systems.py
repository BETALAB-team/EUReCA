"""
File with HVAC systems classes
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


import abc
import os
import io
import logging

import pandas as pd
import numpy as np

from eureca_building.systems_info import systems_info
from eureca_building.fluids_properties import fuels_pci
from eureca_building.config import CONFIG

# including Systems info from system_info json
global systems_info_dict
systems_info_dict = {}
for k,v in systems_info.items():
    systems_info_dict[k] = pd.read_csv(io.StringIO(v), sep = "\t")
    if "Type" in systems_info_dict[k].columns:
        systems_info_dict[k].set_index("Type", drop=True, inplace = True)

class System(metaclass=abc.ABCMeta):
    """
    Father class for HVAC systems classes
    Defines an interface
    """

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0

    @classmethod
    def __subclasshook__(cls, C):
        if cls is System:
            attrs = set(dir(C))
            if set(cls.__abstractmethods__) <= attrs:
                return True

        return NotImplemented

    @abc.abstractmethod
    def __init__(self, *args, **kwargs):
        # raise NotImplementedError(
        #     f"""
        # When creating a class for systems __init__, set, and solve methods must be created and overwritten
        #     """
        #     )
        pass

    @abc.abstractmethod
    def set_system_capacity(self):
        # raise NotImplementedError(
        #     f"""
        # When creating a class for systems __init__, set, and solve methods must be created and overwritten
        #     """
        #     )
        pass

    @abc.abstractmethod
    def solve_system(self):
        # raise NotImplementedError(
        #     f"""
        # When creating a class for systems __init__, set, and solve methods must be created and overwritten
        #     """
        #     )
        pass

    @property
    @abc.abstractmethod
    def electric_consumption(self):
        pass

    @property
    @abc.abstractmethod
    def gas_consumption(self):
        pass

    @property
    @abc.abstractmethod
    def wood_consumption(self):
        pass

    @property
    @abc.abstractmethod
    def oil_consumption(self):
        pass

    @property
    def system_type(self):
        return self._system_type

    @system_type.setter
    def system_type(self, value: str):
        if not isinstance(value, str):
            raise TypeError(
                f"HVAC System, system type not string: {type(value)}"
            )
        self._system_type = value


# %%---------------------------------------------------------------------------------------------------
# %% IdealLoad class

class IdealLoad(System):
    '''
    Class IdealLoad is for the ideal  zone balance. This actualy passes all methods.

    Methods:
        init
        set
        solve
    '''

    def __init__(self, *args, **kwargs):
        pass

    def set_system_capacity(self, design_power, weather):
        pass

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
        pass


# %%---------------------------------------------------------------------------------------------------
# %% CondensingBoiler class
class CondensingBoiler(System):
    '''
    Class CondensingBoiler
    This method considers a generic Condensing Boiler as the heating system
    of the entire building following UNI-TS 11300:2 - 2008

    Methods:
        init
        set
        solve
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0

    def __init__(self, *args, **kwargs):
        '''
        init. Set some attributes for the method

        Parameters
            ----------
            system_type : string
                string that sets the heating system
        Returns
        -------
        None.

        '''

        # Input Data
        self.system_type = "CondensingBoiler"
        self.theta_gn_test_Pn = 70  # [°C]
        self.theta_gn_test_Pint = 35  # [°C]
        self.theta_a_test = 20  # [°C]
        self.PCI_natural_gas = fuels_pci["Natural Gas"]  # Wh/Nm3
        self.system_info = systems_info_dict["CondensingBoiler"]

        # Vectors initialization
        # self.phi_gn_i_Px = np.zeros(l)  # [kW]
        # self.W_aux_Px = np.zeros(l)  # [kW]
        # self.phi_gn_i_Px_1 = np.zeros(l)  # [kW]
        # self.phi_gn_i_Px_2 = np.zeros(l)  # [kW]

    def set_system_capacity(self, design_power, weather):
        '''
        Sets system design power

        Parameters
            ----------
            design_power : float
                Design Heating power  [W]
            Weather : WeatherFile
                WeatherFile object

        Returns
        -------
        None.
        '''

        # Choise of system size based on estimated nominal Power
        Size_0 = 0
        self.design_power = design_power  # [W]
        for i in self.system_info.index:
            if Size_0 < self.design_power <= self.system_info['Size [kW]'][i] * 1000:
                self.ID = i
            Size_0 = self.system_info['Size [kW]'][i] * 1000
        if self.design_power > self.system_info['Size [kW]'][self.system_info.index[-1]] * 1000:
            self.ID = len(self.system_info) - 1
        self.Pint = self.system_info['P_int [%]'][self.ID] * self.design_power  # [W]
        if self.system_info['location [str]'][self.ID] == 'tech_room':
            self.theta_a_gn = 15  # [°C]
        elif self.system_info['location [str]'][self.ID] == 'internal':
            self.theta_a_gn = 20  # [°C]
        elif self.system_info['location [str]'][self.ID] == 'external':
            self.theta_a_gn = weather.general_data['average_out_air_db_temperature_heating_season']  # [°C]

        # Minimum Condensing Boiler efficiency check based on Directive
        # 92/42/CEE'''
        if self.design_power > 400000:
            self.design_power_check = 400000  # [W]
        else:
            self.design_power_check = self.design_power
        self.eta_gn_Pn_min = (91 + 1 * np.log(self.design_power_check / 1000))  # [%]
        self.eta_gn_Pint_min = (97 + 1 * np.log(self.design_power_check / 1000))  # [%]
        if self.system_info['eta_nom [%]'][self.ID] < self.eta_gn_Pn_min or \
                self.system_info['eta_int [%]'][self.ID] < self.eta_gn_Pint_min:
            logging.warning('Warning: Condensing Boiler efficiencies are lower than minimum')

        # No load losses estimation
        self.phi_gn_I_P0 = (self.design_power_check * 4.8 / 100 * (self.design_power_check / 1000) ** (-0.35))  # [W]

        # Auxiliary power input data
        self.W_aux_Pn = (45 * (self.design_power/1000) ** 0.48) # [W]
        self.W_aux_Pint = (15 * (self.design_power/1000) ** 0.48) # [W]
        self.W_aux_P0 = (15)  # [W]
        self.FC_Pint = self.Pint / self.design_power

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
        '''
        This method allows to calculate Condensing Boiler losses following
        the Standard UNI-TS 11300:2 - 2008

        Parameters
            ----------
            heat_flow : float
                required power  [W]
            Weather : WeatherFile
                WeatherFile object
            t : int
                Simulation time step
            T_int : float
                Zone temperature [°]
            RH_int : float
                Zone relative humidity [%]


        Returns
        -------
        None.

        '''
        # Corrected efficiency and losses at nominal power
        self.eta_gn_Pn_cor = self.system_info['eta_nom [%]'][self.ID] + \
                             self.system_info['f_cor_Pn [-]'][self.ID] * (
                                         self.theta_gn_test_Pn - self.system_info['T_gn_w [°C]'][self.ID])
        self.phi_gn_i_Pn_cor = (100 - self.eta_gn_Pn_cor) / self.eta_gn_Pn_cor * self.design_power  # [W]

        # Corrected efficiency and losses at intermadiate power
        self.eta_gn_Pint_cor = self.system_info['eta_int [%]'][self.ID] + \
                               self.system_info['f_cor_Pint [-]'][self.ID] * (
                                           self.theta_gn_test_Pint - self.system_info['T_gn_Pint [°C]'][self.ID])
        self.phi_gn_i_Pint_cor = (100 - self.eta_gn_Pint_cor) / self.eta_gn_Pint_cor * self.Pint  # [W]

        # Corrected no-load losses
        self.phi_gn_i_P0_cor = self.phi_gn_I_P0 * (
                    (self.system_info['T_gn_w [°C]'][self.ID] - self.theta_a_gn) / (self.theta_gn_test_Pn - self.theta_a_test)) ** 1.25

        # Total losses at current power
        if heat_flow <= 0:
            self.phi_gn_i_Px = 0
        else:
            self.phi_gn_i_Px_1 = ((heat_flow/1000) ** 2) * (((self.Pint/1000) * (
                        self.phi_gn_i_Pn_cor/1000 - self.phi_gn_i_P0_cor/1000) - self.design_power/1000 * (
                                                                    self.phi_gn_i_Pint_cor/1000 - self.phi_gn_i_P0_cor/1000)) / (
                                                                   self.design_power/1000 * self.Pint/1000 * (self.design_power/1000 - self.Pint/1000)))*1000 #[W]
            self.phi_gn_i_Px_2 = heat_flow/1000 * ((((self.design_power/1000) ** 2) * (self.phi_gn_i_Pint_cor/1000 - self.phi_gn_i_P0_cor/1000) - (
                    (self.Pint/1000) ** 2) * (self.phi_gn_i_Pn_cor/1000 - self.phi_gn_i_P0_cor/1000)) / (
                                                            self.design_power/1000 * self.Pint/1000 * (self.design_power/1000 - self.Pint/1000))) * 1000 # [W]
            self.phi_gn_i_Px = self.phi_gn_i_Px_1 + self.phi_gn_i_Px_2 + self.phi_gn_i_P0_cor  # [W]

        # Auxiliary power estimation
        if heat_flow <= 0:
            self.FC_Px = 0
            self.W_aux_Px = 0
        else:
            self.FC_Px = heat_flow / self.design_power
            if 0 < self.FC_Px <= self.FC_Pint:
                self.W_aux_Px = self.W_aux_P0 + self.FC_Px / self.FC_Pint * (self.W_aux_Pint - self.W_aux_P0)  # [W]
            elif self.FC_Pint < self.FC_Px <= 1:
                self.W_aux_Px = self.W_aux_Pint + (self.FC_Px - self.FC_Pint) * (self.W_aux_Pn - self.W_aux_Pint) / (
                            1 - self.FC_Pint)  # [W]
            elif self.FC_Px > 1:
                self.FC_Px = 1
                self.W_aux_Px = self.W_aux_Pint + (self.FC_Px - self.FC_Pint) * (self.W_aux_Pn - self.W_aux_Pint) / (
                            1 - self.FC_Pint)  # [W]

        total_energy = heat_flow + self.phi_gn_i_Px
        self.gas_consumption = total_energy / CONFIG.ts_per_hour /self.PCI_natural_gas
        self.electric_consumption = self.W_aux_Px / CONFIG.ts_per_hour

# %%---------------------------------------------------------------------------------------------------
# %% TraditionalBoiler class
class TraditionalBoiler(System):
    '''
    Class TraditionalBoiler
    This method considers a generic Traditional Boiler as the heating system
    of the entire building following UNI-TS 11300:2 - 2008

    Methods:
        init
        set
        solve
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0

    def __init__(self, *args, **kwargs):
        '''
        init. Set some attributes for the method

        Parameters
            ----------
            system_type : string
                string that sets the heating system
        Returns
        -------
        None.

        '''

        # Input Data
        self.system_type = "TraditionalBoiler"
        self.theta_gn_test_Pn = 70  # [°C]
        self.theta_gn_test_Pint = 50  # [°C]
        self.theta_a_test = 20  # [°C]
        self.PCI_natural_gas = fuels_pci["Natural Gas"]  # [Wh/Nm3]
        self.system_info = systems_info_dict["TraditionalBoiler"]

        # Vectors initialization
        # self.phi_gn_i_Px = np.zeros(l)  # [kW]
        # self.W_aux_Px = np.zeros(l)  # [kW]
        # self.phi_gn_i_Px_1 = np.zeros(l)  # [kW]
        # self.phi_gn_i_Px_2 = np.zeros(l)  # [kW]

    def set_system_capacity(self, design_power, weather):
        '''
        Sets system design power

        Parameters
            ----------
            design_power : float
                Design Heating power  [W]
            Weather : WeatherFile
                WeatherFile object

        Returns
        -------
        None.
        '''

        # Choise of system size based on estimated nominal Power
        Size_0 = 0
        self.design_power = design_power  # [W]
        for i in self.system_info.index:
            if Size_0 < self.design_power <= self.system_info['Size [kW]'][i] * 1000:
                self.ID = i
            Size_0 = self.system_info['Size [kW]'][i] * 1000
        if self.design_power > self.system_info['Size [kW]'][self.system_info.index[-1]] * 1000:
            self.ID = len(self.system_info) - 1
        self.Pint = self.system_info['P_int [%]'][self.ID] * self.design_power  # [W]
        if self.system_info['location [str]'][self.ID] == 'tech_room':
            self.theta_a_gn = 15  # [°C]
        elif self.system_info['location [str]'][self.ID] == 'internal':
            self.theta_a_gn = 20  # [°C]
        elif self.system_info['location [str]'][self.ID] == 'external':
            self.theta_a_gn = weather.general_data['average_out_air_db_temperature_heating_season']  # [°C]

        # Minimum Traditional Boiler efficiency check based on Directive
        # 92/42/CEE'''
        if self.design_power > 400000:
            self.design_power_check = 400000  # [W]
        else:
            self.design_power_check = self.design_power
        self.eta_gn_Pn_min = (84 + 2 * np.log(self.design_power_check / 1000))  # [%]
        self.eta_gn_Pint_min = (80 + 3 * np.log(self.design_power_check / 1000))  # [%]
        if self.system_info['eta_nom [%]'][self.ID] < self.eta_gn_Pn_min or \
                self.system_info['eta_int [%]'][self.ID] < self.eta_gn_Pint_min:
            logging.warning('Warning: Traditional Boiler efficiencies are lower than minimum')

        # No load losses estimation
        self.phi_gn_I_P0 = (self.design_power_check / 1000 * 10 * 8.5 * (self.design_power_check / 1000) ** (-0.4))  # [W]

        # Auxiliary power input data
        self.W_aux_Pn = (45 * (self.design_power / 1000) ** 0.48)   # [W]
        self.W_aux_Pint = (15 * (self.design_power/ 1000) ** 0.48)  # [W]
        self.W_aux_P0 = (15)   # [W]
        self.FC_Pint = self.Pint / self.design_power

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
        '''
        This method allows to calculate Traditional Boiler losses following
        the Standard UNI-TS 11300:2 - 2008

        Parameters
            ----------
            heat_flow : float
                required power  [W]
            Weather : WeatherFile
                WeatherFile object
            t : int
                Simulation time step
            T_int : float
                Zone temperature [°]
            RH_int : float
                Zone relative humidity [%]

        Returns
        -------
        None.

        '''
        # Corrected efficiency and losses at nominal power
        self.eta_gn_Pn_cor = self.system_info['eta_nom [%]'][self.ID] + \
                             self.system_info['f_cor_Pn [-]'][self.ID] * (
                                         self.theta_gn_test_Pn - self.system_info['T_gn_w [°C]'][self.ID])
        self.phi_gn_i_Pn_cor = (100 - self.eta_gn_Pn_cor) / self.eta_gn_Pn_cor * self.design_power  # [W]

        # Corrected efficiency and losses at intermadiate power
        self.eta_gn_Pint_cor = self.system_info['eta_int [%]'][self.ID] + \
                               self.system_info['f_cor_Pint [-]'][self.ID] * (
                                           self.theta_gn_test_Pint - self.system_info['T_gn_Pint [°C]'][self.ID])
        self.phi_gn_i_Pint_cor = (100 - self.eta_gn_Pint_cor) / self.eta_gn_Pint_cor * self.Pint  # [W]

        # Corrected no-load losses
        self.phi_gn_i_P0_cor = self.phi_gn_I_P0 * (
                    (self.system_info['T_gn_w [°C]'][self.ID] - self.theta_a_gn) / (self.theta_gn_test_Pn - self.theta_a_test)) ** 1.25

        # Total losses at current power
        if heat_flow <= 0:
            self.phi_gn_i_Px = 0
        else:
            self.phi_gn_i_Px_1 = ((heat_flow/1000) ** 2) * (
                    (
                            (self.Pint/1000) * (
                        self.phi_gn_i_Pn_cor/1000 - self.phi_gn_i_P0_cor/1000) - self.design_power/1000 * (
                                                                    self.phi_gn_i_Pint_cor/1000 - self.phi_gn_i_P0_cor/1000))
                    / (self.design_power/1000 * self.Pint/1000 * (self.design_power/1000 - self.Pint/1000)))*1000 #[W]
            self.phi_gn_i_Px_2 = heat_flow/1000 * \
                                 ((
                                          ((self.design_power/1000) ** 2) *(self.phi_gn_i_Pint_cor/1000 - self.phi_gn_i_P0_cor/1000)
                                          - ((self.Pint/1000) ** 2) * (self.phi_gn_i_Pn_cor/1000 - self.phi_gn_i_P0_cor/1000))
                                  / (self.design_power/1000 * self.Pint/1000 * (self.design_power/1000 - self.Pint/1000))) * 1000 # [W]
            self.phi_gn_i_Px = self.phi_gn_i_Px_1 + self.phi_gn_i_Px_2 + self.phi_gn_i_P0_cor  # [W]

        # Auxiliary power estimation
        if heat_flow <= 0:
            self.FC_Px = 0
            self.W_aux_Px = 0
        else:
            self.FC_Px = heat_flow / self.design_power
            if 0 < self.FC_Px <= self.FC_Pint:
                self.W_aux_Px = self.W_aux_P0 + self.FC_Px / self.FC_Pint * (self.W_aux_Pint - self.W_aux_P0)  # [W]
            elif self.FC_Pint < self.FC_Px <= 1:
                self.W_aux_Px = self.W_aux_Pint + (self.FC_Px - self.FC_Pint) * (self.W_aux_Pn - self.W_aux_Pint) / (
                            1 - self.FC_Pint)  # [W]
            elif self.FC_Px > 1:
                self.FC_Px = 1
                self.W_aux_Px = self.W_aux_Pint + (self.FC_Px - self.FC_Pint) * (self.W_aux_Pn - self.W_aux_Pint) / (
                            1 - self.FC_Pint)  # [W]


        total_energy = heat_flow + self.phi_gn_i_Px
        self.gas_consumption = total_energy / CONFIG.ts_per_hour /self.PCI_natural_gas
        self.electric_consumption = self.W_aux_Px / CONFIG.ts_per_hour

# %%---------------------------------------------------------------------------------------------------
# %% SplitAirCooler class

class SplitAirCooler(System):
    '''
    SplitAirCooler class

    This method considers a generic Split Air Conditioner as the cooling system
    of the entire building following UNI-TS 11300:3 - 2010

    Methods:
        init
        setCondensingBoiler
        solveSystem
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0

    def __init__(self, *args, **kwargs):
        '''
        init. Set some attributes for the method

        Parameters
            ----------
            system_type : string
                string that sets the cooling system

        Returns
        -------
        None.

        '''

        # Check input data type

        # The quality control is implemented in the System class.. Here is not necessary

        # Input Data
        self.C_system_type = "SplitAirCooler"
        self.T_ext_rif_100 = 35  # [°C]
        self.T_ext_rif_75 = 30  # [°C]
        self.T_ext_rif_50 = 25  # [°C]
        self.T_ext_rif_25 = 20  # [°C]
        self.T_int_rif = 27  # [°C]
        self.W_aux_gn = 0.04  # [kW/kW_cond]
        self.Conversion_fac = 2.17
        self.system_info = systems_info_dict["SplitAirCooler"]

        # Vectors initialization
        # self.W_el = np.zeros(l)  # [kWe]
        # self.W_aux = np.zeros(l)  # [kWe]
        # self.H_waste = np.zeros(l)  # [kW]
        # self.eta_mm = np.zeros(l)
        # self.eta_1 = np.zeros(l)

    def set_system_capacity(self, design_power, weather):
        '''
        Choise of system size based on estimated nominal Power

        Parameters
            ----------
            design_power : float
                Design Heating power  [W]
            Weather : WeatherFile
                WeatherFile object

        Returns
        -------
        None.
        '''

        Size_0 = 0.
        self.design_power = -1 * design_power  # [W]
        for i in self.system_info.index:
            if Size_0 <= self.design_power <= self.system_info['Size [kW]'][i] * 1000:
                self.ID = i
            Size_0 = self.system_info['Size [kW]'][i] * 1000
        if self.design_power > self.system_info['Size [kW]'][self.system_info.index[-1]] * 1000:
            self.ID = len(self.system_info) - 1
        self.Q_max = design_power  # [W]

        # EER curve in case of Air-to-Air unit
        self.EER_100 = self.system_info['EER_100 [-]'][self.ID]
        self.EER_75 = self.system_info['EER_75 [-]'][self.ID]
        self.EER_50 = self.system_info['EER_50 [-]'][self.ID]
        self.EER_25 = self.system_info['EER_25 [-]'][self.ID]
        self.EER_20 = self.EER_25 * 0.94
        self.EER_15 = self.EER_25 * 0.85
        self.EER_10 = self.EER_25 * 0.73
        self.EER_5 = self.EER_25 * 0.50
        self.EER_2 = self.EER_25 * 0.26
        self.EER_1 = self.EER_25 * 0.14
        self.x = np.array([1, 0.75, 0.5, 0.25, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01, 0.0])
        self.y = np.array(
            [self.EER_100, self.EER_75, self.EER_50, self.EER_25, self.EER_20, self.EER_15, self.EER_10, self.EER_5,
             self.EER_2, self.EER_1, 0.0])

        self.prosp_C_100 = np.array([[1.634, 1.457, 1.281, 1.105, 0.928, 0.807, 0.685], \
                                     [1.720, 1.518, 1.327, 1.148, 0.979, 0.856, 0.736], \
                                     [1.756, 1.534, 1.332, 1.155, 1.000, 0.871, 0.750], \
                                     [1.782, 1.569, 1.370, 1.187, 1.018, 0.890, 0.767], \
                                     [1.834, 1.639, 1.444, 1.249, 1.054, 0.928, 0.802]])
        self.prosp_C_75 = np.array([[1.457, 1.281, 1.105, 0.928, 0.807, 0.685, 0.585], \
                                    [1.518, 1.327, 1.148, 0.979, 0.856, 0.736, 0.636], \
                                    [1.534, 1.332, 1.155, 1.000, 0.871, 0.750, 0.650], \
                                    [1.569, 1.370, 1.187, 1.018, 0.890, 0.767, 0.667], \
                                    [1.639, 1.444, 1.249, 1.054, 0.928, 0.802, 0.700]])
        self.prosp_C_50 = np.array([[1.281, 1.105, 0.928, 0.807, 0.685, 0.585, 0.505], \
                                    [1.327, 1.148, 0.979, 0.856, 0.736, 0.636, 0.556], \
                                    [1.332, 1.155, 1.000, 0.871, 0.750, 0.650, 0.672], \
                                    [1.370, 1.187, 1.018, 0.890, 0.767, 0.667, 0.587], \
                                    [1.444, 1.249, 1.054, 0.928, 0.802, 0.700, 0.698]])
        self.prosp_C_25 = np.array([[1.062, 0.962, 0.871, 0.788, 0.714, 0.646, 0.585], \
                                    [1.083, 0.981, 0.888, 0.804, 0.728, 0.659, 0.596], \
                                    [1.105, 1.000, 0.905, 0.820, 0.742, 0.672, 0.608], \
                                    [1.126, 1.020, 0.923, 0.836, 0.757, 0.685, 0.620], \
                                    [1.149, 1.040, 0.941, 0.852, 0.771, 0.698, 0.632]])
        self.prosp_C_0 = self.prosp_C_25
        self.prosp_C = np.array([self.prosp_C_100, self.prosp_C_75, self.prosp_C_50, self.prosp_C_25, self.prosp_C_0])
        self.fx = np.array([1, 0.75, 0.5, 0.25, 0])
        self.pcfx = np.array([15, 20, 25, 30, 35, 40, 45])
        self.pcfy = np.array([16, 18, 19, 20, 22])

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
        '''
        This method allows to calculate Split Air Cooler electrical power
        following the Standard UNI-TS 11300:3 - 2010

        Parameters
            ----------
            heat_flow : float
                required power  [W]
            Weather : WeatherFile
                WeatherFile object
            t : int
                Simulation time step
            T_int : float
                Zone temperature [°]
            RH_int : float
                Zone relative humidity [%]

        Returns
        -------
        None.
        '''

        # Check input data type

        # The quality control is implemented in the System class.. Here is not necessary

        self.T_ext_solve = weather.hourly_data["out_air_db_temperature"][t]

        if self.T_ext_solve <= 15:
            self.T_ext_solve = 15.1
        elif self.T_ext_solve > 45:
            self.T_ext_solve = 45

        self.prosp_C_F = np.zeros([5, 7])

        # Wet bulb temperature
        self.T_wb = T_int * np.arctan(0.151977 * (RH_int + 8.313659) ** (0.5)) + np.arctan(T_int + RH_int) - np.arctan(
            RH_int - 1.676331) + 0.00391838 * (RH_int) ** (3 / 2) * np.arctan(0.023101 * RH_int) - 4.686035
        if self.T_wb <= 16:
            self.T_wb = 16.1
        elif self.T_wb > 22:
            self.T_wb = 22

        # Load Factor and EER in real conditions estimation
        if heat_flow >= 0:
            self.W_el = 0.  # [We]
            self.H_waste = 0.  # [W]
            self.W_aux = 0.  # [We]
        else:
            self.F = heat_flow / (self.Q_max)
            xa = self.x[1]
            xb = self.x[0]
            self.EER = (self.F - xb) / (xa - xb) * self.y[0] - (self.F - xa) / (xa - xb) * self.y[0]
            for i in range(len(self.x) - 1):
                xa = self.x[i + 1]
                xb = self.x[i]
                if xa <= self.F <= xb:
                    self.EER = (self.F - xb) / (xa - xb) * self.y[i + 1] - (self.F - xa) / (xa - xb) * self.y[i]
                elif self.F > 1:
                    self.EER = self.EER_100
                else:
                    pass
            for i in range(len(self.fx) - 1):
                n_prosp_a = i + 1
                n_prosp_b = i
                fxa = self.fx[i + 1]
                fxb = self.fx[i]
                if fxa < self.F <= fxb:
                    for n in range(5):
                        for m in range(7):
                            self.prosp_C_F[n][m] = (self.F - fxb) / (fxa - fxb) * self.prosp_C[n_prosp_a][n][m] - (
                                        self.F - fxa) / (fxa - fxb) * self.prosp_C[n_prosp_b][n][m]
                elif self.F > 1:
                    self.prosp_C_F = self.prosp_C_100
                else:
                    pass
            for n in range(4):
                y1 = self.pcfy[n]
                y2 = self.pcfy[n + 1]
                if y1 < self.T_wb <= y2:
                    for m in range(6):
                        x1 = self.pcfx[m]
                        x2 = self.pcfx[m + 1]
                        if x1 < self.T_ext_solve <= x2:
                            self.eta_1= 1 / ((x2 - x1) * (y2 - y1)) * (
                                        self.prosp_C_F[n + 1][m + 1] * (x2 - self.T_ext_solve) * (y2 - self.T_wb) +
                                        self.prosp_C_F[n][m + 1] * (self.T_ext_solve - x1) * (y2 - self.T_wb) +
                                        self.prosp_C_F[n + 1][m] * (x2 - self.T_ext_solve) * (self.T_wb - y1) +
                                        self.prosp_C_F[n][m] * (self.T_ext_solve - x1) * (self.T_wb - y1))

            # Electrical power required in real conditions
            self.eta_mm = self.EER * self.eta_1
            self.W_el = abs(heat_flow) / self.eta_mm # [We]
            self.H_waste = abs(heat_flow) + self.W_el # [W]
            self.W_aux = self.W_aux_gn * (self.H_waste)  # [We]

        self.electric_consumption = (self.W_el + self.W_aux) / CONFIG.ts_per_hour


# %%---------------------------------------------------------------------------------------------------
# %% ChillerAirtoWater class

class ChillerAirtoWater(System):
    '''
    ChillerAirtoWater class

    This method considers a generic Air-to-water Chiller as the cooling system
    of the entire building following UNI-TS 11300:3 - 2010


    Methods:
        init
        setCondensingBoiler
        solveSystem
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0

    def __init__(self, *args, **kwargs):

        '''
        init. Set some attributes for the method

        Parameters
            ----------
            system_type : string
                string that sets the heating system

        Returns
        -------
        None.

        '''

        # Check input data type

        # The quality control is implemented in the System class.. Here is not necessary

        # Input Data
        self.C_system_type = "ChillerAirtoWater"
        self.T_ext_rif_100 = 35  # [°C]
        self.T_ext_rif_75 = 30  # [°C]
        self.T_ext_rif_50 = 25  # [°C]
        self.T_ext_rif_25 = 20  # [°C]
        self.T_w_out_rif = 7  # [°C]
        self.W_aux_gn = 0.04  # [kW/kW_cond]
        self.Conversion_fac = 2.17
        self.system_info = systems_info_dict["ChillerAirtoWater"]

        # Vectors initialization
        # self.W_el = np.zeros(l)  # [kWe]
        # self.W_aux = np.zeros(l)  # [kWe]
        # self.H_waste = np.zeros(l)  # [kW]
        # self.eta_mm = np.zeros(l)
        # self.eta_1 = np.zeros(l)

    def set_system_capacity(self, design_power, weather):
        '''
        Choise of system size based on estimated nominal Power

        Parameters
            ----------
            design_power : float
                Design Heating power  [W]
            Weather : WeatherFile
                WeatherFile object

        Returns
        -------
        None.
        '''

        # Check input data type

        # The quality control is implemented in the System class.. Here is not necessary

        # Choise of system size based on estimated nominal Power
        Size_0 = 0
        self.Pnom = abs(design_power)  # [W]
        for i in self.system_info.index:
            if Size_0 <= self.Pnom <= self.system_info['Size [kW]'][i] * 1000:
                self.ID = i
            Size_0 = self.system_info['Size [kW]'][i] * 1000
        if self.Pnom > self.system_info['Size [kW]'][self.system_info.index[-1]] * 1000:
            self.ID = len(self.system_info) - 1
        self.Q_max = design_power  # [W]

        # EER curve in case of Air-to-Water chiller
        self.EER_100 = self.system_info['EER_100 [-]'][self.ID]
        self.EER_75 = self.system_info['EER_75 [-]'][self.ID]
        self.EER_50 = self.system_info['EER_50 [-]'][self.ID]
        self.EER_25 = self.system_info['EER_25 [-]'][self.ID]
        self.EER_20 = self.EER_25 * 0.95
        self.EER_15 = self.EER_25 * 0.94
        self.EER_10 = self.EER_25 * 0.87
        self.EER_5 = self.EER_25 * 0.71
        self.EER_2 = self.EER_25 * 0.46
        self.EER_1 = self.EER_25 * 0.29
        self.x = np.array([1, 0.75, 0.5, 0.25, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01, 0.0])
        self.y = np.array(
            [self.EER_100, self.EER_75, self.EER_50, self.EER_25, self.EER_20, self.EER_15, self.EER_10, self.EER_5,
             self.EER_2, self.EER_1, 0.0])

        # Assumption: T_water_out = 7°C and delta(T) = 5°C
        self.prosp_C_100 = np.array([1.756, 1.534, 1.332, 1.155, 1.000, 0.871, 0.750])
        self.prosp_C_75 = np.array([1.534, 1.332, 1.155, 1.000, 0.871, 0.750, 0.650])
        self.prosp_C_50 = np.array([1.332, 1.155, 1.000, 0.871, 0.750, 0.650, 0.570])
        self.prosp_C_25 = np.array([1.155, 1.000, 0.871, 0.750, 0.650, 0.570, 0.500])
        self.prosp_C_0 = self.prosp_C_25
        self.prosp_C = np.array([self.prosp_C_100, self.prosp_C_75, self.prosp_C_50, self.prosp_C_25, self.prosp_C_0])
        self.fx = np.array([1, 0.75, 0.5, 0.25, 0])
        self.pcfx = np.array([15, 20, 25, 30, 35, 40, 45])

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
        '''
        This method allows to calculate Air to Water Chiller electrical power
        following the Standard UNI-TS 11300:3 - 2010


        Parameters
            ----------
            heat_flow : float
                required power  [W]
            Weather : WeatherFile
                WeatherFile object
            t : int
                Simulation time step
            T_int : float
                Zone temperature [°]
            RH_int : float
                Zone relative humidity [%]

        Returns
        -------
        None.
        '''

        # Check input data type

        # The quality control is implemented in the System class.. Here is not necessary

        self.T_ext_solve = weather.hourly_data["out_air_db_temperature"][t]
        if self.T_ext_solve <= 15:
            self.T_ext_solve = 15.1
        elif self.T_ext_solve > 45:
            self.T_ext_solve = 45

        self.prosp_C_F = np.zeros(7)

        # Load Factor and EER in real conditions estimation
        if heat_flow >= 0:
            self.W_el = 0.  # [We]
            self.H_waste = 0.  # [W]
            self.W_aux = 0.  # [We]
        else:
            self.F = heat_flow / (self.Q_max)
            xa = self.x[1]
            xb = self.x[0]
            self.EER = (self.F - xb) / (xa - xb) * self.y[0] - (self.F - xa) / (xa - xb) * self.y[0]
            for i in range(len(self.x) - 1):
                xa = self.x[i + 1]
                xb = self.x[i]
                if xa <= self.F <= xb:
                    self.EER = (self.F - xb) / (xa - xb) * self.y[i + 1] - (self.F - xa) / (xa - xb) * self.y[i]
                elif self.F > 1:
                    self.EER = self.EER_100
                else:
                    pass
            for i in range(len(self.fx) - 1):
                fxa = self.fx[i + 1]
                fxb = self.fx[i]
                if fxa < self.F <= fxb:
                    for m in range(7):
                        self.prosp_C_F[m] = (self.F - fxb) / (fxa - fxb) * self.prosp_C[i + 1][m] - (self.F - fxa) / (
                                    fxa - fxb) * self.prosp_C[i][m]
                elif self.F > 1:
                    self.prosp_C_F = self.prosp_C_100
                else:
                    pass
            for m in range(6):
                x1 = self.pcfx[m]
                x2 = self.pcfx[m + 1]
                if x1 < self.T_ext_solve <= x2:
                    self.eta_1 = (self.T_ext_solve - x1) / (x2 - x1) * self.prosp_C_F[m + 1] - (
                                self.T_ext_solve - x2) / (x2 - x1) * self.prosp_C_F[m]

            # Electrical power required in real conditions
            self.eta_mm = self.EER * self.eta_1
            self.W_el = abs(heat_flow) / self.eta_mm # [We]
            self.H_waste = abs(heat_flow) + self.W_el # [W]
            self.W_aux = self.W_aux_gn * (self.H_waste)  # [We]

        self.electric_consumption = (self.W_el + self.W_aux) / CONFIG.ts_per_hour

# %%---------------------------------------------------------------------------------------------------
# %% SplitAirConditioner class
class SplitAirConditioner(System):
    '''
    SplitAirConditioner class

    This method considers a generic Split Air Conditioner as the cooling system
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
        solveSystem
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0

    def __init__(self, *args, **kwargs):
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
            system_type : string
                string that sets the heating system

        Returns
        -------
        None.

        '''

        # Check input data type

        # The quality control is implemented in the System class.. Here is not necessary

        # Input Data
        self.system_type = "SplitAirConditioner"  # system_type è una stringa
        self.T_ext_rif = 35  # [°C]
        self.T_int_rif = 27  # [°C]
        self.W_aux_gn = 0.04  # [kW/kW_cond]
        self.system_info = systems_info_dict["SplitAirConditioner"]
        # self.W_el = np.zeros(l)  # [kWe]
        # self.W_aux = np.zeros(l)  # [kWe]
        # self.H_waste = np.zeros(l)  # [kW]
        # self.eta_mm = np.zeros(l)
        # self.EER = np.zeros(l)
        # self.F = np.zeros(l)

    def set_system_capacity(self, design_power, weather):  # da modificare dentro il file excel di impianti

        '''
        Choise of system size based on estimated nominal Power for air conditioner WITHOUT inverter

        Parameters
            ----------
            design_power : float
                Design Heating power  [W]
            Weather : WeatherFile
                WeatherFile object

        Returns
        -------
        None.
        '''

        # Check input data type

        # The quality control is implemented in the System class.. Here is not necessary

        Size_0 = 0
        self.Pnom = abs(design_power)  # [W]
        for i in self.system_info.index:
            if Size_0 < self.Pnom <= self.system_info['Size [kW]'][i] * 1000:
                self.ID = i
            Size_0 = self.system_info['Size [kW]'][i]  * 1000
        if self.Pnom > self.system_info['Size [kW]'][self.system_info.index[-1]] * 1000:
            self.ID = len(self.system_info) - 1
        self.Q_max = design_power  # [W]

        if self.system_info['inverter [str]'][self.ID] == 'no':
            self.EER_targa = self.system_info['EER_medio [-]'][self.ID]
        elif self.system_info['inverter [str]'][self.ID] == 'yes':  # Choise of system size based on estimated nominal Power for air conditioner WITH inverter
            self.a = 11.2
            self.b = -0.363
            self.c = 2.6 * pow(10, -3)
            self.d = 0.829
            self.e = 6.9 * pow(10, -3)
            self.f = -0.159

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
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
            heat_flow : float
                required power  [W]
            Weather : WeatherFile
                WeatherFile object
            t : int
                Simulation time step
            T_int : float
                Zone temperature [°]
            RH_int : float
                Zone relative humidity [%]

        Returns
        -------
        None.
        '''

        # Check input data type

        # The quality control is implemented in the System class.. Here is not necessary
        T_ext = weather.hourly_data["out_air_db_temperature"][t]
        if heat_flow >= 0:
            self.W_el = 0.  # [We]
            self.H_waste = 0.  # [W]
        else:
            self.F = heat_flow / (self.Q_max)
            if self.system_info['inverter [str]'][self.ID] == 'no':
                self.T_ext_solve = T_ext
                if self.T_ext_solve <= 20:
                    self.T_ext_solve = 20.1
                elif self.T_ext_solve > 35:
                    self.T_ext_solve = 35
                self.EER = -0.14 * (self.T_ext_solve - 35) + self.EER_targa

            elif self.system_info['inverter [str]'][self.ID] == 'yes':  # Choise of system size based on estimated nominal Power for air conditioner WITH inverter'
                self.Pnom_cond = 7000  # [W] P_nom conditioner
                self.number_cond = heat_flow / self.Pnom_cond  # number of air conditioners suppose in each thermal zone
                self.EER = self.a + self.b * T_ext + self.c * T_ext ** 2 + self.d * (
                            heat_flow / self.number_cond / 1000) + self.e * T_ext * (heat_flow / self.number_cond / 1000) + self.f * (
                                          heat_flow / self.number_cond  / 1000) ** 2

            self.H_waste = (abs(heat_flow) * (1 + self.EER)) / self.EER
            self.W_el = self.H_waste - abs(heat_flow)

        self.electric_consumption = (self.W_el) / CONFIG.ts_per_hour

class Heating_EN15316(System):
    '''
    Class TraditionalBoiler
    This method considers a generic Traditional Boiler as the heating system
    of the entire building following UNI-TS 11300:2 - 2008

    Methods:
        init
        set
        solve
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0

    def __init__(self, *args, **kwargs):
        '''
        init. Set some attributes for the method

        Parameters
            ----------
            system_type : string
                string that sets the heating system
        Returns
        -------
        None.

        '''

        self.system_type = kwargs["heating_system_key"]
        info_heating = kwargs["heating_system_key"].split(", ")

        self.generation_type = info_heating[0]
        try:
            self.distribution_type = info_heating[1]
            self.emitter_type = info_heating[2]
            self.emission_control_efficiency = systems_info_dict["EN_15316_emission_control_heating_efficiency"]["Efficiency [-]"][self.emitter_type]
            self.distribution_efficiency = systems_info_dict["EN_15316_distribution_heating_efficiency"]["Efficiency [-]"][self.distribution_type]
        except IndexError:
            self.distribution_type = None
            self.emitter_type = None
            self.emission_control_efficiency = 1.
            self.distribution_efficiency = 1.

        self.generation_efficiency_profile = systems_info_dict["EN_15316_generation_heating_efficiency"].loc[self.generation_type]
        self.generation_auxiliary_electric_load_profile = systems_info_dict["EN_15316_generation_heating_auxiliary_electric_load"].loc[self.generation_type]

        # Input Data
        # self.PCI_natural_gas = fuels_pci["Natural Gas"]  # [Wh/Nm3]

    def set_system_capacity(self, design_power, weather):
        '''
        Sets system design power

        Parameters
            ----------
            design_power : float
                Design Heating power  [W]
            Weather : WeatherFile
                WeatherFile object

        Returns
        -------
        None.
        '''

        loads = self.generation_efficiency_profile.index
        generation_efficiency_profile_capacity = np.array([float(load[:-2]) for load in loads])
        loads = self.generation_auxiliary_electric_load_profile.index
        generation_auxiliary_electric_load_profile_capacity = np.array([float(load[:-2]) for load in loads])

        self.generation_efficiency = np.interp(design_power/1000, generation_efficiency_profile_capacity, self.generation_efficiency_profile.values) # [-]
        self.generation_auxiliary_electric_load = np.interp(design_power/1000, generation_auxiliary_electric_load_profile_capacity, self.generation_auxiliary_electric_load_profile.values)*1000 # [W]

        self.total_efficiency = self.emission_control_efficiency * self.distribution_efficiency * self.generation_efficiency

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
        '''
        This method allows to calculate Traditional Boiler losses following
        the Standard UNI-TS 11300:2 - 2008

        Parameters
            ----------
            heat_flow : float
                required power  [W]
            Weather : WeatherFile
                WeatherFile object
            t : int
                Simulation time step
            T_int : float
                Zone temperature [°]
            RH_int : float
                Zone relative humidity [%]

        Returns
        -------
        None.

        '''
        # Corrected efficiency and losses at nominal power

        total_energy = heat_flow / self.total_efficiency

        if "Oil" in self.generation_type:
            self.oil_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Oil"]
            self.electric_consumption = self.generation_auxiliary_electric_load / CONFIG.ts_per_hour
        elif "Stove" in self.generation_type:
            self.wood_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Wood"]
            self.electric_consumption = self.generation_auxiliary_electric_load / CONFIG.ts_per_hour
        elif "Heat Pump" in self.generation_type:
            self.electric_consumption = total_energy / CONFIG.ts_per_hour
        elif "Gas" in self.generation_type:
            self.gas_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Natural Gas"]
            self.electric_consumption = self.generation_auxiliary_electric_load / CONFIG.ts_per_hour

class Cooling_EN15316(System):
    '''
    Class TraditionalBoiler
    This method considers a generic Traditional Boiler as the heating system
    of the entire building following UNI-TS 11300:2 - 2008

    Methods:
        init
        set
        solve
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0

    def __init__(self, *args, **kwargs):
        '''
        init. Set some attributes for the method

        Parameters
            ----------
            system_type : string
                string that sets the heating system
        Returns
        -------
        None.

        '''

        self.system_type = kwargs["cooling_system_key"]
        info_cooling = kwargs["cooling_system_key"].split(", ")

        self.generation_type = info_cooling[0]
        try:
            self.distribution_type = info_cooling[1]
            self.emitter_type = info_cooling[2]
            self.emission_control_efficiency = systems_info_dict["EN_15316_emission_control_cooling_efficiency"]["Efficiency [-]"][self.emitter_type]
            self.distribution_efficiency = systems_info_dict["EN_15316_distribution_cooling_efficiency"]["Efficiency [-]"][self.distribution_type]
        except IndexError:
            self.distribution_type = None
            self.emitter_type = None
            self.emission_control_efficiency = 1.
            self.distribution_efficiency = 1.

        self.generation_seasonal_performance_factor = systems_info_dict["EN_15316_generation_cooling_seasonal_performance_factor"]["SPF [-]"].loc[self.generation_type]
        self.total_efficiency = self.emission_control_efficiency * self.distribution_efficiency * self.generation_seasonal_performance_factor

    def set_system_capacity(self, design_power, weather):
        '''
        Sets system design power

        Parameters
            ----------
            design_power : float
                Design Heating power  [W]
            Weather : WeatherFile
                WeatherFile object

        Returns
        -------
        None.
        '''
        pass

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
        '''
        This method allows to calculate Traditional Boiler losses following
        the Standard UNI-TS 11300:2 - 2008

        Parameters
            ----------
            heat_flow : float
                required power  [W]
            Weather : WeatherFile
                WeatherFile object
            t : int
                Simulation time step
            T_int : float
                Zone temperature [°]
            RH_int : float
                Zone relative humidity [%]

        Returns
        -------
        None.

        '''
        # Corrected efficiency and losses at nominal power

        total_energy = abs(heat_flow) / self.total_efficiency
        self.electric_consumption = total_energy / CONFIG.ts_per_hour


hvac_heating_systems_classes = {
    "IdealLoad":IdealLoad,
    "CondensingBoiler":CondensingBoiler,
    "TraditionalBoiler":TraditionalBoiler,
    "Traditional Gas Boiler, Centralized, Low Temp Radiator":Heating_EN15316,
    "Traditional Gas Boiler, Single, Low Temp Radiator":Heating_EN15316,
    "Traditional Gas Boiler, Centralized, High Temp Radiator":Heating_EN15316,
    "Traditional Gas Boiler, Single, High Temp Radiator":Heating_EN15316,
    "Traditional Gas Boiler, Centralized, Fan coil":Heating_EN15316,
    "Traditional Gas Boiler, Single, Fan coil":Heating_EN15316,
    "Traditional Gas Boiler, Centralized, Radiant surface":Heating_EN15316,
    "Traditional Gas Boiler, Single, Radiant surface":Heating_EN15316,
    "Condensing Gas Boiler, Centralized, Low Temp Radiator":Heating_EN15316,
    "Condensing Gas Boiler, Single, Low Temp Radiator":Heating_EN15316,
    "Condensing Gas Boiler, Centralized, High Temp Radiator":Heating_EN15316,
    "Condensing Gas Boiler, Single, High Temp Radiator":Heating_EN15316,
    "Condensing Gas Boiler, Centralized, Fan coil":Heating_EN15316,
    "Condensing Gas Boiler, Single, Fan coil":Heating_EN15316,
    "Condensing Gas Boiler, Centralized, Radiant surface":Heating_EN15316,
    "Condensing Gas Boiler, Single, Radiant surface":Heating_EN15316,
    "Oil Boiler, Centralized, High Temp Radiator":Heating_EN15316,
    "Oil Boiler, Single, High Temp Radiator":Heating_EN15316,
    "Stove":Heating_EN15316,
    "A-W Heat Pump, Centralized, Low Temp Radiator":Heating_EN15316,
    "A-W Heat Pump, Single, Low Temp Radiator":Heating_EN15316,
    "A-W Heat Pump, Centralized, Fan coil":Heating_EN15316,
    "A-W Heat Pump, Single, Fan coil":Heating_EN15316,
    "A-W Heat Pump, Centralized, Radiant surface":Heating_EN15316,
    "A-W Heat Pump, Single, Radiant surface":Heating_EN15316,
}
hvac_cooling_systems_classes = {
    "IdealLoad":IdealLoad,
    "SplitAirCooler":SplitAirCooler,
    "ChillerAirtoWater":ChillerAirtoWater,
    "SplitAirConditioner":SplitAirConditioner,
    "A-A split":Cooling_EN15316,
    "A-W chiller, Centralized, Fan coil":Cooling_EN15316,
    "A-W chiller, Centralized, Radiant surface":Cooling_EN15316,
    "A-W chiller, Single, Fan coil":Cooling_EN15316,
    "A-W chiller, Single, Radiant surface":Cooling_EN15316,
}