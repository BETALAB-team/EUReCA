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
from eureca_building.fluids_properties import fuels_pci, water_properties
from eureca_building.config import CONFIG
from eureca_building.solar_thermal_system import SolarThermal_Collector

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
    Defines an interface to force inherited class to implement mandatory methods
    """

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    pellet_consumption = 0
    lpg_consumption = 0
    gasoline_consumption = 0

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
    
    # @property
    # @abc.abstractmethod
    # def solar_gain(self):
    #     pass

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
        
    @property
    def solar_thermal_system(self):
        return self._solar_thermal_system

    @solar_thermal_system.setter
    def solar_thermal_system(self, value: SolarThermal_Collector):
        if not isinstance(value, SolarThermal_Collector):
            raise TypeError(
                f"Solar thermal needs to be defined as a SolarThermal_Collector class"
            )
        self._solar_thermal_system = value
        self.set_solar_gain()
        # print(np.max(self.solar_gain))
        
    
    

    def set_dhw_design_capacity_tank(self, dhw_flow_rate, weather_file):

        T_tank = 55
        T_user = 40
        T_tank_max = 80
        T_tank_min= 45
        T_ground = weather_file.general_data['average_out_air_db_temperature']

        pre_heating_time = 2 # [h]

        max_daily_flow_cons = (dhw_flow_rate[:364*24*CONFIG.ts_per_hour].reshape(364, 24*CONFIG.ts_per_hour)*CONFIG.time_step).sum(axis=1).max() # [m3]
        max_flow_rate = dhw_flow_rate.max() * 3600 # m3/h
        max_flow_rate = 1e-10 if max_flow_rate < 1e-10 else max_flow_rate
        request_time = max_daily_flow_cons/max_flow_rate

        # from UNI 9182 [m3]
        self.dhw_tank_volume = 1.5*max_flow_rate * request_time * (T_user - T_ground) / \
                               (pre_heating_time + request_time) * pre_heating_time / (T_tank - T_ground) # [m3]
        self.dhw_tank_volume = np.round(self.dhw_tank_volume,decimals = 1)
        self.dhw_tank_volume  = 1e-10 if self.dhw_tank_volume  < 1e-10 else self.dhw_tank_volume

        self.dhw_design_load = 1.5*max_flow_rate * 1000 * request_time * (T_user - T_ground) / \
                               (pre_heating_time + request_time) * 1.163 # [W]

        # DHW volume
        self.dhw_tank_design_delta_T = T_tank - T_user
        self.dhw_tank_design_charge = self.dhw_tank_volume * \
                               water_properties["density"] * \
                               water_properties["specific_heat"] * \
                               self.dhw_tank_design_delta_T / 3600    # [Wh]
        self.dhw_tank_maximum_charge=self.dhw_tank_volume * \
                               water_properties["density"] * \
                               water_properties["specific_heat"] * \
                               (T_tank_max-T_user) / 3600          
        self.dhw_tank_minimum_charge=self.dhw_tank_volume * \
                               water_properties["density"] * \
                               water_properties["specific_heat"] * \
                               (T_tank_min-T_user) / 3600          
        self.dhw_tank_design_charge  = 1e-10 if self.dhw_tank_design_charge  < 1e-10 else self.dhw_tank_design_charge

        self.dhw_tank_current_charge = self.dhw_tank_design_charge / 2
        self.dhw_tank_current_charge_perc = 0.5
        self.charging_mode = 0.
        # self.discharging_mode =0.
        self.dhw_capacity_to_tank = 0.
        self.dhw_tank_current_charge_perc=self.dhw_tank_minimum_charge/self.dhw_tank_design_charge
        self.losses_discharging_rate = 0.02 # %/h

    def dhw_tank_solver(self, dhw_demand, weather,timestep,**kwargs):
        # if hasattr(self,"solar_gain"):
        #     solar_thermal_gain=self.solar_gain[timestep]/CONFIG.ts_per_hour
        # else:
        #     solar_thermal_gain=0
        self.charging_mode = 0
        self.discharging_mode = 0
        self.solar_gain_out = self.solar_gain.iloc[timestep]/CONFIG.ts_per_hour if hasattr(self,"solar_gain") else 0
        solar_gain=self.solar_gain_out
        self.tank_discharge=0
        self.dhw_capacity_to_tank=0
        loss_rate=self.losses_discharging_rate*max(1,self.dhw_tank_current_charge_perc)
        self.storage_tank_loss=self.dhw_tank_design_charge * loss_rate/CONFIG.ts_per_hour/100

        self.dhw_tank_current_charge=self.dhw_tank_current_charge+solar_gain-dhw_demand/CONFIG.ts_per_hour
        self.dhw_tank_current_charge=self.dhw_tank_current_charge-self.storage_tank_loss
        self.dhw_tank_current_charge_perc = self.dhw_tank_current_charge / self.dhw_tank_design_charge *100
        if (self.dhw_tank_current_charge<self.dhw_tank_minimum_charge):
            self.charging_mode = 1
            self.dhw_capacity_to_tank=min(self.dhw_tank_design_charge-self.dhw_tank_current_charge,self.dhw_design_load/CONFIG.ts_per_hour)
            self.dhw_tank_current_charge=self.dhw_tank_current_charge+self.dhw_capacity_to_tank
        if (self.dhw_tank_current_charge>self.dhw_tank_maximum_charge):
            self.discharging_mode = 1
            # self.charging_mode = 0
            self.tank_discharge=self.dhw_tank_current_charge-self.dhw_tank_maximum_charge
            self.dhw_tank_current_charge=self.dhw_tank_current_charge-self.tank_discharge
            

        # self.dhw_capacity_from_tank =  dhw_demand+ self.dhw_tank_design_charge * self.losses_discharging_rate/CONFIG.ts_per_hour
        # self.dhw_tank_current_charge=min(1.1*self.dhw_tank_design_charge,self.dhw_tank_current_charge+self.dhw_capacity_to_tank-self.dhw_capacity_from_tank)
        self.dhw_tank_current_charge_perc = self.dhw_tank_current_charge / self.dhw_tank_design_charge *100 
        # self.dhw_tank_current_charge_perc = self.dhw_tank_current_charge

    def set_solar_gain(self):
        if hasattr(self,"solar_thermal_system"):
            self.solar_gain = self.solar_thermal_system.gained_heat



# %%---------------------------------------------------------------------------------------------------
# %% IdealLoad class

class IdealLoad(System):
    '''Class IdealLoad is for the ideal  zone balance. This class passes all methods.
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0
    pellet_consumption = 0
    lpg_consumption = 0
    gasoline_consumption = 0

    def __init__(self, *args, **kwargs):
        """IdealLoad init method. No input needed

        Parameters
        ----------
        args
        kwargs
        """

        self.convective_fraction = 0.65
        self.sigma = {
            "1C" : (1-self.convective_fraction, self.convective_fraction),
            "2C" : (
            (1-self.convective_fraction)/2, # Radiative IW
            (1-self.convective_fraction)/2, # Radiative AW
            self.convective_fraction)       # Convective
        }

    def set_system_capacity(self, design_power, weather):
        """Set system capacity method.

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        weather : eureca_building.weather.WeatherFile
            WeatherFile object
        """
        pass

    def solve_system(self, heat_flow, weather, t, T_int, RH_int, *args):
        """Solve the system power for the time step

        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]

        """
        pass


# %%---------------------------------------------------------------------------------------------------
# %% CondensingBoiler class
class CondensingBoiler(System):
    '''Class CondensingBoiler. This method considers a generic Condensing Boiler as the heating system
    of the entire building following UNI-TS 11300:2 - 2008
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0
    pellet_consumption = 0
    lpg_consumption = 0
    gasoline_consumption = 0

    def __init__(self, *args, **kwargs):
        '''init method. Set some attributes for the method are initialized

        Parameters
        ----------
        args
        kwargs
        '''

        # Input Data
        self.system_type = "CondensingBoiler"
        self.theta_gn_test_Pn = 70  # [°C]
        self.theta_gn_test_Pint = 35  # [°C]
        self.theta_a_test = 20  # [°C]
        self.PCI_natural_gas = fuels_pci["Natural Gas"]  # Wh/Nm3
        self.system_info = systems_info_dict["CondensingBoiler"]

        self.convective_fraction = 0.5
        self.sigma = {
            "1C" : (1-self.convective_fraction, self.convective_fraction),
            "2C" : (
            (1-self.convective_fraction)/2, # Radiative IW
            (1-self.convective_fraction)/2, # Radiative AW
            self.convective_fraction)       # Convective
        }

        # if ST exists:
        #     self.ST_system = kwargs["ST system"]
        # else:
        #     self.ST_system = None

        # Vectors initialization
        # self.phi_gn_i_Px = np.zeros(l)  # [kW]
        # self.W_aux_Px = np.zeros(l)  # [kW]
        # self.phi_gn_i_Px_1 = np.zeros(l)  # [kW]
        # self.phi_gn_i_Px_2 = np.zeros(l)  # [kW]

    def set_system_capacity(self, design_power, weather):
        '''Sets system design power

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
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

        # if ST exists:
        #     self.ST_system.calc_production(weather file)

    def solve_system(self, heat_flow, dhw_flow, weather, t, T_int, RH_int):
        '''This method allows to calculate Condensing Boiler power and losses following
         the Standard UNI-TS 11300:2 - 2008

         Parameters
         ----------
         heat_flow : float
             required power  [W]
         Weather : eureca_building.weather.WeatherFile
             WeatherFile object
         t : int
             Simulation time step
         T_int : float
             Zone temperature [°]
         RH_int : float
             Zone relative humidity [%]
         '''

        self.dhw_tank_solver(dhw_flow, weather,t)
        heat_flow += self.dhw_capacity_to_tank

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
    '''Class TraditionalBoiler. This class considers a generic Traditional Boiler as the heating system
    of the entire building following UNI-TS 11300:2 - 2008
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0
    pellet_consumption = 0
    lpg_consumption = 0
    gasoline_consumption = 0

    def __init__(self, *args, **kwargs):
        '''init method. Set some attributes for the method

        Parameters
        ----------
        args
        kwargs
        '''

        # Input Data
        self.system_type = "TraditionalBoiler"
        self.theta_gn_test_Pn = 70  # [°C]
        self.theta_gn_test_Pint = 50  # [°C]
        self.theta_a_test = 20  # [°C]
        self.PCI_natural_gas = fuels_pci["Natural Gas"]  # [Wh/Nm3]
        self.system_info = systems_info_dict["TraditionalBoiler"]

        self.convective_fraction = 0.5
        self.sigma = {
            "1C" : (1-self.convective_fraction, self.convective_fraction),
            "2C" : (
            (1-self.convective_fraction)/2, # Radiative IW
            (1-self.convective_fraction)/2, # Radiative AW
            self.convective_fraction)       # Convective
        }

        self.dhw_tank_volume = 0.2 # [m3]

        # Vectors initialization
        # self.phi_gn_i_Px = np.zeros(l)  # [kW]
        # self.W_aux_Px = np.zeros(l)  # [kW]
        # self.phi_gn_i_Px_1 = np.zeros(l)  # [kW]
        # self.phi_gn_i_Px_2 = np.zeros(l)  # [kW]

    def set_system_capacity(self, design_power, weather):
        '''Sets system design power

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
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

    def solve_system(self, heat_flow, dhw_flow, weather, t, T_int, RH_int):
        '''This method allows to calculate Traditional Boiler losses following
        the Standard UNI-TS 11300:2 - 2008

        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]

        '''

        self.dhw_tank_solver(dhw_flow, weather,t)
        heat_flow += self.dhw_capacity_to_tank

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
    '''SplitAirCooler class. This class considers a generic Split Air Conditioner as the cooling system
    of the entire building following UNI-TS 11300:3 - 2010
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0

    def __init__(self, *args, **kwargs):
        '''init method. Set some attributes for the method

        Parameters
        ----------
        args
        kwargs
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

        self.convective_fraction = 1.
        self.sigma = {
            "1C" : (1-self.convective_fraction, self.convective_fraction),
            "2C" : (
            (1-self.convective_fraction)/2, # Radiative IW
            (1-self.convective_fraction)/2, # Radiative AW
            self.convective_fraction)       # Convective
        }

        # Vectors initialization
        # self.W_el = np.zeros(l)  # [kWe]
        # self.W_aux = np.zeros(l)  # [kWe]
        # self.H_waste = np.zeros(l)  # [kW]
        # self.eta_mm = np.zeros(l)
        # self.eta_1 = np.zeros(l)

    def set_system_capacity(self, design_power, weather):
        '''Choice of system size based on estimated nominal Power

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
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
        '''This method allows to calculate Split Air Cooler electrical power
        following the Standard UNI-TS 11300:3 - 2010

        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]
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
    '''ChillerAirtoWater class. This method considers a generic Air-to-water Chiller as the cooling system
    of the entire building following UNI-TS 11300:3 - 2010
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0

    def __init__(self, *args, **kwargs):
        '''init method. Set some attributes for the method

        Parameters
        ----------
        args
        kwargs
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

        self.convective_fraction = 0.5
        self.sigma = {
            "1C" : (1-self.convective_fraction, self.convective_fraction),
            "2C" : (
            (1-self.convective_fraction)/2, # Radiative IW
            (1-self.convective_fraction)/2, # Radiative AW
            self.convective_fraction)       # Convective
        }

        # Vectors initialization
        # self.W_el = np.zeros(l)  # [kWe]
        # self.W_aux = np.zeros(l)  # [kWe]
        # self.H_waste = np.zeros(l)  # [kW]
        # self.eta_mm = np.zeros(l)
        # self.eta_1 = np.zeros(l)

    def set_system_capacity(self, design_power, weather):
        '''Choice of system size based on estimated nominal Power

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
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
        '''This method allows to calculates Air to Water Chiller electrical power
        following the Standard UNI-TS 11300:3 - 2010


        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]
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
    '''SplitAirConditioner class

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
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0


    def __init__(self, *args, **kwargs):
        '''init method. Set some attributes for the method

        Parameters
        ----------
        args
        kwargs
        '''

        # Check input data type

        # The quality control is implemented in the System class.. Here is not necessary

        # Input Data
        self.system_type = "SplitAirConditioner"  # system_type è una stringa
        self.T_ext_rif = 35  # [°C]
        self.T_int_rif = 27  # [°C]
        self.W_aux_gn = 0.04  # [kW/kW_cond]
        self.system_info = systems_info_dict["SplitAirConditioner"]

        self.convective_fraction = 1.
        self.sigma = {
            "1C" : (1-self.convective_fraction, self.convective_fraction),
            "2C" : (
            (1-self.convective_fraction)/2, # Radiative IW
            (1-self.convective_fraction)/2, # Radiative AW
            self.convective_fraction)       # Convective
        }

        # self.W_el = np.zeros(l)  # [kWe]
        # self.W_aux = np.zeros(l)  # [kWe]
        # self.H_waste = np.zeros(l)  # [kW]
        # self.eta_mm = np.zeros(l)
        # self.EER = np.zeros(l)
        # self.F = np.zeros(l)

    def set_system_capacity(self, design_power, weather):  # da modificare dentro il file excel di impianti
        '''Choice of system size based on estimated nominal Power for air conditioner WITHOUT inverter

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
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
        '''This method allows to calculate from literature: Split Air Cooler condensing power

        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]
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
    '''Class Heating_EN15316. This method considers a generic heating system as the heating system
    of the entire building following EN 15316.
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0
    pellet_consumption = 0
    lpg_consumption = 0
    gasoline_consumption = 0

    def __init__(self, *args, **kwargs):
        '''init method. Set some attributes are set
        The heating_system_key label must be passed as kwargs. Example:
        Heating_EN15316(heating_system_key = "Traditional Gas Boiler, Single, Low Temp Radiator")

        Parameters
        ----------
        args : list
            list of optional arguments
        kwargs : dict
            kwargs must include {heating_system_key : string_of_heating_system}

        '''

        self.system_type = kwargs["heating_system_key"]
        info_heating = kwargs["heating_system_key"].split(", ")

        self.generation_type = info_heating[0]
        try:
            self.distribution_type = info_heating[1]
            self.emitter_type = info_heating[2]
            self.emission_control_efficiency = systems_info_dict["EN_15316_emission_control_heating_efficiency"]["Efficiency [-]"][self.emitter_type]
            self.distribution_efficiency = systems_info_dict["EN_15316_distribution_heating_efficiency"]["Efficiency [-]"][self.distribution_type]
            self.convective_fraction = systems_info_dict["EN_15316_emission_control_heating_efficiency"]["Convective fraction [-]"][self.emitter_type]
        except IndexError:
            # Stove
            self.distribution_type = None
            self.emitter_type = None
            self.emission_control_efficiency = 1.
            self.distribution_efficiency = 1.
            self.convective_fraction = 0.5

        self.generation_efficiency_profile = systems_info_dict["EN_15316_generation_heating_efficiency"].loc[self.generation_type]
        self.generation_auxiliary_electric_load_profile = systems_info_dict["EN_15316_generation_heating_auxiliary_electric_load"].loc[self.generation_type]

        self.sigma = {
            "1C" : (1-self.convective_fraction, self.convective_fraction),
            "2C" : (
            (1-self.convective_fraction)/2, # Radiative IW
            (1-self.convective_fraction)/2, # Radiative AW
            self.convective_fraction)       # Convective
        }

        self.dhw_tank_volume = 0.2 # [m3]

        # Input Data
        # self.PCI_natural_gas = fuels_pci["Natural Gas"]  # [Wh/Nm3]

    def set_system_capacity(self, design_power, weather):
        ''''Choice of system size based on estimated nominal Power

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        '''

        loads = self.generation_efficiency_profile.index
        generation_efficiency_profile_capacity = np.array([float(load[:-2]) for load in loads])
        loads = self.generation_auxiliary_electric_load_profile.index
        generation_auxiliary_electric_load_profile_capacity = np.array([float(load[:-2]) for load in loads])

        self.generation_efficiency = np.interp(design_power/1000, generation_efficiency_profile_capacity, self.generation_efficiency_profile.values) # [-]
        self.generation_auxiliary_electric_load = np.interp(design_power/1000, generation_auxiliary_electric_load_profile_capacity, self.generation_auxiliary_electric_load_profile.values)*1000 # [W]

        self.total_efficiency = self.emission_control_efficiency * self.distribution_efficiency * self.generation_efficiency

    def solve_system(self, heat_flow, dhw_flow, weather, t, T_int, RH_int,**kwargs):
        '''This method allows to calculate the system power for each time step

        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]
        kwargs
        '''


        self.dhw_tank_solver(dhw_flow, weather,t)
        heat_flow += self.dhw_capacity_to_tank

        total_energy = heat_flow / self.total_efficiency

        if "Oil" in self.generation_type:
            self.oil_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Oil"]
            self.electric_consumption = self.generation_auxiliary_electric_load / CONFIG.ts_per_hour
        elif "Coal" in self.generation_type:
            self.coal_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Coal"]
            self.electric_consumption = self.generation_auxiliary_electric_load / CONFIG.ts_per_hour
        elif "District Heating" in self.generation_type:
            self.DH_consumption = total_energy / CONFIG.ts_per_hour
            self.electric_consumption = self.generation_auxiliary_electric_load / CONFIG.ts_per_hour
        elif "Stove" in self.generation_type:
            self.wood_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Wood"]
            self.electric_consumption = self.generation_auxiliary_electric_load / CONFIG.ts_per_hour
        elif "Heat Pump" in self.generation_type:
            self.electric_consumption = total_energy / CONFIG.ts_per_hour
        elif "Gas" in self.generation_type:
            self.gas_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Natural Gas"]
            self.electric_consumption = self.generation_auxiliary_electric_load / CONFIG.ts_per_hour
        elif "Electric Heater" in self.generation_type:
            self.electric_consumption = total_energy / CONFIG.ts_per_hour + self.generation_auxiliary_electric_load / CONFIG.ts_per_hour

    def solve_quasi_steady_state(self, heat_flow, dhw_flow):
        '''This method allows to calculate the system power for each time step

        Parameters
        ----------
        heat_flow : float
            required power  [Wh]
        dhw_flow : float
            required power  [Wh]
        '''
        # Corrected efficiency and losses at nominal power

        self.oil_consumption = 0
        self.coal_consumption = 0
        self.DH_consumption = 0
        self.wood_consumption = 0
        self.pellet_consumption = 0
        self.electric_consumption = 0
        self.gas_consumption = 0
        self.gasoline_consumption = 0
        self.lpg_consumption = 0

        total_energy = (heat_flow + dhw_flow) / self.total_efficiency # Wh

        if "Oil" in self.generation_type:
            self.oil_consumption = total_energy / fuels_pci["Oil"]
        elif "Coal" in self.generation_type:
            self.coal_consumption = total_energy / fuels_pci["Coal"]
        elif "District Heating" in self.generation_type:
            self.DH_consumption = total_energy
        elif "Stove" in self.generation_type:
            self.wood_consumption = total_energy / fuels_pci["Wood"]
        elif "Heat Pump" in self.generation_type:
            self.electric_consumption = total_energy
        elif "Gas" in self.generation_type:
            self.gas_consumption = total_energy / fuels_pci["Natural Gas"]
        elif "Electric Heater" in self.generation_type:
            self.electric_consumption = total_energy


class Cooling_EN15316(System):
    '''Class Cooling_EN15316. This method considers a generic cooling system as the heating system
    of the entire building following EN 15316.
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0

    def __init__(self, *args, **kwargs):
        '''init method. Set some attributes are set
        The cooling_system_key label must be passed as kwargs. Example:
        Cooling_EN15316(cooling_system_key = "A-W chiller, Centralized, Radiant surface")

        Parameters
        ----------
        args
        kwargs
            kwargs must include {cooling_system_key : string_of_cooling_system}
        '''

        self.system_type = kwargs["cooling_system_key"]
        info_cooling = kwargs["cooling_system_key"].split(", ")

        self.generation_type = info_cooling[0]
        try:
            self.distribution_type = info_cooling[1]
            self.emitter_type = info_cooling[2]
            self.emission_control_efficiency = systems_info_dict["EN_15316_emission_control_cooling_efficiency"]["Efficiency [-]"][self.emitter_type]
            self.distribution_efficiency = systems_info_dict["EN_15316_distribution_cooling_efficiency"]["Efficiency [-]"][self.distribution_type]
            self.convective_fraction = systems_info_dict["EN_15316_emission_control_cooling_efficiency"]["Convective fraction [-]"][self.emitter_type]
        except IndexError:
            self.distribution_type = None
            self.emitter_type = None
            self.emission_control_efficiency = 1.
            self.distribution_efficiency = 1.
            self.convective_fraction = 1.

        self.generation_seasonal_performance_factor = systems_info_dict["EN_15316_generation_cooling_seasonal_performance_factor"]["SPF [-]"].loc[self.generation_type]
        self.total_efficiency = self.emission_control_efficiency * self.distribution_efficiency * self.generation_seasonal_performance_factor

        self.sigma = {
            "1C" : (1-self.convective_fraction, self.convective_fraction),
            "2C" : (
            (1-self.convective_fraction)/2, # Radiative IW
            (1-self.convective_fraction)/2, # Radiative AW
            self.convective_fraction)       # Convective
        }

    def set_system_capacity(self, design_power, weather):
        ''''Choice of system size based on estimated nominal Power

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        '''
        pass

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
        '''This method allows to calculate the system power for each time step

        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]
        '''
        # Corrected efficiency and losses at nominal power

        total_energy = abs(heat_flow) / self.total_efficiency
        self.electric_consumption = total_energy / CONFIG.ts_per_hour

    def solve_quasi_steady_state(self, heat_flow):
        '''This method allows to calculate the system power for each month

        Parameters
        ----------
        heat_flow : float
            required power  [Wh]
        '''
        # Corrected efficiency and losses at nominal power

        total_energy = abs(heat_flow) / self.total_efficiency

        self.electric_consumption = total_energy

class HP_Staffell(System):
    '''Class Heating HP Staffel.
    '''

    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0
    pellet_consumption = 0
    lpg_consumption = 0
    gasoline_consumption = 0

    def __init__(self, *args, **kwargs):
        '''init method. Set some attributes are set
        The heating_system_key label must be passed as kwargs. Example:
        Heating_EN15316(heating_system_key = "Traditional Gas Boiler, Single, Low Temp Radiator")

        Parameters
        ----------
        args : list
            list of optional arguments
        kwargs : dict
            kwargs must include {heating_system_key : string_of_heating_system}

        '''

        self.system_type = kwargs["heating_system_key"]
        info_heating = kwargs["heating_system_key"].split(", ")

        self.generation_type = info_heating[0]
        try:
            self.distribution_type = info_heating[1]
            self.emitter_type = info_heating[2]
            self.emission_control_efficiency = systems_info_dict["EN_15316_emission_control_heating_efficiency"]["Efficiency [-]"][self.emitter_type]
            self.distribution_efficiency = systems_info_dict["EN_15316_distribution_heating_efficiency"]["Efficiency [-]"][self.distribution_type]
            self.convective_fraction = systems_info_dict["EN_15316_emission_control_heating_efficiency"]["Convective fraction [-]"][self.emitter_type]
        except IndexError:
            # Stove
            self.distribution_type = None
            self.emitter_type = None
            self.emission_control_efficiency = 1.
            self.distribution_efficiency = 1.
            self.convective_fraction = 0.5

        self.T_emitter = {
            "Low Temp Radiator":60.,
            "High Temp Radiator":70.,
            "Radiant surface":35.,
            "Fan coil":40,
        }[self.emitter_type]

        self.COP_fun = {
            "A-W HP Staffel":lambda t_ext : 6.81 - 0.121 * (self.T_emitter - t_ext) + 0.000630 * (self.T_emitter - t_ext)**2,
            "G-W HP Staffel":lambda t_g : 8.77 - 0.150 * (self.T_emitter - t_g) + 0.000734 * (self.T_emitter - t_g)**2,
        }[self.generation_type]

        self.sigma = {
            "1C" : (1-self.convective_fraction, self.convective_fraction),
            "2C" : (
            (1-self.convective_fraction)/2, # Radiative IW
            (1-self.convective_fraction)/2, # Radiative AW
            self.convective_fraction)       # Convective
        }

        self.dhw_tank_volume = 0.2 # [m3]

        # Input Data
        # self.PCI_natural_gas = fuels_pci["Natural Gas"]  # [Wh/Nm3]

    def set_system_capacity(self, design_power, weather):
        self.total_efficiency = self.emission_control_efficiency * self.distribution_efficiency

    def solve_system(self, heat_flow, dhw_flow, weather, t, T_int, RH_int,**kwargs):
        '''This method allows to calculate the system power for each time step

        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]
        kwargs
        '''

        self.dhw_tank_solver(dhw_flow, weather, t)
        heat_flow += self.dhw_capacity_to_tank

        T_source = {
            "A-W HP Staffel": weather.hourly_data["out_air_db_temperature"][t],
            "G-W HP Staffel": weather.general_data['average_out_air_db_temperature'] - 2,
        }[self.generation_type]

        COP_ts = self.COP_fun(T_source)

        total_energy = heat_flow / self.total_efficiency / COP_ts

        self.electric_consumption = total_energy / CONFIG.ts_per_hour

    # def solve_quasi_steady_state(self, heat_flow, dhw_flow):
    #     '''This method allows to calculate the system power for each time step
    #
    #     Parameters
    #     ----------
    #     heat_flow : float
    #         required power  [Wh]
    #     dhw_flow : float
    #         required power  [Wh]
    #     '''
    #     # Corrected efficiency and losses at nominal power
    #
    #     self.oil_consumption = 0
    #     self.coal_consumption = 0
    #     self.DH_consumption = 0
    #     self.wood_consumption = 0
    #     self.pellet_consumption = 0
    #     self.electric_consumption = 0
    #     self.gas_consumption = 0
    #     self.gasoline_consumption = 0
    #     self.lpg_consumption = 0
    #
    #     total_energy = (heat_flow + dhw_flow) / self.total_efficiency # Wh
    #
    #     if "Oil" in self.generation_type:
    #         self.oil_consumption = total_energy / fuels_pci["Oil"]
    #     elif "Coal" in self.generation_type:
    #         self.coal_consumption = total_energy / fuels_pci["Coal"]
    #     elif "District Heating" in self.generation_type:
    #         self.DH_consumption = total_energy
    #     elif "Stove" in self.generation_type:
    #         self.wood_consumption = total_energy / fuels_pci["Wood"]
    #     elif "Heat Pump" in self.generation_type:
    #         self.electric_consumption = total_energy
    #     elif "Gas" in self.generation_type:
    #         self.gas_consumption = total_energy / fuels_pci["Natural Gas"]
    #     elif "Electric Heater" in self.generation_type:
    #         self.electric_consumption = total_energy


class HeatingFromParams(System):
    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0
    pellet_consumption = 0

    def __init__(self, *args, **kwargs):

        # kwargs =
        #       heating_system_params = {
        #     "name":...,
        #     "description":....,
        #
        #     "SH emission system": ...,
        #     "SH emission target temperature [°C]": ...,
        #     "SH emission convective fraction [-]": ...,
        #     "SH emission efficiency [-]": ...,
        #     "SH distribution efficiency [-]": ...,
        #     "SH regulation efficiency [-]": ...,
        #     "SH generation efficiency [-]": ...,
        #     "SH COP [-]": ...,
        #     "SH fuel": ...,
        #
        #     "DHW emission efficiency [-]": ...,
        #     "DHW distribution efficiency [-]": ...,
        #     "DHW regulation efficiency [-]": ...,
        #     "DHW generation efficiency [-]": ...,
        #     "DHW COP [-]": ...,
        #     "DHW fuel": ...,
        # }

        try:
            self.system_info = kwargs["heating_system_params"]
        except KeyError:
            raise KeyError(
                "When using HeatingFromParams system class, a heating_system_params dict must be provided, with heating system performance parameters")

        self.name = self.system_info["name"]
        self.description = self.system_info["description"]

        self.emission_control_efficiency = self.system_info["SH emission efficiency [-]"] * self.system_info[
            "SH regulation efficiency [-]"]
        self.emission_temp = self.system_info["SH emission target temperature [°C]"]
        self.distribution_efficiency = self.system_info["SH distribution efficiency [-]"]
        self.convective_fraction = self.system_info["SH emission convective fraction [-]"]
        self.generation_efficiency = self.system_info["SH generation efficiency [-]"]

        self.sigma = {
            "1C": (1 - self.convective_fraction, self.convective_fraction),
            "2C": (
                (1 - self.convective_fraction) / 2,  # Radiative IW
                (1 - self.convective_fraction) / 2,  # Radiative AW
                self.convective_fraction)  # Convective
        }

        self.dhw_emission_control_efficiency = self.system_info["DHW emission efficiency [-]"] * self.system_info[
            "DHW regulation efficiency [-]"]
        self.dhw_distribution_efficiency = self.system_info["DHW distribution efficiency [-]"]
        self.dhw_generation_efficiency = self.system_info["DHW generation efficiency [-]"]

        self.fuel_type = self.system_info["SH fuel"]
        self.dhw_fuel_type = self.system_info["DHW fuel"]

        # In case of HP needs to be changed
        self.total_efficiency = self.emission_control_efficiency * self.distribution_efficiency * self.generation_efficiency
        self.dhw_total_efficiency = self.dhw_emission_control_efficiency * self.dhw_distribution_efficiency * self.dhw_generation_efficiency

        # Input Data
        # self.PCI_natural_gas = fuels_pci["Natural Gas"]  # [Wh/Nm3]

        self.charging_mode = np.nan
        self.dhw_tank_current_charge_perc = np.nan
        self.dhw_capacity_to_tank = np.nan

        if "Electric" in self.fuel_type:
            self.COP = 1 if np.isnan(self.system_info["SH COP [-]"]) else self.system_info["SH COP [-]"]
        if "Electric" in self.dhw_fuel_type:
            self.dhw_COP = 1 if np.isnan(self.system_info["DHW COP [-]"]) else self.system_info["DHW COP [-]"]

    def set_system_capacity(self, design_power, weather):
        ''''Choice of system size based on estimated nominal Power

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        '''

        pass

    def solve_system(self, heat_flow, dhw_flow, weather, t, T_int, RH_int):
        '''This method allows to calculate the system power for each time step

        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]
        '''
        # Corrected efficiency and losses at nominal power

        self.oil_consumption = 0
        self.coal_consumption = 0
        self.DH_consumption = 0
        self.wood_consumption = 0
        self.pellet_consumption = 0
        self.electric_consumption = 0
        self.gas_consumption = 0
        self.gasoline_consumption = 0
        self.lpg_consumption = 0

        self.dhw_tank_solver(dhw_flow, weather, t)

        total_energy = heat_flow / self.total_efficiency
        dhw_total_energy = self.dhw_capacity_to_tank / self.dhw_total_efficiency

        if "Oil" in self.fuel_type:
            self.oil_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Oil"]
        elif "Gasoline" in self.fuel_type:
            self.gasoline_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Gasoline"]
        elif "Coal" in self.fuel_type:
            self.coal_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Coal"]
        elif "LPG" in self.fuel_type:
            self.lpg_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["LPG"]
        elif "Natural" in self.fuel_type and "Gas" in self.fuel_type :
            self.gas_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Natural Gas"]
        elif "Wood" in self.fuel_type:
            self.wood_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Wood"]
        elif "Pellet" in self.fuel_type:
            self.pellet_consumption = total_energy / CONFIG.ts_per_hour / fuels_pci["Pellets"]
        elif "Electric" in self.fuel_type:
            self.electric_consumption = total_energy / CONFIG.ts_per_hour / self.COP
        elif "DH" in self.fuel_type:
            self.DH_consumption = total_energy / CONFIG.ts_per_hour

        if "Oil" in self.dhw_fuel_type:
            self.oil_consumption += dhw_total_energy / CONFIG.ts_per_hour / fuels_pci["Oil"]
        elif "Gasoline" in self.dhw_fuel_type:
            self.gasoline_consumption += dhw_total_energy / CONFIG.ts_per_hour / fuels_pci["Gasoline"]
        elif "Coal" in self.dhw_fuel_type:
            self.coal_consumption += dhw_total_energy / CONFIG.ts_per_hour / fuels_pci["Coal"]
        elif "LPG" in self.dhw_fuel_type:
            self.lpg_consumption += dhw_total_energy / CONFIG.ts_per_hour / fuels_pci["LPG"]
        elif "Natural" in self.dhw_fuel_type and "Gas" in self.dhw_fuel_type:
            self.gas_consumption += dhw_total_energy / CONFIG.ts_per_hour / fuels_pci["Natural Gas"]
        elif "Wood" in self.dhw_fuel_type:
            self.wood_consumption += dhw_total_energy / CONFIG.ts_per_hour / fuels_pci["Wood"]
        elif "Wood" in self.dhw_fuel_type:
            self.pellet_consumption += dhw_total_energy / CONFIG.ts_per_hour / fuels_pci["Pellets"]
        elif "Electric" in self.dhw_fuel_type:
            self.electric_consumption += dhw_total_energy / CONFIG.ts_per_hour / self.dhw_COP
        elif "DH" in self.dhw_fuel_type:
            self.DH_consumption += dhw_total_energy / CONFIG.ts_per_hour

    def solve_quasi_steady_state(self, heat_flow, dhw_flow):
        '''This method allows to calculate the system power for each month

        Parameters
        ----------
        heat_flow : float
            required power  [Wh]
        dhw_flow : float
            required power  [Wh]
        '''
        # Corrected efficiency and losses at nominal power

        self.oil_consumption = 0
        self.coal_consumption = 0
        self.DH_consumption = 0
        self.wood_consumption = 0
        self.pellet_consumption = 0
        self.electric_consumption = 0
        self.gas_consumption = 0
        self.gasoline_consumption = 0
        self.lpg_consumption = 0

        total_energy = heat_flow / self.total_efficiency  # Wh

        if "Oil" in self.fuel_type:
            self.oil_consumption = total_energy / fuels_pci["Oil"]
        elif "Gasoline" in self.fuel_type:
            self.gasoline_consumption = total_energy / fuels_pci["Gasoline"]
        elif "Coal" in self.fuel_type:
            self.coal_consumption = total_energy / fuels_pci["Coal"]
        elif "LPG" in self.fuel_type:
            self.lpg_consumption = total_energy / fuels_pci["LPG"]
        elif "NaturalGas" in self.fuel_type:
            self.gas_consumption = total_energy / fuels_pci["Natural Gas"]
        elif "Wood" in self.fuel_type:
            self.wood_consumption = total_energy / fuels_pci["Wood"]
        elif "Pellet" in self.fuel_type:
            self.pellet_consumption = total_energy / fuels_pci["Pellets"]
        elif "Electric" in self.fuel_type:
            self.electric_consumption = total_energy / self.COP
        elif "DH" in self.fuel_type:
            self.DH_consumption = total_energy / CONFIG.ts_per_hour

        dhw_total_energy = dhw_flow / self.dhw_total_efficiency

        if "Oil" in self.dhw_fuel_type:
            self.oil_consumption += dhw_total_energy / fuels_pci["Oil"]
        elif "Gasoline" in self.dhw_fuel_type:
            self.gasoline_consumption += dhw_total_energy / fuels_pci["Gasoline"]
        elif "Coal" in self.dhw_fuel_type:
            self.coal_consumption += dhw_total_energy / fuels_pci["Coal"]
        elif "LPG" in self.dhw_fuel_type:
            self.lpg_consumption += dhw_total_energy / fuels_pci["LPG"]
        elif "NaturalGas" in self.dhw_fuel_type:
            self.gas_consumption += dhw_total_energy / fuels_pci["Natural Gas"]
        elif "Wood" in self.dhw_fuel_type:
            self.wood_consumption += dhw_total_energy / fuels_pci["Wood"]
        elif "Pellet" in self.dhw_fuel_type:
            self.pellet_consumption += dhw_total_energy / fuels_pci["Pellets"]
        elif "Electric" in self.dhw_fuel_type:
            self.electric_consumption += dhw_total_energy / self.dhw_COP
        elif "DH" in self.dhw_fuel_type:
            self.DH_consumption += dhw_total_energy / CONFIG.ts_per_hour


class CoolingFromParams(System):
    gas_consumption = 0
    electric_consumption = 0
    wood_consumption = 0
    oil_consumption = 0
    coal_consumption = 0
    DH_consumption = 0
    pellet_consumption = 0

    def __init__(self, *args, **kwargs):

        # kwargs =
        #       cooling_system_params = {
        #     "name":...,
        #     "description":....,
        #
        #     "SC emission system": ...,
        #     "SC emission target temperature [°C]": ...,
        #     "SC emission convective fraction [-]": ...,
        #     "SC emission efficiency [-]": ...,
        #     "SC distribution efficiency [-]": ...,
        #     "SC regulation efficiency [-]": ...,
        #     "SC generation efficiency [-]": ...,
        #     "SC fuel": ...,
        # }

        # try:
        self.system_info = kwargs["cooling_system_params"]
        # except KeyError:
        #     raise KeyError("When using CoolingFromParams system class, a cooling_system_params dict must be provided, with cooling system performance parameters")

        self.emission_temp = self.system_info["SC emission target temperature [°C]"]

        self.convective_fraction = self.system_info["SC emission convective fraction [-]"]

        self.emission_control_efficiency = self.system_info["SC emission efficiency [-]"] * self.system_info[
            "SC regulation efficiency [-]"]
        self.distribution_efficiency = self.system_info["SC distribution efficiency [-]"]
        self.convective_fraction = self.system_info["SC emission convective fraction [-]"]
        self.generation_efficiency = self.system_info["SC generation efficiency [-]"]

        self.total_efficiency = self.emission_control_efficiency * self.distribution_efficiency * self.generation_efficiency

        self.sigma = {
            "1C": (1 - self.convective_fraction, self.convective_fraction),
            "2C": (
                (1 - self.convective_fraction) / 2,  # Radiative IW
                (1 - self.convective_fraction) / 2,  # Radiative AW
                self.convective_fraction)  # Convective
        }

        self.fuel_type = self.system_info["SC fuel"]

    def set_system_capacity(self, design_power, weather):
        ''''Choice of system size based on estimated nominal Power

        Parameters
        ----------
        design_power : float
            Design Heating power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        '''

        pass

    def solve_system(self, heat_flow, weather, t, T_int, RH_int):
        '''This method allows to calculate the system power for each time step

        Parameters
        ----------
        heat_flow : float
            required power  [W]
        Weather : eureca_building.weather.WeatherFile
            WeatherFile object
        t : int
            Simulation time step
        T_int : float
            Zone temperature [°]
        RH_int : float
            Zone relative humidity [%]
        '''
        # Corrected efficiency and losses at nominal power

        self.oil_consumption = 0
        self.coal_consumption = 0
        self.DH_consumption = 0
        self.wood_consumption = 0
        self.pellet_consumption = 0
        self.electric_consumption = 0
        self.gas_consumption = 0
        self.gasoline_consumption = 0
        self.lpg_consumption = 0

        total_energy = heat_flow / self.total_efficiency

        if "Electric" in self.fuel_type:
            self.electric_consumption = -1 * total_energy / CONFIG.ts_per_hour
        elif "DH" in self.fuel_type:
            self.DH_consumption = total_energy / CONFIG.ts_per_hour

    def solve_quasi_steady_state(self, heat_flow):
        '''This method allows to calculate the system power for each month

        Parameters
        ----------
        heat_flow : float
            required power  [Wh]
        '''
        # Corrected efficiency and losses at nominal power

        self.oil_consumption = 0
        self.coal_consumption = 0
        self.DH_consumption = 0
        self.wood_consumption = 0
        self.pellet_consumption = 0
        self.electric_consumption = 0
        self.gas_consumption = 0
        self.gasoline_consumption = 0
        self.lpg_consumption = 0

        total_energy = heat_flow / self.total_efficiency  # Wh

        if "Electric" in self.fuel_type:
            self.electric_consumption = -1 * total_energy
        elif "DH" in self.fuel_type:
            self.DH_consumption = total_energy

hvac_heating_systems_classes = {
    "IdealLoad":IdealLoad,
    "CondensingBoiler":CondensingBoiler,
    "TraditionalBoiler":TraditionalBoiler,
    "A-W HP Staffel, Centralized, Low Temp Radiator":HP_Staffell,
    "G-W HP Staffel, Centralized, Low Temp Radiator":HP_Staffell,
    "A-W HP Staffel, Centralized, High Temp Radiator":HP_Staffell,
    "G-W HP Staffel, Centralized, High Temp Radiator":HP_Staffell,
    "A-W HP Staffel, Centralized, Fan coil":HP_Staffell,
    "G-W HP Staffel, Centralized, Fan coil":HP_Staffell,
    "A-W HP Staffel, Centralized, Radiant surface":HP_Staffell,
    "G-W HP Staffel, Centralized, Radiant surface":HP_Staffell,
    "A-W HP Staffel, Single, Low Temp Radiator":HP_Staffell,
    "G-W HP Staffel, Single, Low Temp Radiator":HP_Staffell,
    "A-W HP Staffel, Single, High Temp Radiator":HP_Staffell,
    "G-W HP Staffel, Single, High Temp Radiator":HP_Staffell,
    "A-W HP Staffel, Single, Fan coil":HP_Staffell,
    "G-W HP Staffel, Single, Fan coil":HP_Staffell,
    "A-W HP Staffel, Single, Radiant surface":HP_Staffell,
    "G-W HP Staffel, Single, Radiant surface":HP_Staffell,
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
    "Coal Heater, Centralized, High Temp Radiator":Heating_EN15316,
    "Coal Heater, Single, High Temp Radiator":Heating_EN15316,
    "District Heating, Centralized, High Temp Radiator":Heating_EN15316,
    "District Heating, Centralized, Low Temp Radiator":Heating_EN15316,
    "District Heating, Centralized, Fan coil":Heating_EN15316,
    "District Heating, Centralized, Radiant surface":Heating_EN15316,
    "Stove":Heating_EN15316,
    "A-W Heat Pump, Centralized, Low Temp Radiator":Heating_EN15316,
    "A-W Heat Pump, Single, Low Temp Radiator":Heating_EN15316,
    "A-W Heat Pump, Centralized, Fan coil":Heating_EN15316,
    "A-W Heat Pump, Single, Fan coil":Heating_EN15316,
    "A-W Heat Pump, Centralized, Radiant surface":Heating_EN15316,
    "A-W Heat Pump, Single, Radiant surface":Heating_EN15316,
    "Electric Heater":Heating_EN15316,
    "From manual parameters":HeatingFromParams,
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
    "From manual parameters":CoolingFromParams,
}