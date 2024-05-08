
"""
This module includes the class to manage the Air Handling Unit
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


'''IMPORTING MODULES'''


import sys
import logging
import copy
import numpy as np

from eureca_building.fluids_properties import air_properties, vapour_properties
from eureca_building.schedule import Schedule
from eureca_building.schedule_properties import ventilation_prop
from eureca_building.ventilation import MechanicalVentilation
from eureca_building.weather import WeatherFile
from eureca_building.exceptions import (
    ScheduleOutsideBoundaryCondition,
    PropertyOutsideBoundaries
)

class _BaseAirHandlingUnit:
    '''This class is a class to be inherited by different Mechanical Ventilation Unit objects.
    Some generic methods are stored
    '''

    # Class Variables
    cp_air = air_properties["specific_heat"]         # [J/(kg K)]
    rho_air = air_properties["density"]              # [kg/m3]
    p_atm = air_properties["atmospheric_pressure"]   # [Pa]
    r_0 = vapour_properties["latent_heat"]           # [J/kg]
    cpv = vapour_properties["specific_heat"]         # [J/(kg K)]

    def checkSatCond(self, temp, x, p):
        '''
        Check Saturation Condition

        This function takes as inputs temperature [°C] and humidity ratio
        [kg_vap/kg_as] to check if a point is outside saturation conditions

        Parameters
        ----------
        temp : float
            Temperature [°C]
        x : float
            Specific Humidity [kg_vap/kg_as].
        p : float
            Pressure [Pa].

        Returns
        -------
        tuple
            boolean (wheter saturation is reached), and Saturation Pressure [Pa].

        '''

        # Check input data type

        if not isinstance(temp, float):
            raise TypeError(f'ERROR input T is not an interger: T {temp}')
        if not isinstance(x, float):
            raise TypeError(f'ERROR input x is not an interger: x {x}')

        # Control input data quality
        # if temp < -40 or temp > 70:
        #     logging.warning(
        #         f"WARNING CheckSatCond function, input temperature outside limit boundary [-15,60]: T {temp}"
        #     )
        # if x < 0.0005 or x > 0.040:
        #     logging.warning(f"WARNING CheckSatCond function, input humidity outside limit boundary [0.0005,0.04]: x {x}")

        # Is or not outside the saturation condition? True/False
        pp = p * x / (0.622 + x)
        if temp < 0:
            psat = 610.5 * np.exp((21.875 * temp) / (265.5 + temp))
        else:
            psat = 610.5 * np.exp((17.269 * temp) / (237.3 + temp))
        if pp - psat > 0.01:
            sat_cond = False
        else:
            sat_cond = True
        return sat_cond, psat

    def electric_consumption_based_on_mass_flow_rate(self, mass_flow):
        """
        ######################## Warning #############################
        Typical electric consumption from catalogue data.
        This method calculates a very rough estimation
        of the electric consumption based on the mass flow rate.
        It does not account for real pressure drops or other phenomena

        Parameters
        ----------
        mass_flow : numpy.array
            ventialation mass flow rate [kg/s]

        Returns
        -------
        np.array
            consumption [W].
        """

        vol_flow_m3_h = mass_flow / self.rho_air * 3600
        # consumption_W = -3E-05 * vol_flow_m3_h ** 2 + 0.5251 * vol_flow_m3_h + 8.0258
        # consumption_W = 0.4176*vol_flow_m3_h + 54.377
        consumption_W = 0.4434 * vol_flow_m3_h
        return consumption_W

    def properties(self):
            """ Just a function to print the memorized conditions
            """
            return f"""
    HR :\tT {self._chart_T_hr:.1f} °C,\tx {self._chart_x_hr:.5f} kg/kg,\th {self.h_hr:.1f} J/kg
    MIX:\tT {self._chart_T_mix:.1f} °C,\tx {self._chart_x_mix:.5f} kg/kg,\th {self.h_mix:.1f} J/kg
    PRE:\tT {self._chart_T_preh_deu:.1f} °C,\tx {self._chart_x_preh_deu:.5f} kg/kg,\th {self.h_ph:.1f} J/kg
    AS :\tT {self._chart_T_as:.1f} °C,\tx {self._chart_x_as:.5f} kg/kg,\th {self.h_as:.1f} J/kg
    SUP:\tT {self._chart_T_posth:.1f} °C,\tx {self._chart_x_posth:.5f} kg/kg,\th {self.h_sup:.1f} J/kg

    AHU_SENS:\t{self.AHU_demand_sens} W
    AHU_LAT:\t{self.AHU_demand_lat} W
    AHU_TOT:\t{self.AHU_demand} W
            """

        # %%---------------------------------------------------------------------------------------------------
        # %%

    def _psychro_plot(self):
        """ Just a function to get a psychrometric chart (internal use only)
        """
        try:
            import matplotlib.pyplot as plt
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                "To run the AirHandlingUnit._psychro_plot you need to install matplotlib package")
        fig, ax = plt.subplots()
        ax.set_ylabel("Specific Humidity [" + "$g_{v}/kg_{da}$" + "]")
        ax.set_xlabel("Temperature [" + "$°C$" + "]")
        ax.set_xlim(-10, 45)
        ax.set_ylim(0, 30)
        ax.set_title("Psychrometric chart")
        t = np.arange(-10, 46, 1)
        p_sat = 6.1094 * np.exp(17.625 * t / (t + 243.04)) * 100
        for ur in np.arange(0.1, 1.1, 0.1):
            p = p_sat * ur
            sh = 0.622 * (p / (101325 - p)) * 1000
            ax.plot(t, sh, 'k-', linewidth=0.3)
            x_text = min([t[-1] - 5, 35])
            y_text = min([sh[-1] - 2, 28])

        ax.text(27.7, 26.9, f"100%", backgroundcolor="white", fontsize=6, ma="center")
        ax.text(30.5, 23.7, f"80%", backgroundcolor="white", fontsize=6, ma="center")
        ax.text(32.5, 20., f"60%", backgroundcolor="white", fontsize=6, ma="center")
        ax.text(34.5, 14.5, f"40%", backgroundcolor="white", fontsize=6, ma="center")
        ax.text(35.8, 8, f"20%", backgroundcolor="white", fontsize=6, ma="center")

        x_text, y_text = -9., 3.5

        for h in np.arange(0., 200., 10.):
            x = (h - 1.006 * t) / (1.86 * t + 2501) * 1000
            ax.plot(t, x, 'k:', linewidth=0.3)
            if y_text < 30.:
                ax.text(x_text, y_text, f"{h:.0f}" + " [" + "$kJ/kg_{da}$" + "]", backgroundcolor="white",
                        fontsize=6)
            x_text += 3
            y_text += 2.7

        self._psychro_chart = (fig, ax)

    def print_psychro_chart(self):
        """ Just a function to print the psychrometric chart with current transformations
        """
        if not hasattr(self, '_psychro_chart'):
            self._psychro_plot()

        fig, ax = self._psychro_chart

        self.values = np.array([
            [self._chart_T_ext, self._chart_x_ext * 1000],
            [self._chart_T_hr, self._chart_x_hr * 1000],
            [self._chart_T_mix, self._chart_x_mix * 1000],
            [self._chart_T_preh_deu, self._chart_x_preh_deu * 1000],
            [self._chart_T_as, self._chart_x_as * 1000],
            [self._chart_T_posth, self._chart_x_posth * 1000]])

        self.values_tz = np.array([
            [self._chart_T_zone, self._chart_x_zone * 1000],
            [self._chart_T_mix, self._chart_x_mix * 1000]
        ])

        ax.plot(self.values[:, 0], self.values[:, 1], 'r-o', fillstyle="none")
        ax.plot(self.values_tz[:, 0], self.values_tz[:, 1], 'k--o', linewidth=0.6, fillstyle="none")
        try:
            import matplotlib.pyplot as plt
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                "To run the AirHandlingUnit._psychro_plot you need to install matplotlib package")
        plt.show()

#%%---------------------------------------------------------------------------------------------------
# AHU class
class AirHandlingUnit(_BaseAirHandlingUnit):
    '''This class manages the air handling unit.
    Some general variables are set as class variables while the __init__ memorizes the inputs
    '''
    
    def __init__(self,
                 name: str,
                 mechanical_vent: MechanicalVentilation,
                 supply_temperature: Schedule,
                 supply_specific_humidity: Schedule,
                 ahu_operation: Schedule,
                 humidity_control: bool,
                 sensible_heat_recovery_eff: float,
                 latent_heat_recovery_eff: float,
                 outdoor_air_ratio: float,
                 weather: WeatherFile,
                 thermal_zone,

                 tag: str = None,
                 ):
        """Air Handling Unit Constructor: creates the AHU object and memorizes the attributes (using properties set methods tho check types)


        Parameters
        ----------
        name : str
            name of the Air Handling Unit
        mechanical_vent : eureca_building.ventilation.MechanicalVentilation
            ventialation object to define air flow rate
        supply_temperature : eureca_building.schedule.Schedule
            Schedule object
        supply_specific_humidity : eureca_building.schedule.Schedule
            Schedule object
        ahu_operation : eureca_building.schedule.Schedule
            Schedule object to define opeartion (-1 cooling, 1 heating, 0 fan mode)
        humidity_control : bool
            whether do humidification/dehumidification
        sensible_heat_recovery_eff : float
            sensible heat recovery efficiency, must be between 0 and 1
        latent_heat_recovery_eff : float
            sensible heat recovery efficiency, must be between 0 and 1
        outdoor_air_ratio : float
            outdoor air fraction, must be between 0 and 1
        weather : eureca_building.weather.WeatherFile
            Weather object
        thermal_zone : eureca_building.thermal_zone.ThermalZone
            ThermalZone object
        tag : str
            possible tags

        Raises
        ------
        TypeError
            checks the input type
        ValueError
            checks the input type
        """
        
        # Check input data type

        if not isinstance(name, str):
            raise TypeError(f'ERROR AHU inizialization, name must be an string: Name {name}')
        
        # Inizialization
        self.ahu_name = name
        self.mechanical_ventilation = mechanical_vent
        self.supply_temperature = supply_temperature
        self.supply_specific_humidity = supply_specific_humidity
        self.ahu_operation = (ahu_operation, weather)
        self.humidity_control = humidity_control
        self.sensible_heat_recovery_eff = sensible_heat_recovery_eff
        self.latent_heat_recovery_eff = latent_heat_recovery_eff
        self.outdoor_air_ratio = outdoor_air_ratio
        self.tag = tag


        self.air_flow_rate_kg_S, self.vapour_flow_rate_kg_S = self.mechanical_ventilation.get_flow_rate(weather, volume = thermal_zone._volume, area = thermal_zone._net_floor_area)
        self.electric_consumption_W = self.electric_consumption_based_on_mass_flow_rate(self.air_flow_rate_kg_S)

        # Association of AHU to thermal zone
        try:
            thermal_zone.add_air_handling_unit(self, weather)
        except AttributeError:
            raise TypeError(f'ERROR AHU inizialization, thermal zone must be an ThermalZone object: thermal_zone {type(thermal_zone)}')


    @property
    def supply_temperature(self):
        return self._supply_temperature

    @supply_temperature.setter
    def supply_temperature(self, value):
        if not isinstance(value, Schedule):
            raise ValueError(f"Air Handling Unit object {self.ahu_name}, supply_temperature not Schedule: {type(value)}")
        if np.any(np.less(value.schedule, ventilation_prop["mechanical"]['temperature limit'][0])):
            logging.warning(
                f"Air Handling Unit object {self.ahu_name}, supply temperature schedule goes below {ventilation_prop['mechanical']['temperature limit'][0]} °C")
        if np.any(np.greater(value.schedule, ventilation_prop["mechanical"]['temperature limit'][1])):
            logging.warning(
                f"Air Handling Unit object {self.ahu_name}, supply temperature schedule goes above {ventilation_prop['mechanical']['temperature limit'][1]} °C")
        self._supply_temperature = copy.deepcopy(value)

    @property
    def supply_specific_humidity(self):
        return self._supply_specific_humidity

    @supply_specific_humidity.setter
    def supply_specific_humidity(self, value):
        if not isinstance(value, Schedule):
            raise ValueError(f"Air Handling Unit object {self.ahu_name}, supply_specific_humidity not Schedule: {type(value)}")
        if np.any(np.less(value.schedule, ventilation_prop["mechanical"]['specific humidity limit'][0])):
            logging.warning(
                f"Air Handling Unit object {self.ahu_name}, supply specific humidity schedule goes below {ventilation_prop['mechanical']['specific humidity limit'][0]} kg_v/kg_as")
        if np.any(np.greater(value.schedule, ventilation_prop["mechanical"]['specific humidity limit'][1])):
            logging.warning(
                f"Air Handling Unit object {self.ahu_name}, supply specific humidity schedule goes above {ventilation_prop['mechanical']['specific humidity limit'][1]} kg_v/kg_as")
        self._supply_specific_humidity = copy.deepcopy(value)

    @property
    def ahu_operation(self):
        return self._ahu_operation

    @ahu_operation.setter
    def ahu_operation(self, value):
        sched = value[0]
        weather = value[1]
        if not isinstance(sched, Schedule):
            raise ValueError(f"Air Handling Unit object {self.ahu_name}, ahu_operation not Schedule: {type(sched)}")
        if np.any(np.not_equal(sched.schedule,1)*np.not_equal(sched.schedule,0)*np.not_equal(sched.schedule,-1)):
            raise ScheduleOutsideBoundaryCondition(
                f"Air Handling Unit object {self.ahu_name}, ahu_operation schedule with values not allowed: allowed values = [-1,0,1]")
        self._ahu_operation = sched

        # Needed to solve correctly thermal zone in free cooling
        ext_air_db_temp = weather.hourly_data['out_air_db_temperature']
        ext_air_spec_hum = weather.hourly_data['out_air_specific_humidity']
        filt = self._ahu_operation.schedule == 0
        self.supply_temperature.schedule[filt] = ext_air_db_temp[filt]
        self.supply_specific_humidity.schedule[filt] = ext_air_spec_hum[filt]


    @property
    def humidity_control(self):
        return self._humidity_control

    @humidity_control.setter
    def humidity_control(self, value):
        if not isinstance(value, bool):
            raise ValueError(f"Air Handling Unit object {self.ahu_name}, humidity control not bool: {type(value)}")
        self._humidity_control = value

    @property
    def sensible_heat_recovery_eff(self):
        return self._sensible_heat_recovery_eff

    @sensible_heat_recovery_eff.setter
    def sensible_heat_recovery_eff(self, value):
        if not isinstance(value, float):
            raise ValueError(f"Air Handling Unit object {self.ahu_name}, sensible_heat_recovery_eff not float: {type(value)}")
        if value > 1. or value < 0.:
            raise PropertyOutsideBoundaries(
                f"Air Handling Unit object {self.ahu_name}, sensible_heat_recovery_eff in range [0.,1]: {type(value)}")

        self._sensible_heat_recovery_eff = value

    @property
    def latent_heat_recovery_eff(self):
        return self._latent_heat_recovery_eff

    @latent_heat_recovery_eff.setter
    def latent_heat_recovery_eff(self, value):
        if not isinstance(value, float):
            raise ValueError(
                f"Air Handling Unit object {self.ahu_name}, latent_heat_recovery_eff not float: {type(value)}")
        if value > 1. or value < 0.:
            raise PropertyOutsideBoundaries(
                f"Air Handling Unit object {self.ahu_name}, latent_heat_recovery_eff in range [0.,1]: {type(value)}")

        self._latent_heat_recovery_eff = value

    @property
    def outdoor_air_ratio(self):
        return self._outdoor_air_ratio

    @outdoor_air_ratio.setter
    def outdoor_air_ratio(self, value):
        if not isinstance(value, float):
            raise ValueError(
                f"Air Handling Unit object {self.ahu_name}, outdoor_air_ratio not float: {type(value)}")
        if value > 1. or value < 0.:
            raise PropertyOutsideBoundaries(
                f"Air Handling Unit object {self.ahu_name}, outdoor_air_ratio in range [0.,1]: {type(value)}")

        self._outdoor_air_ratio = value
        
    
    def air_handling_unit_calc(self,
                               t,
                               weather,
                               T_int,
                               x_int,
                               ):
        """Solution of the time step calculation. It uses outdoor conditions (from WeatherFile), and zone conditions (from zone)

        Parameters
        ----------
        t : int
            timestep: int [-]
        weather : eureca_building.weather.WeatherFile
            WeatherFile object
        T_int : float
            zone internal temperature: float [°C]
        x_int : float
            zone internal specific humidity: float [kg_v/kg_da]


        """

        # Check input data type 
        
        if not isinstance(t, int):
            raise TypeError(f'ERROR AHUCalc, bd {self.ahu_name}, time step {t}, input t is not an interger: t {t}')

        T_ext = weather.hourly_data['out_air_db_temperature'][t]
        x_ext = weather.hourly_data['out_air_specific_humidity'][t]

        self._chart_T_ext, self._chart_x_ext = T_ext, x_ext
        self._chart_T_zone, self._chart_x_zone = T_int, x_int

        AHU_operation = self.ahu_operation.schedule[t]
        self.T_sup = self.supply_temperature.schedule[t]
        self.x_sup = self.supply_specific_humidity.schedule[t]
        m_vent = self.air_flow_rate_kg_S[t]

        OutAirRatio = self.outdoor_air_ratio

        # Saturation conditions check
        sat_cond, psat = self.checkSatCond(T_ext,x_ext,self.p_atm)
        if sat_cond == False:
            x_ext = 0.99*0.622*psat/(self.p_atm-0.99*psat)
            if AHU_operation == 0:
                self.x_sup = x_ext
        sat_cond, psat = self.checkSatCond(T_int,x_int,self.p_atm)
        if sat_cond == False:
            logging.warning(f'ERROR  AHUCalc method, bd {self.ahu_name}, time step {t}, Zone conditions outside saturation limit')
        sat_cond, psat = self.checkSatCond(self.T_sup,self.x_sup,self.p_atm)
        if sat_cond == False:
            logging.warning(f'ERROR:  AHUCalc method, bd {self.ahu_name}, time step {t}, Supply conditions outside saturation limit')
        
        
        # Pre-processing on Heat Recovery and Mixer
        # Enthalpy in J/(kg)
        self.h_ext = self.cp_air*T_ext + (self.r_0 + self.cpv*T_ext)*x_ext
        self.h_z = self.cp_air*T_int + (self.r_0 + self.cpv*T_int)*x_int
        
        
        # Heat Recovery Bypass in condition of free heating or cooling
        if (T_int-T_ext)*AHU_operation > 0:
            self.T_hr = T_ext + self.sensible_heat_recovery_eff*(T_int-T_ext)
            self.T_out = T_int - self.sensible_heat_recovery_eff*(T_int-T_ext)
        else:
            self.T_hr = T_ext
            self.T_out = T_int
            # This is the case of free floating (sets the outdoor air ratio to 1)
            OutAirRatio = 1
        if (x_int-x_ext)*AHU_operation > 0:
            self.x_hr = x_ext + self.latent_heat_recovery_eff*(x_int-x_ext)
        else:
            self.x_hr = x_ext
        # Enthalpy in J/(kg)
        self.h_hr = self.cp_air*self.T_hr +(self.r_0 + self.cpv*self.T_hr)*self.x_hr
        
        # Mixer
        # Enthalpy in J/(kg)
        self.h_mix = OutAirRatio*self.h_hr + (1 - OutAirRatio)*self.h_z
        self.x_mix = OutAirRatio*self.x_hr + (1 - OutAirRatio)*x_int
        self.T_mix = (self.h_mix - self.r_0*self.x_mix)/(self.cp_air + self.cpv*self.x_mix)
        
        
        # BATTERIES DEMAND CALCULATION
        
        # SENSIBLE AND LATENT CONTROL MODE
        if self.humidity_control == True:
        
            # Heating mode
            if AHU_operation == 1:
                                
                
                # Adiabatic Saturator
                self.x_as = self.x_sup
                p_as = self.p_atm*self.x_as/(0.622+self.x_as)
                if p_as >= 610.5:
                    self.T_as = 237.3*np.log(p_as/610.5)/(17.269-np.log(p_as/610.5))
                else:
                    self.T_as = 265.5*np.log(p_as/610.5)/(21.875-np.log(p_as/610.5))
                # Enthalpy in J/(kg)
                self.h_as = self.cp_air*self.T_as + (self.r_0 + self.cpv*self.T_as)*self.x_as
                
                
                # Pre-Heater control
                if self.h_mix < self.h_as:
                    self.h_ph = self.h_as
                    self.x_ph = self.x_mix
                    self.T_ph = (self.h_ph-self.r_0*self.x_ph)/(self.cp_air+self.cpv*self.x_ph)
                else:
                    # In this case the Pre-Heater is by-passed
                    self.h_ph = self.h_mix
                    self.x_ph = self.x_mix
                    self.T_ph = self.T_mix
                # Enthalpy in J/(kg)
                self.h_as = self.h_ph
                
                
                # Humidification need check
                if self.x_ph > self.x_sup:
                    self.x_as = self.x_ph
                    self.x_sup = self.x_as
                    self.T_as = (self.h_as - self.r_0*self.x_as)/(self.cp_air+self.cpv*self.x_as)
                
                
                # Controlling T_sup > T_as
                if self.T_sup < self.T_as:
                    self.T_sup = self.T_as
                # Enthalpy in J/(kg)
                self.h_sup = self.cp_air*self.T_sup+(self.r_0+self.cpv*self.T_sup)*self.x_sup
                
                
                # Batteries Demand [W]
                # Pre-Heater
                self.preh_deu_Dem = m_vent*(self.h_ph-self.h_mix)
                self.preh_deu_Dem_sens = m_vent*self.cp_air*(self.T_ph-self.T_mix)
                self.preh_deu_Dem_lat = m_vent*self.r_0*(self.x_ph-self.x_mix)
                # Adiabatic Saturator
                self.sat_Dem = m_vent*(self.h_as-self.h_ph)
                self.sat_Dem_lat = m_vent*self.r_0*(self.x_as-self.x_ph)
                self.sat_Dem_sens = self.sat_Dem-self.sat_Dem_lat
                # Post-Heater
                self.posth_Dem = m_vent*(self.h_sup-self.h_as)
                self.posth_Dem_sens = m_vent*self.cp_air*(self.T_sup-self.T_as)
                self.posth_Dem_lat = m_vent*self.r_0*(self.x_sup - self.x_as)
                # Total Demand
                self.AHU_demand = self.preh_deu_Dem + self.sat_Dem + self.posth_Dem
                self.AHU_demand_sens = self.preh_deu_Dem_sens + self.sat_Dem_sens + self.posth_Dem_sens
                self.AHU_demand_lat = self.preh_deu_Dem_lat + self.sat_Dem_lat + self.posth_Dem_lat

                # Conditions for chart
                self._chart_T_hr,       self._chart_x_hr =          self.T_hr, self.x_hr
                self._chart_T_mix,      self._chart_x_mix =         self.T_mix, self.x_mix
                self._chart_T_preh_deu, self._chart_x_preh_deu =    self.T_ph, self.x_ph
                self._chart_T_as,       self._chart_x_as =          self.T_as, self.x_as
                self._chart_T_posth,    self._chart_x_posth =       self.T_sup, self.x_sup
                
                
                # Cooling mode
            elif AHU_operation == -1:
                
                # Dehumidificator
                self.x_de = self.x_sup
                p_de = self.p_atm*self.x_de/(0.622+self.x_de)
                if p_de >= 610.5:
                    self.T_de = 237.3*np.log(p_de/610.5)/(17.269-np.log(p_de/610.5))
                else:
                    self.T_de = 265.5*np.log(p_de/610.5)/(21.875-np.log(p_de/610.5))
                
                
                # Dehumidifcator and post-heater control
                if self.x_de >= self.x_mix:
                    self.x_de = self.x_mix
                    self.x_sup = self.x_de
                    if self.T_sup <= self.T_mix:
                        self.T_de = self.T_sup
                    else:
                        self.T_de = self.T_mix
                        # self.T_sup = self.T_de  # Why???????
                self.h_de = self.cp_air*self.T_de+(self.r_0+self.cpv*self.T_de)*self.x_de
                self.h_sup = self.cp_air*self.T_sup+(self.r_0+self.cpv*self.T_sup)*self.x_sup
                
                
                # Batteries Demand [kW]
                # Dehumidificator
                self.preh_deu_Dem = m_vent*(self.h_de - self.h_mix)
                self.preh_deu_Dem_sens = m_vent*self.cp_air*(self.T_de - self.T_mix)
                self.preh_deu_Dem_lat = m_vent*self.r_0*(self.x_de - self.x_mix)
                self.sat_Dem = 0.
                self.sat_Dem_lat = 0.
                self.sat_Dem_sens = 0.
                # Post-Heater
                self.posth_Dem = m_vent*(self.h_sup - self.h_de)
                self.posth_Dem_sens = m_vent*self.cp_air*(self.T_sup - self.T_de)
                self.posth_Dem_lat = m_vent*self.r_0*(self.x_sup - self.x_de)
                # Total Demand
                self.AHU_demand = self.preh_deu_Dem + self.posth_Dem
                self.AHU_demand_sens = self.preh_deu_Dem_sens + self.posth_Dem_sens
                self.AHU_demand_lat = self.preh_deu_Dem_lat + self.posth_Dem_lat

                # Conditions for chart
                self._chart_T_hr,       self._chart_x_hr =          self.T_hr, self.x_hr
                self._chart_T_mix,      self._chart_x_mix =         self.T_mix, self.x_mix
                self._chart_T_preh_deu, self._chart_x_preh_deu =    self.T_de, self.x_de
                self._chart_T_as,       self._chart_x_as =          self.T_de, self.x_de
                self._chart_T_posth,    self._chart_x_posth =       self.T_sup, self.x_sup
                

            elif AHU_operation == 0:
                
                # Plant Off
                self.T_hr, self.T_mix, self.T_sup = T_ext, T_ext, T_ext
                self.x_hr, self.x_mix, self.x_sup = x_ext, x_ext, x_ext
                self.h_hr, self.h_mix, self.h_sup = self.h_ext, self.h_ext, self.h_ext

                # Pre-Heater
                self.preh_deu_Dem = 0.
                self.preh_deu_Dem_sens = 0.
                self.preh_deu_Dem_lat = 0.
                # Adiabatic Saturator
                self.sat_Dem = 0.
                self.sat_Dem_lat = 0.
                self.sat_Dem_sens = 0.
                # Post-Heater
                self.posth_Dem = 0.
                self.posth_Dem_sens = 0.
                self.posth_Dem_lat = 0.

                # Batteries Demand [kW]
                self.AHU_demand = 0
                self.AHU_demand_sens = 0
                self.AHU_demand_lat = 0

                # Conditions for chart
                self._chart_T_hr,       self._chart_x_hr =          self.T_hr, self.x_hr
                self._chart_T_mix,      self._chart_x_mix =         self.T_mix, self.x_mix
                self._chart_T_preh_deu, self._chart_x_preh_deu =    self.T_mix, self.x_mix
                self._chart_T_as,       self._chart_x_as =          self.T_mix, self.x_mix
                self._chart_T_posth,    self._chart_x_posth =       self.T_mix, self.x_mix
            
            else:
                print('AHUOnOff value not allowed')
        
        
        # SENSIBLE CONTROL MODE
        if self.humidity_control == False:
            
            
            # Heating mode
            if AHU_operation == 1:
                
                
                # Pre-Heater and Adiabatic Saturator doesn't exist!!
                self.h_ph = self.h_mix
                self.x_ph = self.x_mix
                self.T_ph = self.T_mix
                self.h_as = self.h_ph
                self.x_as = self.x_ph
                self.T_as = self.T_ph
                
                
                # Post-Heater
                if self.T_sup <self.T_as:
                    self.T_sup = self.T_as
                
                self.x_sup = self.x_as
                self.h_sup = self.cp_air*self.T_sup+(self.r_0+self.cpv*self.T_sup)*self.x_sup

                # Pre-Heater
                self.preh_deu_Dem = 0.
                self.preh_deu_Dem_sens = 0.
                self.preh_deu_Dem_lat = 0.
                # Adiabatic Saturator
                self.sat_Dem = 0.
                self.sat_Dem_lat = 0.
                self.sat_Dem_sens = 0.
                # Post-Heater
                self.posth_Dem = m_vent * self.cp_air * (self.T_sup - self.T_as)
                self.posth_Dem_sens = m_vent * self.cp_air * (self.T_sup - self.T_as)
                self.posth_Dem_lat = 0.


                # Batteries Demand [kW]
                self.AHU_demand = m_vent*self.cp_air*(self.T_sup-self.T_as)
                self.AHU_demand_sens = self.AHU_demand
                self.AHU_demand_lat = 0

                # Conditions for chart
                self._chart_T_hr,       self._chart_x_hr =          self.T_hr, self.x_hr
                self._chart_T_mix,      self._chart_x_mix =         self.T_mix, self.x_mix
                self._chart_T_preh_deu, self._chart_x_preh_deu =    self.T_ph, self.x_ph
                self._chart_T_as,       self._chart_x_as =          self.T_as, self.x_as
                self._chart_T_posth,    self._chart_x_posth =       self.T_sup, self.x_sup

    
            elif AHU_operation == -1:
                # Cooling mode

                # Pre-Heater and Adiabatic Saturator doesn't exist!!
                self.h_ph = self.h_mix
                self.x_ph = self.x_mix
                self.T_ph = self.T_mix
                self.h_as = self.h_ph
                self.x_as = self.x_ph
                self.T_as = self.T_ph
                
                if self.T_sup > self.T_mix:
                    self.T_sup = self.T_mix
                    self.x_sup = self.x_mix
                    self.h_sup = self.h_mix
                else:
                    SatCond, psat = self.checkSatCond(self.T_sup,self.x_mix,self.p_atm)
                    if SatCond == False:
                        # print('Info: supply temperature is too low -- changed')

                        if self.T_sup < 0:
                            p_sat = 610.5*np.exp((21.875*self.T_sup)/(265.5+self.T_sup))
                        else:
                            p_sat = 610.5*np.exp((17.269*self.T_sup)/(237.3+self.T_sup))
                        self.x_sup = 0.99*0.622*p_sat/(self.p_atm-0.99*p_sat)
                        self.h_sup = self.cp_air*self.T_sup+(self.r_0+self.cpv*self.T_sup)*self.x_sup
                    else:
                        self.x_sup= self.x_mix
                        self.h_sup = self.cp_air*self.T_sup+(self.r_0+self.cpv*self.T_sup)*self.x_sup

                # Pre-Heater
                self.preh_deu_Dem = 0.
                self.preh_deu_Dem_sens = 0.
                self.preh_deu_Dem_lat = 0.
                # Adiabatic Saturator
                self.sat_Dem = 0.
                self.sat_Dem_lat = 0.
                self.sat_Dem_sens = 0.
                # Post-Heater
                self.posth_Dem = m_vent * self.cp_air * (self.T_sup - self.T_as)
                self.posth_Dem_sens = m_vent * self.cp_air * (self.T_sup - self.T_as)
                self.posth_Dem_lat = 0.
                
                #  Batteries Demand [kW]
                self.AHU_demand = m_vent*(self.h_sup - self.h_mix)
                self.AHU_demand_sens = self.AHU_demand
                self.AHU_demand_lat = 0

                # Conditions for chart
                self._chart_T_hr,       self._chart_x_hr =          self.T_hr, self.x_hr
                self._chart_T_mix,      self._chart_x_mix =         self.T_mix, self.x_mix
                self._chart_T_preh_deu, self._chart_x_preh_deu =    self.T_de, self.x_de
                self._chart_T_as,       self._chart_x_as =          self.T_de, self.x_de
                self._chart_T_posth,    self._chart_x_posth =       self.T_sup, self.x_sup
                
                
            elif AHU_operation == 0:
                # Plant OFF
                
                self.T_hr, self.T_mix, self.T_sup = T_ext, T_ext, T_ext
                self.x_hr, self.x_mix, self.x_sup = x_ext, x_ext, x_ext
                self.h_hr, self.h_mix, self.h_sup = self.h_ext, self.h_ext, self.h_ext

                # Pre-Heater
                self.preh_deu_Dem = 0.
                self.preh_deu_Dem_sens = 0.
                self.preh_deu_Dem_lat = 0.
                # Adiabatic Saturator
                self.sat_Dem = 0.
                self.sat_Dem_lat = 0.
                self.sat_Dem_sens = 0.
                # Post-Heater
                self.posth_Dem = 0.
                self.posth_Dem_sens = 0.
                self.posth_Dem_lat = 0.

                # Batteries Demand [kW]
                self.AHU_demand = 0
                self.AHU_demand_sens = 0
                self.AHU_demand_lat = 0

                # Conditions for chart
                self._chart_T_hr,       self._chart_x_hr =          self.T_hr, self.x_hr
                self._chart_T_mix,      self._chart_x_mix =         self.T_mix, self.x_mix
                self._chart_T_preh_deu, self._chart_x_preh_deu =    self.T_mix, self.x_mix
                self._chart_T_as,       self._chart_x_as =          self.T_mix, self.x_mix
                self._chart_T_posth,    self._chart_x_posth =       self.T_mix, self.x_mix
            
            else:
                sys.exit('AHUOnOff value not allowed at time step: '+str(t))

# %%---------------------------------------------------------------------------------------------------
# AHU class
class HeatRecoveryUnit(_BaseAirHandlingUnit):
    '''This class manages a heat recovery unit (only fans and heat recovery, without active heating/cooling of air.
    '''

    def __init__(self,
                 name: str,
                 mechanical_vent: MechanicalVentilation,
                 sensible_heat_recovery_eff: float,
                 latent_heat_recovery_eff: float,
                 weather: WeatherFile,
                 thermal_zone,

                 tag: str = None,
                 ):
        """Heat recovery unit Constructor: creates the HRU object and memorizes the attributes (using properties set methods tho check types)


        Parameters
        ----------
        name : str
            name of the Air Handling Unit
        mechanical_vent : eureca_building.ventilation.MechanicalVentilation
            ventialation object to define air flow rate
        sensible_heat_recovery_eff : float
            sensible heat recovery efficiency, must be between 0 and 1
        latent_heat_recovery_eff : float
            sensible heat recovery efficiency, must be between 0 and 1
        weather : eureca_building.weather.WeatherFile
            Weather object
        thermal_zone : eureca_building.thermal_zone.ThermalZone
            ThermalZone object
        tag : str
            possible tags

        Raises
        ------
        TypeError
            checks the input type
        ValueError
            checks the input type
        """

        # Check input data type

        if not isinstance(name, str):
            raise TypeError(f'ERROR AHU inizialization, name must be an string: Name {name}')

        # Inizialization
        self.ahu_name = name
        self.mechanical_ventilation = mechanical_vent
        self.sensible_heat_recovery_eff = sensible_heat_recovery_eff
        self.latent_heat_recovery_eff = latent_heat_recovery_eff
        self.tag = tag

        self.air_flow_rate_kg_S, self.vapour_flow_rate_kg_S = self.mechanical_ventilation.get_flow_rate(weather,
                                                                                                        volume=thermal_zone._volume,
                                                                                                        area=thermal_zone._net_floor_area)
        self.electric_consumption_W = self.electric_consumption_based_on_mass_flow_rate(self.air_flow_rate_kg_S)

        # Association of AHU to thermal zone
        try:
            thermal_zone.add_air_handling_unit(self, weather)
        except AttributeError:
            raise TypeError(
                f'ERROR AHU inizialization, thermal zone must be an ThermalZone object: thermal_zone {type(thermal_zone)}')

    @property
    def sensible_heat_recovery_eff(self):
        return self._sensible_heat_recovery_eff

    @sensible_heat_recovery_eff.setter
    def sensible_heat_recovery_eff(self, value):
        if not isinstance(value, float):
            raise ValueError(
                f"Air Handling Unit object {self.ahu_name}, sensible_heat_recovery_eff not float: {type(value)}")
        if value > 1. or value < 0.:
            raise PropertyOutsideBoundaries(
                f"Air Handling Unit object {self.ahu_name}, sensible_heat_recovery_eff in range [0.,1]: {type(value)}")

        self._sensible_heat_recovery_eff = value

    @property
    def latent_heat_recovery_eff(self):
        return self._latent_heat_recovery_eff

    @latent_heat_recovery_eff.setter
    def latent_heat_recovery_eff(self, value):
        if not isinstance(value, float):
            raise ValueError(
                f"Air Handling Unit object {self.ahu_name}, latent_heat_recovery_eff not float: {type(value)}")
        if value > 1. or value < 0.:
            raise PropertyOutsideBoundaries(
                f"Air Handling Unit object {self.ahu_name}, latent_heat_recovery_eff in range [0.,1]: {type(value)}")

        self._latent_heat_recovery_eff = value

    def air_handling_unit_calc(self,
                               t,
                               weather,
                               T_int,
                               x_int,
                               ):
        """Solution of the time step calculation. It uses outdoor conditions (from WeatherFile), and zone conditions (from zone)

        Parameters
        ----------
        t : int
            timestep: int [-]
        weather : eureca_building.weather.WeatherFile
            WeatherFile object
        T_int : float
            zone internal temperature: float [°C]
        x_int : float
            zone internal specific humidity: float [kg_v/kg_da]


        """

        # Check input data type

        if not isinstance(t, int):
            raise TypeError(f'ERROR AHUCalc, bd {self.ahu_name}, time step {t}, input t is not an interger: t {t}')

        T_ext = weather.hourly_data['out_air_db_temperature'][t]
        x_ext = weather.hourly_data['out_air_specific_humidity'][t]

        self._chart_T_ext, self._chart_x_ext = T_ext, x_ext
        self._chart_T_zone, self._chart_x_zone = T_int, x_int

        m_vent = self.air_flow_rate_kg_S[t]

        # Saturation conditions check
        sat_cond, psat = self.checkSatCond(T_ext, x_ext, self.p_atm)
        if sat_cond == False:
            x_ext = 0.99 * 0.622 * psat / (self.p_atm - 0.99 * psat)
        sat_cond, psat = self.checkSatCond(T_int, x_int, self.p_atm)
        if sat_cond == False:
            logging.warning(
                f'ERROR  HRUCalc method, bd {self.ahu_name}, time step {t}, Zone conditions outside saturation limit')

        # Pre-processing on Heat Recovery and Mixer
        # Enthalpy in J/(kg)
        self.h_ext = self.cp_air * T_ext + (self.r_0 + self.cpv * T_ext) * x_ext
        self.h_z = self.cp_air * T_int + (self.r_0 + self.cpv * T_int) * x_int

        # Heat Recovery Bypass in condition of free heating or cooling
        self.T_hr = T_ext + self.sensible_heat_recovery_eff * (T_int - T_ext)
        self.T_out = T_int - self.sensible_heat_recovery_eff * (T_int - T_ext)
        self.x_hr = x_ext + self.latent_heat_recovery_eff * (x_int - x_ext)
        # Enthalpy in J/(kg)
        self.h_hr = self.cp_air * self.T_hr + (self.r_0 + self.cpv * self.T_hr) * self.x_hr

        # Mixer
        # Enthalpy in J/(kg)
        self.h_mix = self.h_hr
        self.x_mix = self.x_hr
        self.T_mix = self.T_hr

        self.h_ph = self.h_hr
        self.x_ph = self.x_hr
        self.T_ph = self.T_hr

        self.h_as = self.h_hr
        self.x_as = self.x_hr
        self.T_as = self.T_hr

        self.h_sup = self.h_hr
        self.x_sup = self.x_hr
        self.T_sup = self.T_hr

        # Batteries Demand [W]
        # Pre-Heater
        self.preh_deu_Dem = m_vent * (self.h_ph - self.h_mix)
        self.preh_deu_Dem_sens = m_vent * self.cp_air * (self.T_ph - self.T_mix)
        self.preh_deu_Dem_lat = m_vent * self.r_0 * (self.x_ph - self.x_mix)
        # Adiabatic Saturator
        self.sat_Dem = m_vent * (self.h_as - self.h_ph)
        self.sat_Dem_lat = m_vent * self.r_0 * (self.x_as - self.x_ph)
        self.sat_Dem_sens = self.sat_Dem - self.sat_Dem_lat
        # Post-Heater
        self.posth_Dem = m_vent * (self.h_sup - self.h_as)
        self.posth_Dem_sens = m_vent * self.cp_air * (self.T_sup - self.T_as)
        self.posth_Dem_lat = m_vent * self.r_0 * (self.x_sup - self.x_as)
        # Total Demand
        self.AHU_demand = self.preh_deu_Dem + self.sat_Dem + self.posth_Dem
        self.AHU_demand_sens = self.preh_deu_Dem_sens + self.sat_Dem_sens + self.posth_Dem_sens
        self.AHU_demand_lat = self.preh_deu_Dem_lat + self.sat_Dem_lat + self.posth_Dem_lat

        # Conditions for chart
        self._chart_T_hr, self._chart_x_hr = self.T_hr, self.x_hr
        self._chart_T_mix, self._chart_x_mix = self.T_mix, self.x_mix
        self._chart_T_preh_deu, self._chart_x_preh_deu = self.T_ph, self.x_ph
        self._chart_T_as, self._chart_x_as = self.T_as, self.x_as
        self._chart_T_posth, self._chart_x_posth = self.T_sup, self.x_sup


