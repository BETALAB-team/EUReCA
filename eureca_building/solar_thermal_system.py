"""
Solar Thermal Calculation
---
Created by: Mohamad
Betalab - DII, University of Padua
---
the collector is sized based on the BOSCH manual for albany, NY
"""
from eureca_building.weather import WeatherFile
import numpy as np


'''The Class defining a flat plate collector coupled with thermal storage tank'''
class SolarThermal_Collector():

    def __init__(self,
                 name: str,
                 surface_list: list,
                 dhw:float,
                 weatherobject: WeatherFile,
                 Fluid_inlet_temperature=12,
                 Fluid_design_max_outlet_temperature=90,
                 max_coverage_factor=0.05,
                 mount_surfaces=["Roof"],
                 Collector_parameters={'efficiency_slope':-0.013,
                                       'efficiency_intercept':0.8}
                 ): 
        """
        Models a flat-plate solar thermal collector system based on BOSCH manual sizing.

        Attributes
        ----------
        name : str
            Name of the system.
        surface_list : list
            List of building surfaces to mount the collectors.
        dhw : float
            Daily domestic hot water demand [kWh/day].
        weatherobject : WeatherFile
            Weather object with irradiance and temperature data.
        ...
        """
        
        self.name=name
        self.max_coverage_factor=max_coverage_factor
        self._surfaces=[s for s in surface_list if s.surface_type in mount_surfaces]
        self.weather=weatherobject._epw_hourly_data
        self.weather_md=weatherobject._epw_general_data
        self.Collector_parameters=Collector_parameters
        city_irrad= weatherobject.general_data['yearly_solar_irradiation']/1000
        Air_temperature=self.weather['temp_air']
        Albany_city_convert=1569/city_irrad
        dhw=dhw
        slope_sizing=0.0104*Albany_city_convert
        sized_area=slope_sizing*dhw



        self.gained_heat=0
        tot_area=0
        for Surface in self._surfaces:
            tot_area=tot_area+Surface._area
        self.sized_coverage_factor=sized_area/tot_area
        coverage_factor=min(self.sized_coverage_factor,max_coverage_factor)
        for Surface in self._surfaces: 
            tilt=Surface._height_round
            orient=Surface._azimuth_round
            poa_global = weatherobject.hourly_data_irradiances[orient][tilt]['global']
            poa_direct = weatherobject.hourly_data_irradiances[orient][tilt]['direct']
            global_design=np.quantile(poa_global[poa_global>0],0.9)
            a=self.Collector_parameters["efficiency_intercept"] 
            b=self.Collector_parameters["efficiency_slope"]
            Efficiency_design=a+b*(0.5*Fluid_inlet_temperature+0.5*Fluid_design_max_outlet_temperature-Air_temperature)
            nominator=a+b*(Fluid_inlet_temperature-Air_temperature)
            
            poa_diffuse = poa_global - poa_direct
            shading_vector= Surface.shading_coefficient if hasattr(Surface, "shading_coefficient") else 1.
            poa_direct = poa_direct*shading_vector
            poa_global = poa_direct+poa_diffuse
            denominator=1-0.5*b*poa_global/(global_design*Efficiency_design)*(Fluid_design_max_outlet_temperature-Fluid_inlet_temperature)
            Efficiency=nominator/denominator
            # fluid_outlet_temperature=Fluid_inlet_temperature\
            #     +(Fluid_design_max_outlet_temperature)/global_design*poa_global
            # fluid_mean_temperature=0.5*(fluid_outlet_temperature+Fluid_inlet_temperature)
            # poa_global_standardization=poa_global

            # Efficiency= self.Collector_parameters["efficiency_slope"]*(fluid_mean_temperature-Air_temperature)\
            #     + self.Collector_parameters["efficiency_intercept"] 


            
            max_Area=Surface._area*coverage_factor
            Absorbed_heat=max_Area*poa_global*Efficiency # W 
            
            self.gained_heat+=Absorbed_heat
            
 
        
        



