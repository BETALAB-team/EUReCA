"""
Solar Thermal Calculation
---
Created by: Mohamad
Betalab - DII, University of Padua
---

"""
from eureca_building.weather import WeatherFile
import numpy as np


'''The Class defining a flat plate collector coupled with thermal storage tank'''
class SolarThermal_Collector():
    '''
    Reference: The model available from [book]
    '''
    def __init__(self,
                 name: str,
                 surface_list: list,
                 weatherobject: WeatherFile,
                 Fluid_inlet_temperature=12,
                 Fluid_design_max_outlet_temperature=90,
                 coverage_factor=0.05,
                 mount_surfaces=["Roof"],
                 Collector_parameters={'efficiency_slope':-0.013,
                                       'efficiency_intercept':0.8}
                 ): 
        
        
        self.name=name
        self.coverage_factor=coverage_factor
        self._surfaces=[s for s in surface_list if s.surface_type in mount_surfaces]
        self.weather=weatherobject._epw_hourly_data
        self.weather_md=weatherobject._epw_general_data
        self.Collector_parameters=Collector_parameters
 
        Air_temperature=self.weather['temp_air']
     
        self.gained_heat=0

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


            
            Area=Surface._area*coverage_factor
            Absorbed_heat=Area*poa_global*Efficiency # W 
            
            self.gained_heat+=Absorbed_heat
            
 
        
        



