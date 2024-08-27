"""
Solar Thermal Calculation
---
Created by: Mohamad
Betalab - DII, University of Padua
---

"""
from eureca_building.weather import WeatherFile



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
                 coverage_factor=0.8,
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
        Efficiency= self.Collector_parameters["efficiency_slope"]*(Fluid_inlet_temperature-Air_temperature)\
            + self.Collector_parameters["efficiency_intercept"]      
        self.gained_heat=0
        for Surface in self._surfaces: 
            tilt=Surface._height_round
            orient=Surface._azimuth_round
            poa_global = weatherobject.hourly_data_irradiances[orient][tilt]['global']
            poa_direct = weatherobject.hourly_data_irradiances[orient][tilt]['direct']
            poa_diffuse = poa_global - poa_global
            shading_vector= Surface.shading_coefficient if hasattr(Surface, "shading_coefficient") else 1.
            poa_direct = poa_direct*shading_vector
            poa_global = poa_direct+poa_diffuse
            
            Area=Surface._area
            Absorbed_heat=poa_global*Area*Efficiency # W 

            self.gained_heat+= Absorbed_heat
            
 
        
        



