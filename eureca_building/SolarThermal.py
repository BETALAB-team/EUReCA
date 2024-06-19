"""
Solar Thermal Calculation
---
Created by: Mohamad
Betalab - DII, University of Padua
---

"""
import pandas as pd
import numpy as np
import math
import pvlib
from eureca_building.weather import WeatherFile
from eureca_building.config import CONFIG


'''The Class defining a flat plate collector coupled with thermal storage tank'''
class Flate_plate_collector_with_storage():
    def __init__(self,
                 name: str,
                 weatherobject: WeatherFile,
                 surface_list: list,
                 mount_surfaces=["Roof"],
                 coverage_factor=0.8,
                 Collector_parameters={'Plate_Extinction_coefficient':0.99924,
                                 'Plate_thickness':-5.49097,
                                 'Plate_refraction_coefficient':0.01918,
                                 'Normal_Absorption_coefficient':0.06999,
                                 'Normal_Collector_efficiency_factor':0.75,
                                 'Collector_Heat_loss_coefficient':6.5,
                                 'Collector_Working_fluid_Heating_Capacity':4.218,
                                 'Reflection_coefficient_ground':0.26144},
                 Storage_parameters={}): #Coverage factor for the area that is covered by the PV
        
        
        self.name=name
        self.coverage_factor=coverage_factor
        self._surfaces=[s for s in surface_list if s.surface_type in mount_surfaces]
        self.weather=weatherobject._epw_hourly_data
        self.weather_md=weatherobject._epw_general_data
        self.FPC_parameters=Collector_parameters
 
        
 
        irradiances=self.weather.hourly_data_irradiance
        self.Surface_spec_absorptions={}
        for Surface in self._surfaces:    
            SunData=pd.DataFrame()
            Absorption_per_area_for_surface=pd.DataFrame()
            Surface_Irradiance_data=irradiances[float(Surface._azimuth_round)][float(Surface._height_round)]
            SunData['I_Beam']=Surface_Irradiance_data['direct']
            SunData['I_diffuse']=Surface_Irradiance_data['global']-Surface_Irradiance_data['direct']
            SunData['Incidence_Angle']= Surface_Irradiance_data['AOI']
            Absorption_per_area_for_surface = SunData.apply(lambda row: self.Absorbed_Heat_Flat_Plate(Incidence_Angle=row['AOI'],
                                                                                                      Slope_Angle=Surface._height_round,
                                                                                                      I_beam=row['I_Beam'],
                                                                                                      I_diffuse=row['I_diffuse']), axis=1)
            self.Surface_spec_absorptions[Surface]=Absorption_per_area_for_surface
        
        for Surface in self._surfaces:    
            WeatherData=pd.DataFrame()
            ##Absorption_per_area_for_surface=pd.DataFrame()
            Surface_Irradiance_data=irradiances[float(Surface._azimuth_round)][float(Surface._height_round)]
            WeatherData['Air_Temperature']=self.weather['temp_air']
            Exposure_loss_coefficient=(1-math.exp(-collector_heat_loss_coefficient*Area/\
                                                      Fluid_mass_flow_rate/Fluid_thermal_capacity*\
                                                          collector_efficiency_factor))
            self.Surface_spec_absorptions[Surface]=Absorption_per_area_for_surface
        
        
            
        
        def Absorbed_Heat_Flat_Plate(self,
                                     Incidence_Angle,Slope_Angle,
                                     I_beam,I_diffuse,
                                     Extinction=self.FPC_parameters['Plate_Extinction_coefficient'],
                                     Thickness=self.FPC_parameters['Plate_thickness'],
                                     Refraction_coefficient=self.FPC_parameters['Plate_refraction_coefficient'],
                                     Normal_Absorption_coeffictient=self.FPC_parameters['Normal_Absorption_coefficient'],
                                     Reflection_coefficient_ground=self.FPC_parameters['Reflection_coefficient_ground']):
            
            
            '''A function that calculates effective absorption coefficient for a 
            given effective angle'''
            def tau_alpha(theta,
                          Normal_Absorption=Normal_Absorption_coeffictient,
                          Refraction=Refraction_coefficient,
                          Extinct=Extinction,
                          Thick=Thickness):
                
                Absorption_fac=1 + (2.0345) * 10 ** (-3) * theta ** (1) +\
                                (-1.993) * 10 ** (-4) * theta ** (2) +\
                                (5.324) * 10 ** (-6) * theta ** (3) +\
                                (-4.799) * 10 ** (-8) * theta ** (4) 
                Alpha=Absorption_fac*Normal_Absorption
                theta_rad=max(math.radians(theta),0.01)
                theta_refract=math.asin(math.sin(theta_rad)/Refraction)
                refract_parallel=(math.tan(theta_refract-theta_rad))**2/(math.tan(theta_refract+theta_rad))**2
                refract_prependicular=(math.sin(theta_refract-theta_rad))**2/(math.sin(theta_refract+theta_rad))**2
                tau_r= 0.5*((1-refract_parallel)/(1+refract_parallel)+(1-refract_prependicular)/(1+refract_prependicular))
                tau_alpha=math.exp(-Extinct*Thick/math.cos(theta_refract))
                tau_tot=tau_r*tau_alpha
                tau_alpha_non_reflected=tau_tot*Alpha
                tau_alpha_reflected=tau_alpha_non_reflected*1.01
                return tau_alpha_reflected
            
            ## Calculate Effective Angles for Diffuse And Ground
            Diffuse_effective_angle=59.68-0.1388*Slope_Angle+0.001497*Slope_Angle**2  
            Ground_effective_angle=90-0.5789*Slope_Angle+0.002693*Slope_Angle**2
            Absorbed_Beam=I_beam*tau_alpha(Incidence_Angle)
            Absorbed_Diffuse=I_diffuse*tau_alpha(Diffuse_effective_angle)*\
                                (1+math.cos(Slope_Angle))/2
            Absorbed_Ground=Reflection_coefficient_ground*(I_beam+I_diffuse)*\
                                tau_alpha(Ground_effective_angle)*\
                                (1-math.cos(Slope_Angle))/2
            Absorbed_total_power=Absorbed_Beam+Absorbed_Diffuse+Absorbed_Ground
            return Absorbed_total_power 
        
        

Fluid_mass_flow_rate=2
Fluid_inlet_temperature=12
Air_depth=0.1
Area=100

Slope=30
I_beam=100
I_diffuse=20
Incidence=0
Reflect =0.25
Air_Temperature=16
Fluid_thermal_capacity=4218
collector_efficiency_factor=0.9-0.6*(math.asin(Area/(Area+1))+1.76)/3.14
collector_heat_loss_coefficient=6.5

Exposure_loss_coefficient=(1-math.exp(-collector_heat_loss_coefficient*Area/\
                                          Fluid_mass_flow_rate/Fluid_thermal_capacity*\
                                              collector_efficiency_factor))
Temp_rise_sun=Absorbed_heat/collector_heat_loss_coefficient
Temp_rise_air=Air_Temperature-Fluid_inlet_temperature
Temp_rise=max((-Temp_rise_air+Temp_rise_sun)*Exposure_loss_coefficient,0)
Utilized_Heat_Gain=(Fluid_mass_flow_rate*Fluid_thermal_capacity*Temp_rise)
UU=Utilized_Heat_Gain/(Area*I_beam)



