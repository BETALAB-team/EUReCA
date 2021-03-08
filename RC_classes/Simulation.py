''' IMPORTING MODULES '''

import pandas as pd
import os
import numpy as np
import time as tm
from RC_classes.Envelope import loadEnvelopes
from RC_classes.EndUse import loadArchetype, loadSimpleArchetype
from RC_classes.CityJSON import JsonCity
from RC_classes.WeatherData import  Weather
from RC_classes.BuildingsPlants import loadPlants


        
#%% ---------------------------------------------------------------------------------------------------
#%% Sim class     


class Sim():
    '''
    This class manages the simulation of the district
    
    Methods:
        init
        set_input_from_dictionary
        preprocessing
        city_creation
        urban_shading
        buildings_params_and_loads
        plants_design_and_creation
        simulation
        output
        
    '''
    
    
    
    def __init__(self, name = 'city'):
        '''
        init method.
        It sets only the name 

        Parameters
        ----------
        name : str, optional
            the city name. The default is 'city'.

        Returns
        -------
        None.

        '''

        self.name = 'city'
        
    def set_input_from_text_file (self, input_path):
        '''
        This method read an input file with three dictionaries 
        Input_files : dictionary
                    Dictionary that includes these keys and values:
                            default_input_path : str
                            epw : str
                            envelopes : str
                            end-uses : str
                            plants : str
                            json : str
                            json_mod : str
                            end-uses_mod : str
        Sim_input : dictionary
                    Dictionary that includes these keys and values:
                               year : int 
                               first_day : int 
                               tz : str
                               ts : int 
                               hours : int 
                               model : str
                               DD_boundaries : np.array
                               Time_to_regime : int 
                               DesignDays_calc : bool
                               Plant_calc : bool
                               OutputReport : bool
                               
                               azSubdiv : int 
                               hSubdiv : int 
                               SolarCalc : bool
                               Shading_calc : bool
                               toll_az : float
                               toll_dist : float
                               toll_theta : float
                               R_f : float
            
        UWG_data : dictionary
                      Dictionary that includes Urban data for Urban Weather Gen
                      For instance:
                              # Urban Canyon Data
                                UWG_data = {
                                    'UWG_calc' : bool(False), 
                                    'Area' : float(538000),
                                    'Diameter': float(800),
                                    'Perimeter': float(3280),
                                    'Ortogonal_length': float(800),
                                    'h_layer_inf': 2,
                                    'h_layer_sup': 150,
                                    'z_i_day': 1000,
                                    'z_i_night': 50,
                                    'z_i_m': 10,
                                    'z_i_ref': 150,
                                    'ExtWall_radiative_coef': float(5),
                                    'ExtWall_convective_coef': float(17),
                                    'ExtWall_emissivity': float(0.85),
                                    'ExtWall_reflection': float(0.15),
                                    'Average_U_windows': float(2),
                                    'Vegetation_density': float(0.16),
                                    'Road_radiative_coef': float(5),
                                    
                                    'Road_layers': np.array([
                                        [0.05, 1.9*10**6, 273 + 10],
                                        [0.2 , 2*10**6,   273 + 12],
                                        [1   , 1.4*10**6, 273 + 13],
                                        ], dtype = float),                          # first column thickness [m], second capacity [J/m3K], first row asphalt, second stones, third gravel 
                                        
                                    'Soil_layers': np.array([
                                        [3 , 1.4*10**6, 273 + 10],
                                        [3 , 1.4*10**6, 273 + 12],
                                        [3 , 1.4*10**6, 273 + 13],
                                        ], dtype = float),                          # first column thickness [m], second capacity [J/m3K]
                                    'T_deep': float(13+273),
                                    
                                    
                                    'cars_number_over_1000_px': float(39545000/60360),
                                    'motorcycles_number_over_1000_px': float(6896000/60360),
                                    'othervehicle_number_over_1000_px': float(5364000/60360),
                                    'emission_factor_cars_over_1000_px': float(	24.74 ),
                                    'emission_factor_motorcycles_over_1000_px': float(13.34),
                                    'emission_factor_othervehicle_over_1000_px': float(104.95),
                                    'vehicle_multiplying_factor': float(1),
                                    'distance_per_hour': float(64000),
                                    'population_density': float(2283.2),
                                    'first_day': Sim_input['first_day'],
                                    'Hw_week': np.array([0.9,0.5,0.9,2.1,3.3,4.2,5.8,7.1,6.3,5.8,5.4,5.4,5.8,6.3,6.7,7.1,7.5,5.8,4.2,3.3,2.5,1.7,1.3,0.9])/100,
                                    'Hw_end': np.array([1.7,1.2,1.,1.7,2.6,3.3,4.2,5,5.8,6.2,6.7,6.7,6.7,6.7,6.3,5.8,5,4.6,4.2,3.7,3.3,2.9,2.5,2.1])/100    
                                    }

        Parameters
        ----------
        input_path : str
            path of the input file to read 

        Returns
        -------
        None.

        '''
        
        # Loading Input_files
        
        exec(open(input_path).read(), globals())
        
        
        try:
            self.input_folder = str(Input_files['default_input_path'])    
            self.epw_name = str(Input_files['epw'])    
            self.envelopes_name = str(Input_files['envelopes'])  
            self.end_uses_name = str(Input_files['end-uses'])
            self.plants_name = str(Input_files['plants'])   
            self.json_name = str(Input_files['json'])
            
            self.json_mode = str(Input_files['json_mod'])   
            self.end_uses_mode = str(Input_files['end-uses_mod'])  
            
        except KeyError:
            raise KeyError(f"""
                           Input files dictionary is incomplete, you must insert: 
                               default_input_path
                               epw
                               envelopes
                               end-uses
                               plants
                               json
                               json_mod
                               end-uses_mod
                               
                            The actual dictionary keys are: {Input_files.keys()} 
                           """)
        except ValueError:
            raise ValueError(f"""
                           There is a wrong data type in the Input_files dictionary, you must insert: 
                               default_input_path : str
                               epw : str
                               envelopes : str
                               end-uses : str
                               plants : str
                               json : str
                               json_mod : str
                               end-uses_mod : str
                               
                           The actual data are: {Input_files.values()} 
                           """)
    
        # Loading Sim_input
        
        try:
            self.year = int(Sim_input['year']) 
            self.first_day = int(Sim_input['first_day'])
            self.tz = str(Sim_input['tz'])
            self.ts = int(Sim_input['ts'])
            self.hours = int(Sim_input['hours'])
            self.model = str(Sim_input['model'])
            self.DD_boundaries = Sim_input['DD_boundaries']
            self.time_to_regime = int(Sim_input['Time_to_regime'])
            self.DoDDcalc = bool(Sim_input['DesignDays_calc'])
            self.DoPlantCalc = bool(Sim_input['Plant_calc'])
            self.StampOutputReport = bool(Sim_input['OutputReport'])
            
            
            self.azSubdiv = int(Sim_input['azSubdiv'])
            self.hSubdiv = int(Sim_input['hSubdiv'])
            self.SolarCalc = bool(Sim_input['SolarCalc'])
            
            
            self.ShadingCalc = bool(Sim_input['Shading_calc'])
            self.azToll = float(Sim_input['toll_az'])
            self.distToll = float(Sim_input['toll_dist'])
            self.thetaToll = float(Sim_input['toll_theta'])
            self.R_f = float(Sim_input['R_f'])
            
        except KeyError:
            raise KeyError(f"""
                           Sim Input dictionary is incomplete, you must insert: 
                               year
                               first_day
                               tz
                               ts
                               hours
                               model
                               DD_boundaries
                               Time_to_regime
                               DesignDays_calc
                               Plant_calc
                               OutputReport
                               
                               azSubdiv
                               hSubdiv
                               SolarCalc
                               
                               Shading_calc
                               toll_az
                               toll_dist
                               toll_theta
                               R_f

                               
                            The actual dictionary keys are: {Input_files.keys()} 
                           """)
        except ValueError:
            raise ValueError(f"""
                           There is a wrong data type in the Sim_input dictionary, you must insert: 
                               year : int 
                               first_day : int 
                               tz : str
                               ts : int 
                               hours : int 
                               model : str
                               DD_boundaries : np.array
                               Time_to_regime : int 
                               DesignDays_calc : bool
                               Plant_calc : bool
                               OutputReport : bool
                               
                               azSubdiv : int 
                               hSubdiv : int 
                               SolarCalc : bool
                               Shading_calc : bool
                               toll_az : float
                               toll_dist : float
                               toll_theta : float
                               R_f : float
                           The actual data are: {Sim_input.values()} 
                           """)

        try:
            self.UWGCalc =  bool(UWG_data['UWG_calc'])
        except:
            raise ValueError(f"UWG_data dictionary should include a UWG_calc key with a boolean attribute")
            
        # This dictionary is not checked because it is already checked in the Climate file and urban canyon class  
        self.UWG_data = UWG_data
        
        # Dictionsary with sim times
        self.times = {}

    
    def preprocessing(self):
        '''
        Pre-processing

        Returns
        -------
        None.

        '''
        
        # Cleaning warning file
        
        if os.path.isfile(os.path.join('.','OutputReport','warnings.txt')):
            os.remove(os.path.join('.','OutputReport','warnings.txt'))
        
        start = tm.time()
        
        # Creation of the weather obj
        
        self.weather = Weather(self.epw_name, input_path = self.input_folder, tz = self.tz, 
                         year = self.year, ts = self.ts, hours = self.hours, n_years = 1, 
                         irradiances_calc = self.SolarCalc, 
                         azSubdiv = self.azSubdiv, hSubdiv = self.hSubdiv, 
                         shad_tol = [self.azToll,
                                     self.distToll,
                                     self.thetaToll])
        
         
        # Loading Envelope and Schedule Data 
        
        self.envelopes = loadEnvelopes(os.path.join(self.input_folder,self.envelopes_name))                                            # Envelope file loading
        
        if self.end_uses_mode == 'Daily':
            PlantDays=[2520,3984,6192,6912]                                            # 15th April, 15th June, 15th September, 15th October
            self.sched = loadSimpleArchetype(os.path.join(self.input_folder,self.end_uses_name),
                                                     np.arange(8760),self.first_day,self.ts,PlantDays)
        elif self.end_uses_mode == 'Yearly':
            self.sched = loadArchetype(os.path.join(self.input_folder,self.end_uses_name),
                                                     np.arange(8760),self.ts)
        else:
            raise ValueError('Set a proper schedule inporting methodology (Daily or Yearly): SchedMethod {self.end_uses_mode}')
        
        # Plants List 
        
        self.Plants_list = loadPlants(os.path.join(self.input_folder,self.plants_name))
        
        end = tm.time()
        
        self.times['preprocessing'] = end - start
        
        print('Pre-processing:       ', end - start)
    
   
    def city_creation(self):
        '''
        City object creation

        Returns
        -------
        None.

        '''
        # District initialization'
        start = tm.time()
        
        self.city = JsonCity(os.path.join(self.input_folder,self.json_name),
                                      self.envelopes, self.weather, model = self.model,mode = self.json_mode)
        
        self.city.create_urban_canyon([self.weather.ts,self.weather.hours],self.UWGCalc,self.UWG_data)
        
        end = tm.time()
        
        self.times['city creation'] = end - start
        
        print('Jsoncity:             ', end - start)
        
    
    
    def urban_shading(self):
        '''
        Urban shading calculation

        Returns
        -------
        None.

        '''
        
        # Mutual shading effect evaluation
        
        start = tm.time()
        if self.ShadingCalc:
            self.city.shading_effect(self.weather.SolarPosition,
                                 mode = self.json_mode,
                                 toll_az = self.weather.Az_toll,
                                 toll_dist = self.weather.Dist_toll,
                                 toll_theta = self.weather.Theta_toll,
                                 R_f = self.R_f)
    
        end = tm.time()
                
        self.times['urban shading'] = end - start
        
        print('Shading effect TOT:   ', end - start)
        
        
    def buildings_params_and_loads(self):
        '''
        Calculation of the parameters and loads of the buildings

        Returns
        -------
        None.

        '''
        
        # Parameters and loads calculation
        
        start = tm.time()
        
        self.city.paramsandloads(self.envelopes, self.sched, self.weather, mode = self.json_mode)
        
        end = tm.time()
                
        self.times['buildings params and loads'] = end - start
        
        print('Paramscalc:           ', end - start)
        
    def plants_design_and_creation(self):
        '''
        Design power calculation of the buildings and plants size definition

        Returns
        -------
        None.

        '''
        
        # Check inputs
        
        if not self.DoDDcalc and self.DoPlantCalc:
            print('------------------------------------------------------------------------------'+
                  '\nWARNING: To do plant calculation Design Days calculation must be done!!\nSimulation will run doing the calculation of the design days\n'+
                  '------------------------------------------------------------------------------')
            self.DoDDcalc = True
            
    
        # Design days and Buildings Plants evaluation
        
        start = tm.time()
        design_days = [np.arange(self.DD_boundaries[0,0]*self.ts,self.DD_boundaries[0,1]*self.ts),
                       np.arange(self.DD_boundaries[1,0]*self.ts,self.DD_boundaries[1,1]*self.ts)] 
        if self.DoDDcalc:
            self.city.designdays(self.DoPlantCalc,self.time_to_regime,design_days,self.weather)
        
        end = tm.time()
        
        self.times['Design Days'] = end - start
        
        print('DesignDaysCalc:       ', end - start)
           
        start = tm.time()
        if self.DoPlantCalc:
            self.city.cityplants(self.Plants_list,self.weather)
        else:
            pass
        end = tm.time()
                
        self.times['Building plants'] = end - start
        
        print('Building Plant:       ', end - start)        
       
    def simulation(self):
        '''
        Simulation

        Returns
        -------
        None.

        '''      
        
        # Simulation
        sim_time = np.arange(len(self.weather.Text)) 
        start = tm.time()
        self.city.citysim(sim_time,self.weather,self.DoPlantCalc,self.Plants_list)
        end = tm.time()
                
        self.times['Simulation'] = end - start
        
        print('Simulation:           ', end - start)
    
    def output(self):
        '''
        POST-PROCESSING

        Returns
        -------
        None.

        '''
        
        start = tm.time()
        year = (8760-1)*self.ts+1
        nb = len(self.city.buildings.values())
        
        Col_Bui = ['Name', 'Age class', 'End use', 'H plant', 'C plant',
                    'BuiHeight','Footprint','ExtWallArea','nFloor','TotalArea','Volume','PDesH','PDesC',
                    'TotOpaquaA','TotWinA','TotUA',
                    'ExtWallU','RoofU','GroundU','WinU',
                    'Htr_is','Htr_w','Htr_ms','Htr_em','Cm',
                    'RrestAW','R1AW','RalphaStarAW','RalphsStarIL','RaplphaStarIW','R1IW','C1AW','C1IW']
        
        BuiInfo = pd.DataFrame(index = self.city.buildings.keys() ,columns = Col_Bui)
        
        for bd in self.city.buildings.keys():
            
            building = self.city.buildings[bd]
            z = building.zones['Zone']
            
            data = [building.name, building.archId, building.end_use, building.H_plant_type, building.C_plant_type,
                    building.buildingHeight, building.footprint, building.extWallArea, building.nFloors, building.total_area, building.Volume,
                    building.Pnom_H_BD, building.Pnom_C_BD,
              z.Tot_opaque_area, z.Tot_glazed_area, z.UA_tot,
              z.Strat['ExtWall'].U, z.Strat['Roof'].U, z.Strat['GroundFloor'].U, z.Strat['Window'].U,
              z.Htr_is, z.Htr_w, z.Htr_ms, z.Htr_em, z.Cm,
              z.RrestAW,z.R1AW,z.RalphaStarAW,z.RalphaStarIL,z.RalphaStarIW,z.R1IW,z.C1AW,z.C1IW]
            
            data_dic = dict(zip(Col_Bui,data))
            BuiInfo.loc[bd] = data_dic
        
        # Output vectors inizialization 
        HF = np.zeros([nb,year])
        LF = np.zeros([nb,year])
        AHUD = np.zeros([nb,year])
        T = np.zeros([nb,year])
        RH = np.zeros([nb,year])
        ElEn = np.zeros([nb,year])
        GasCon = np.zeros([nb,year])
        Final = np.zeros([nb,year])
        
        i=0
        columns = []
        for bd in self.city.buildings.keys():
            columns.append(bd)
            HF[i]=self.city.buildings[bd].heatFlowBD
            LF[i]=self.city.buildings[bd].latentFlowBD
            AHUD[i]=self.city.buildings[bd].AHUDemandBD
            T[i] = self.city.buildings[bd].zones['Zone'].Air_temp
            RH[i] = self.city.buildings[bd].zones['Zone'].RH_i
            ElEn[i] = self.city.buildings[bd].BDPlant.Electrical_energy_consumption
            GasCon[i] = self.city.buildings[bd].BDPlant.Gas_consumption
            Final[i] = self.city.buildings[bd].BDPlant.Final_energy_demand   
            i += 1
            
        H_F = HF.transpose()
        L_F = LF.transpose()
        AHU_D = AHUD.transpose()
        T = T.transpose()
        RH = RH.transpose()
        ElEn = ElEn.transpose()
        GasCon = GasCon.transpose()
        Final = Final.transpose()
        
        # Output 
        
        if self.StampOutputReport:
            pd.DataFrame(data = H_F, columns = columns).to_csv(os.path.join('OutputReport',self.model,'HeatFlow.csv'))
            pd.DataFrame(data = L_F, columns = columns).to_csv(os.path.join('OutputReport',self.model,'LantentFlow.csv'))
            pd.DataFrame(data = AHU_D, columns = columns).to_csv(os.path.join('OutputReport',self.model,'AHU.csv'))
            pd.DataFrame(data = T, columns = columns).to_csv(os.path.join('OutputReport',self.model,'Temp.csv'))
            pd.DataFrame(data = RH, columns = columns).to_csv(os.path.join('OutputReport',self.model,'RH.csv'))
            pd.DataFrame(data = ElEn, columns = columns).to_csv(os.path.join('OutputReport',self.model,'ElectricalEnergy.csv'))
            pd.DataFrame(data = GasCon, columns = columns).to_csv(os.path.join('OutputReport',self.model,'GasConsumption.csv'))
            # pd.DataFrame(data = Final, columns = columns).to_csv(os.path.join('OutputReport',model,'FinalEnergy.csv'))
            
            BuiInfo.to_csv(os.path.join('OutputReport',self.model,'BuildingsParams.csv'))
        
        '''
        if mode == 'geojson':
            Padua.complexmerge()
        '''
        
        end = tm.time()
                
        self.times['Postprocessing'] = end - start
        
        print('Post-processing:      ', end - start)

#%% TEST METHOD

if __name__ == '__main__':

    # Setting Input Data 
    
    Input_files = {'epw' : 'ITA_Venezia-Tessera.161050_IGDG.epw',
                   'envelopes' : 'Envelopes.xlsx',
                   'end-uses' : 'ScheduleSemp.xlsx',
                   'plants' : 'PlantsList.xlsx',
                   'json' : 'PaduaRestricted.json',
                   'default_input_path' : os.path.join('.','Input'),
                   'end-uses_mod' : 'Daily',
                   'json_mod' : 'cityjson'
                   }
    
    Sim_input = {'year' : 2020,                                 # the year, it is usedin order to create the schedule set, but it's almost useless
                 'first_day' : 7,                               # first day of the year
                 'tz' : 'Europe/Rome',                          # Thermal zone
                 'ts' : 1,                                      # Number of time-steps per hour
                 'hours' : 8760,                                # number of hours (actually only 8760 is implemented)
                 'model' : '2C',                                # model to be used (2C, 1C)
                 'DD_boundaries' : np.array(
                                 [[167,504], [4681,5017]],
                                 dtype = int),                  # Heating and Cooling Design Days Periods
                 'Time_to_regime' : 168,                        # Time needed to reach a regime condition for Design Days Calculation
                 'DesignDays_calc' : True,                      # Select 'YES' or 'NO' to calculate or not design days demand
                 'Plant_calc' : True,                           # Select 'YES' or 'NO' to calculate or not buildings plant
                 'OutputReport' : True,  
                 
                 'azSubdiv' : 8,                                # Azimuth subdivision
                 'hSubdiv' : 3,                                 # Tilt subdivision
                 'SolarCalc' : False,                           # Need of a Solar Calculation?
 
                 'Shading_calc' : True,                         # Select 'YES' or 'NO' to take into consideration the shading effect
                 'toll_az' : float(80),                         # Semi-tollerance on azimuth of shading surface [°]
                 'toll_dist' : float(100),                      # Tollerance on distance of shading surface [m]
                 'toll_theta' : float(80),                      # Semi-tollerance on position of the shading surfaces [°]
                 'R_f' : float(0)                               # Reduction factor of the direct solar radiation due to the shading effect [0-1]

                }
    
    
    # Urban Canyon Data
    
    UWG_data = {
        'UWG_calc' : bool(False), 
        'Area' : float(538000),
        'Diameter': float(800),
        'Perimeter': float(3280),
        'Ortogonal_length': float(800),
        'h_layer_inf': 2,
        'h_layer_sup': 150,
        'z_i_day': 1000,
        'z_i_night': 50,
        'z_i_m': 10,
        'z_i_ref': 150,
        'ExtWall_radiative_coef': float(5),
        'ExtWall_convective_coef': float(17),
        'ExtWall_emissivity': float(0.85),
        'ExtWall_reflection': float(0.15),
        'Average_U_windows': float(2),
        'Vegetation_density': float(0.16),
        'Road_radiative_coef': float(5),
        
        'Road_layers': np.array([
            [0.05, 1.9*10**6, 273 + 10],
            [0.2 , 2*10**6,   273 + 12],
            [1   , 1.4*10**6, 273 + 13],
            ], dtype = float),                          # first column thickness [m], second capacity [J/m3K], first row asphalt, second stones, third gravel 
            
        'Soil_layers': np.array([
            [3 , 1.4*10**6, 273 + 10],
            [3 , 1.4*10**6, 273 + 12],
            [3 , 1.4*10**6, 273 + 13],
            ], dtype = float),                          # first column thickness [m], second capacity [J/m3K]
        'T_deep': float(13+273),
        
        
        'cars_number_over_1000_px': float(39545000/60360),
        'motorcycles_number_over_1000_px': float(6896000/60360),
        'othervehicle_number_over_1000_px': float(5364000/60360),
        'emission_factor_cars_over_1000_px': float(	24.74 ),
        'emission_factor_motorcycles_over_1000_px': float(13.34),
        'emission_factor_othervehicle_over_1000_px': float(104.95),
        'vehicle_multiplying_factor': float(1),
        'distance_per_hour': float(64000),
        'population_density': float(2283.2),
        'first_day': Sim_input['first_day'],
        'Hw_week': np.array([0.9,0.5,0.9,2.1,3.3,4.2,5.8,7.1,6.3,5.8,5.4,5.4,5.8,6.3,6.7,7.1,7.5,5.8,4.2,3.3,2.5,1.7,1.3,0.9])/100,
        'Hw_end': np.array([1.7,1.2,1.,1.7,2.6,3.3,4.2,5,5.8,6.2,6.7,6.7,6.7,6.7,6.3,5.8,5,4.6,4.2,3.7,3.3,2.9,2.5,2.1])/100    
        }
    
    '''
    Reference H_traffic
    Italian vehicles number: http://www.opv.aci.it/WEBDMCircolante/
    emission factors
    https://link.springer.com/article/10.1007/s00704-008-0086-5/tables/4
    with 64 km/h
    
    '''

    city = Sim()
    city.set_input_from_dictionary(Input_files,Sim_input,UWG_data)
    city.preprocessing()
    city.city_creation()
    city.urban_shading()
    city.buildings_params_and_loads()
    city.plants_design_and_creation()
    city.simulation()
    city.output()
