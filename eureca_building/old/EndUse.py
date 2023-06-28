'''IMPORTING MODULES'''

import sys
import pandas as pd
import os
import numpy as np
from RC_classes.auxiliary_functions import wrn

#%% ---------------------------------------------------------------------------------------------------
#%% Useful functions to create the schedule archetypes

def loadArchetype(path,timeIndex,ts):
    '''
    Archetype loading in case you use loadSchedComp
    This function takes the path of the excel file constining the schedules and loads it
    Works for the yearly schedule (see file ScheduleComp.xlsx in /Input/
    
    Parameters
    ----------
    path : string
        Path containing the string of the file_schedule.xlsx
    timeIndex : np array of int
        This is the array containing the index of the simulation time steps .
    ts : int
        Number of time steps per hour.

    Returns
    -------
    archetypes: dictionary with archetype_key/Archetype(object) data
    '''
    
    # Check input data type  

    if not isinstance(path, str):
        raise TypeError(f'ERROR input path is not a string: path {path}') 
    if not isinstance(timeIndex, np.ndarray):
        raise TypeError(f'ERROR input timeIndex is not a np.array: timeIndex {timeIndex}') 
    if not isinstance(ts, int):
        raise TypeError(f'ERROR input ts is not an integer: ts {ts}')         
    
    # Control input data quality
    
    if ts > 4:
        wrn(f"WARNING loadSimpleArchetype function, input ts is higher than 4, this means more than 4 time steps per hours were set: ts {ts}")
         
    try:
        sched = pd.read_excel(path,sheet_name="Schedule",header=2,index_col=[0]).set_index(timeIndex)
        arch = pd.read_excel(path,sheet_name="Archetype",header=1,index_col=[0])
    except FileNotFoundError:
        raise FileNotFoundError(f'ERROR Failed to open the schedule xlsx file {path}... Insert a proper path')

    # Reset some indexes set the right index
    
    sched = sched.reset_index().drop(columns=["Time"])
    sched.index=timeIndex
    
    # Creation of the archetypes' dictionary
        
    archetypes = dict()
    for i in arch.index:
        if i != 'Archetype':
            archetypes[i] = Archetype(i)
            archetypes[i].loadSchedComp(arch.loc[i],sched)
            archetypes[i].rescale_df(ts)
            archetypes[i].create_np()
    return archetypes

#%%

def loadSimpleArchetype(path,timeIndex,first_day = 1,ts = 1, PlantDays = [2520,3984,6192,6912]):
    '''
    Archetype loading in case you use loadSchedSemp
    This function takes the path of the excel file constining the schedules and loads it
    Works for the daily schedule (see file ScheduleSemp.xlsx in /Input/
    
    Parameters
    ----------
    path : string
        Path containing the string of the file_schedule.xlsx
    timeIndex : np array of int
        This is the array containing the index of the simulation time steps .
    first_day : int
        First day of the week (1 monday, 2 tuesday, ..... 7 sunday)
    ts : int
        Number of time steps per hour.
    PlantDays : list
        List of integers that sets the 4 time steps of heating and cooling season start and stops
        [last heating timestep,
         first cooling timestep,
         last cooling timstep,
         fist heating time step]

    Returns
    -------
    archetypes: dictionary with archetype_key/Archetype(object) data
    '''
    
    # Check input data type  
    
    if not isinstance(path, str):
        raise TypeError(f'ERROR input path is not a string: path {path}') 
    if not isinstance(timeIndex, np.ndarray):
        raise TypeError(f'ERROR input timeIndex is not a np.array: timeIndex {timeIndex}') 
    if not isinstance(first_day, int):
        raise TypeError(f'ERROR input first_day is not an integer: first_day {first_day}')   
    if not isinstance(ts, int):
        raise TypeError(f'ERROR input ts is not an integer: ts {ts}')
    if not isinstance(PlantDays, list):
        raise TypeError(f'ERROR input PlantDays is not a list: PlantDays {PlantDays}\nThis parameter must be a list of integers (default [2520,3984,6192,6912])')
    if not isinstance(PlantDays[0], int) or not isinstance(PlantDays[1], int) or not isinstance(PlantDays[2], int) or not isinstance(PlantDays[3], int):
        raise TypeError(f'ERROR input PlantDays in not a list of integers: PlantDays {PlantDays}')   
    
    # Control input data quality
    
    if first_day > 7 or first_day < 1:
            wrn(f"WARNING loadSimpleArchetype function, input fisrt_day should be in the range [0,7]: first_day {first_day}")
    if ts > 4:
            wrn(f"WARNING loadSimpleArchetype function, input ts is higher than 4, this means more than 4 time steps per hours were set: ts {ts}")
        
    try:
        ex = pd.ExcelFile(path)
    except FileNotFoundError:
        raise FileNotFoundError(f'ERROR Failed to open the schedule xlsx file {path}... Insert a proper path')
    
    # Creation of the archetypes' dictionary
    
    archetypes = dict()
    last_H_day = PlantDays[0]
    first_C_day = PlantDays[1]
    last_C_day = PlantDays[2]
    first_H_day = PlantDays[3]
    
    # read archetype names from the excel sheet
    archetype_list = pd.read_excel(path,sheet_name='GeneralData',header=[0],index_col=[0], skiprows = 1)
    names = archetype_list.index
    
    for use in names:
        archetypes[use] = Archetype(use)
        year = pd.DataFrame()
        schedule = pd.read_excel(path,sheet_name=use,header=[0,1,2],index_col=[0])
        schedule = schedule.iloc[1:25]*schedule.loc['Valore nominale']
        week = pd.concat([schedule['Weekday']]*5+[schedule['Weekend']]*2)
        
        # The heating cooling availabilty is create considering the PlantDays variable
        
        first_week = week.iloc[(first_day-1)*24 :]
        year = pd.concat([first_week]+[week]*53).iloc[:8760].set_index(timeIndex)
        year.loc[last_H_day:first_C_day,('Plant Availability','[-]')]=year.loc[last_H_day:first_C_day,('Plant Availability','[-]')]*0
        year.loc[first_C_day:last_C_day,('Plant Availability','[-]')]=year.loc[first_C_day:last_C_day,('Plant Availability','[-]')]*-1
        year.loc[last_C_day:first_H_day,('Plant Availability','[-]')]=year.loc[last_C_day:first_H_day,('Plant Availability','[-]')]*0
        archetypes[use].loadSchedSemp(year, archetype_list.loc[use])
        archetypes[use].rescale_df(ts)
        archetypes[use].create_np()

    return archetypes

#%%--------------------------------------------------------------------------------------------------- 
#%% Archetype class

class Archetype:
    '''
    This class manages the end use with its schedules

    init method: sets only the name and creats a Dataframe:
        name: a string with the name
        
    loadSchedComp: loads the complex (annual) method for schedules:
        arch: list containing the name of the schedules of the archetype
        sched: dictionary with all the annual schedules
        
    loadSchedSemp: loads the simple method:
        takes a yearly dataframe with the archetype schedules
        
    rescale_df: rescales the dataframe with respect to the ts parameter (time steps per hour)
    
    create_np: convert every schedule in numpy arrays
        
    Methods:
        init
        loadSchedComp
        loadSchedSemp
        rescale_df
        create_np
    ''' 
    
    def __init__(self,name):
    
        '''
        Initializes the vectors of the AHU
        
        Parameters
            ----------
            name : string
                name of the archetype
                
        Returns
        -------
        None.
        
        '''    
        
        # Check input data type
        
        if not isinstance(name, str):
            raise TypeError(f'ERROR Archetype initialization, name must be a string: name {name}')
        
        # Inizialization
        
        self.name = name
        self.sched_df = pd.DataFrame()
        self.scalar_data = {}
        '''
        Take care of units!!!
        All formula referes to the complex schedule file units, 
        Simple schedule file has different units which are converted in the loadSchedSemp method
        (vapuor flow rates)
        '''
    
    def loadSchedComp(self,arch,sched):
    
        '''
        Used for ScheduleComp.xlsx Excel file
        
        Parameters
            ----------
            arch : pandas dataframe
                This must include all the schedules' keys
            
            sched : pandas dataframe
                This dataframe includes all the yearly schedules
                
        Returns
        -------
        None.
        
        '''
            
        # Each schedule is set to a different attribute of the class
        
        try:
            self.sched_df['appliances'] = sched[arch['Appliances']]
            self.sched_df['lighting'] = sched[arch['Lighting']]
            self.sched_df['people'] = sched[arch['People (Sensible)']]
            self.sched_df['vapour'] = sched[arch['Vapour']]
            self.sched_df['heatingTSP'] = sched[arch['HeatTSP']]
            self.sched_df['coolingTSP'] = sched[arch['CoolTSP']]
            self.sched_df['HeatingRHSP'] = sched[arch['HeatRHSP']]
            self.sched_df['CoolingRHSP'] = sched[arch['CoolRHSP']]
            self.sched_df['ventFlowRate'] = sched[arch['VentFlowRate']]
            self.sched_df['infFlowRate'] = sched[arch['InfFlowRate']]
            self.sched_df['plantOnOffSens'] = sched[arch['PlantONOFFSens']]
            self.sched_df['plantOnOffLat'] = sched[arch['PlantONOFFLat']]
            self.sched_df['AHUOnOff'] = sched[arch['AHUONOFF']]
            self.sched_df['AHUHUM'] = sched[arch['AHUHUM']]
            self.sched_df['AHUTSupp'] = sched[arch['AHUTSupp']]
            self.sched_df['AHUxSupp'] = sched[arch['AHUxSupp']]
            
            self.scalar_data['conFrac'] = float(arch['ConvFrac'])
            self.scalar_data['AHUHUM'] = bool(arch['AHUHum'])
            self.scalar_data['sensRec'] = float(arch['SensRec'])
            self.scalar_data['latRec'] = float(arch['LatRec'])
            self.scalar_data['outdoorAirRatio'] = float(arch['OutAirRatio'])
        except KeyError:
            raise KeyError(f'ERROR Archetype object {self.name}: can not find all schedules')
        
    def loadSchedSemp(self,year_df,supplementary_data):
        '''
        Used for ScheduleSemp.xlsx Excel file
        
        Parameters
            ----------
            arch : pandas dataframe
                This includes the yearly schedule of a single archetype
            supplementary_data : pd.series
                This series includes some additional data about the archetype (Sensible and 
                                                                               Latent AHU recovery, 
                                                                               Convective fraction of internal gains)
            
        Returns
        -------
        None.
        
        '''    
    
        # Each schedule is set to a different attribute of the class
    
        try:
            self.sched_df['appliances'] = year_df['Appliances','[W/m²]']
            self.sched_df['lighting'] = year_df['Lighting','[W/m²]']
            self.sched_df['people'] = year_df['Occupancy (Sensible)','[W/m²]']
            self.sched_df['vapour'] = year_df['Vapour FlowRate','[g/(m² s)]']/1000 #--> kg conversion (a)
            self.sched_df['heatingTSP'] = year_df['HeatSP','[°C]']
            self.sched_df['coolingTSP'] = year_df['CoolSP','[°C]']
            self.sched_df['HeatingRHSP'] = year_df['HumSP','[%]']/100
            self.sched_df['CoolingRHSP'] = year_df['DehumSP','[%]']/100
            self.sched_df['ventFlowRate'] = year_df['Ventilation FlowRate','[m³/(s m²)]']
            self.sched_df['infFlowRate'] = year_df['Infiltration FlowRate','[Vol/h]']
            self.sched_df['plantOnOffSens'] = year_df['Plant Availability','[-]']
            self.sched_df['plantOnOffLat'] = year_df['Plant Availability','[-]']
            self.sched_df['AHUOnOff'] = year_df['Plant Availability','[-]']
        except KeyError:
            raise KeyError(f'ERROR Archetype object {self.name}: can not find all schedules')
        
        # Following parameters are not set in the Excel file
        
        self.sched_df['AHUTSupp'] = pd.Series([22.]*8760,index=self.sched_df['appliances'].index)
        self.sched_df['AHUxSupp'] = pd.Series([0.005]*8760,index=self.sched_df['appliances'].index)
        self.sched_df['AHUxSupp'].iloc[3984:6192] = 0.007
        
        try:
            self.scalar_data['conFrac'] = float(supplementary_data.loc['ConvFrac'])
            self.scalar_data['AHUHUM'] = bool(supplementary_data.loc['AHUHum'])
            self.scalar_data['sensRec'] = float(supplementary_data.loc['SensRec'])
            self.scalar_data['latRec']= float(supplementary_data.loc['LatRec'])
            self.scalar_data['outdoorAirRatio'] = float(supplementary_data.loc['OutAirRatio'])
        except KeyError:
            raise KeyError(f"ERROR Loading end use {self.name}. GeneralData does not have the correct columns names: ConvFrac, AHUHum, SensRec, LatRec, OutAirRatio")
        except ValueError:
            raise ValueError(f"""ERROR 
                             Loading end use {self.name}. GeneralData
                             I'm not able to parse the General data. 
                                 ConvFrac should be a float {supplementary_data['ConvFrac']}
                                 AHUHum should be a boolean {supplementary_data['AHUHum']}
                                 SensRec should be a float {supplementary_data['SensRec']}
                                 LatRec  should be a float {supplementary_data['LatRec']}
                                 OutAirRatio   should be a float {supplementary_data['OutAirRatio']}
                             """)
        
        # Check the quality of input data
        if not 0. <= self.scalar_data['conFrac'] <= 1.:
            wrn(f"WARNING Loading end use {self.name}. Convective fraction of the heat gain outside boundary condition [0-1]: ConvFrac {self.scalar_data['conFrac']}")
        if not 0. <= self.scalar_data['sensRec'] <= 1.:
            wrn(f"WARNING Loading end use {self.name}. Sensible recovery of the AHU outside boundary condition [0-1]: sensRec {self.scalar_data['sensRec']}")
        if not 0. <= self.scalar_data['latRec'] <= 1.:
            wrn(f"WARNING Loading end use {self.name}. Latent recovery of the AHU outside boundary condition [0-1]: sensRec {self.scalar_data['latRec']}")
        if not 0. <= self.scalar_data['outdoorAirRatio'] <= 1.:
            wrn(f"WARNING Loading end use {self.name}. Outdoor air ratio of the AHU outside boundary condition [0-1]: outdoorAirRatio {self.scalar_data['outdoorAirRatio']}")
        
    def rescale_df(self,ts):
        '''
        rescale the archetype dataframe with respect to the number of time steps per hour
        
        Parameters
            ----------
            ts : int
                Number of time steps per hour
                
        Returns
        -------
        None.        
        '''   
        
        # Check input data type
        
        if not isinstance(ts, int):
            raise TypeError(f'ERROR input ts is not an integer: ts {ts}')         
    
        # Rescale 
        
        m = str(60/ts) + 'min'
        time = pd.date_range('2020-01-01', periods=8760, freq='1h')
        self.sched_df.set_index(time, inplace=True) 
        self.sched_df = self.sched_df.resample(m).pad()                               # Steps interpolation Resampling
        #Boundary = Boundary_0.resample(str(ts)+'S').interpolate(method='linear')        # Linear interpolation Resampling 
        # There are several upsample methods: pad(), bfill(), mean(), interpolate(), apply custum function...
        # For info look:
        #   https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.resample.html
        #   https://machinelearningmastery.com/resample-interpolate-time-series-data-python/
        
    def create_np(self):
    
        '''
        Creates the np.array attributes from the self.sched_df dataframe
        
        Parameters
            ----------
                
        Returns
        -------
        None.        
        '''      
    
        self.appliances = self.sched_df['appliances'].to_numpy(dtype = np.float_)
        self.lighting = self.sched_df['lighting'].to_numpy(dtype = np.float_)
        self.people = self.sched_df['people'].to_numpy(dtype = np.float_)
        self.vapour = self.sched_df['vapour'].to_numpy(dtype = np.float_)
        self.heatingTSP = self.sched_df['heatingTSP'].to_numpy(dtype = np.float_)
        self.coolingTSP = self.sched_df['coolingTSP'].to_numpy(dtype = np.float_)
        self.HeatingRHSP = self.sched_df['HeatingRHSP'].to_numpy(dtype = np.float_)
        self.CoolingRHSP = self.sched_df['CoolingRHSP'].to_numpy(dtype = np.float_)
        self.ventFlowRate = self.sched_df['ventFlowRate'].to_numpy(dtype = np.float_)
        self.infFlowRate = self.sched_df['infFlowRate'].to_numpy(dtype = np.float_)
        self.plantOnOffSens = self.sched_df['plantOnOffSens'].to_numpy(dtype = np.int_)
        self.plantOnOffLat = self.sched_df['plantOnOffLat'].to_numpy(dtype = np.int_)
        self.AHUOnOff = self.sched_df['AHUOnOff'].to_numpy(dtype = np.int_)
        #self.AHUHUM = self.sched_df['AHUHUM'].to_numpy(dtype = np.bool_)
        self.AHUTSupp = self.sched_df['AHUTSupp'].to_numpy(dtype = np.float_)
        self.AHUxSupp = self.sched_df['AHUxSupp'].to_numpy(dtype = np.float_)
        #self.conFrac = self.sched_df['conFrac'].to_numpy(dtype = np.float_)
        #self.AHUHumidistat = self.sched_df['AHUHum']
        #self.sensRec = self.sched_df['sensRec'].to_numpy(dtype = np.float_)
        #self.latRec = self.sched_df['latRec'].to_numpy(dtype = np.float_)
        #self.outdoorAirRatio = self.sched_df['outdoorAirRatio'].to_numpy(dtype = np.float_)        

#%%--------------------------------------------------------------------------------------------------- 
#%%
'''
TEST METHOD
'''
if __name__=='__main__':
    schedulepath = os.path.join('..','Input','Schedule_ICA-RC.xlsx')
    time=pd.date_range('2012-01-01', periods=8760, freq='1H')
    archetypes = loadArchetype(schedulepath,time)
    schedulepath2 = os.path.join('..','Input','Schedule_template.xlsx')
    archetypes2 = loadSimpleArchetype(schedulepath2,time) 