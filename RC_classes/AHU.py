'''IMPORTING MODULES'''

import sys
import numpy as np
from RC_classes.auxiliary_functions import wrn

#%% ---------------------------------------------------------------------------------------------------

def checkSatCond(T,x,p):
    '''
    Check Saturation Condition
    
    This function takes as inputs temperature [°C] and humidity ratio
    [kg_vap/kg_as] to check if a point is outside saturation conditions

    Parameters
    ----------
    T : float
        Temperature [°C]
    x : float
        Specific Humidity [kg_vap/kg_as].
    p : float
        Pressure [Pa].

    Returns
    -------
    Boolean, Saturation Pressure [Pa].
    
    '''
    
    # Check input data type
    
    if not isinstance(T, float):
        raise TypeError(f'ERROR input T is not an interger: T {T}')    
    if not isinstance(x, float):
        raise TypeError(f'ERROR input x is not an interger: x {x}')
    
    # Control input data quality
    if T < -15 or T > 60:
        wrn(f"WARNING CheckSatCond function, input temperature outside limit boundary [-15,60]: T {T}")
    if x < 0.0005 or x > 0.040:
        wrn(f"WARNING CheckSatCond function, input humidity outside limit boundary [0.0005,0.04]: x {x}")
        
    # Is or not outide the saturation condition? True/False
    pp = p*x/(0.622+x)
    if T < 0:
        psat = 610.5*np.exp((21.875*T)/(265.5+T))
    else:
        psat = 610.5*np.exp((17.269*T)/(237.3+T))        
    if pp - psat > 0.01:
        SatCond = False
    else:
        SatCond = True
    return SatCond, psat


#%%--------------------------------------------------------------------------------------------------- 
# AHU class
class AirHandlingUnit:
    '''
    This class manages the air handling unit.
    Some general variables are set as class variables while the __init__ method
    creates some useful vectors

    The method AHUCalc solves the AHU balance for a single time step. It takes:
        t: timestep [-]
        G_da_vent: mass flow rate [kg/s]
        AHUOnOff: AHU availabilty [-]
        AHUHUM: Humidistat availability [-]
        SensRecoveryEff: Sensible heat recovery efficiency [-]
        LatRecoveryEff: Latent heat recovery efficiency [-]
        OutAirRatio: Outdoor air ratio [-]
        T_ext: external temperature [°C]
        x_ext: external specific humidity [kg_v/kg_da]
        T_int: zone internal temperature [°C]
        x_int: zone internal specific humidity [kg_v/kg_da]
        T_sup: supply temperature [°C]
        x_sup: supplyspecific humidity [kg_v/kg_da]
    
    Methods:
        init
        AHUCalca   
    '''   
    
    # Class Variables
    cp_air = 1          # [kJ/(kg K)]
    p_atm = 101325      # [Pa]
    r_0 = 2501          # [kJ/kg]
    cpv = 1.875         # [kJ/(kg K)]
    T_out = 15          # [°C]                            
    
    def __init__(self,l, bdName = 'None'):
        
        '''
        Initializes the vectors of the AHU
        
        Parameters
            ----------
            l : int
                number of simulation time steps [-]
            bdName : string 
                name of the building.. just to print the warnings 
        Returns
        -------
        None.
        
        '''
        
        # Check input data type
        
        if not isinstance(l, int):
            raise TypeError(f'ERROR AHU inizialization, l must be an integer: l {l}')
        if not isinstance(bdName, str):
            raise TypeError(f'ERROR AHU inizialization, bdName must be an string: bdName {bdName}')
        
        # Inizialization
        
        self.AHUDemand = np.zeros(l)
        self.AHUDemand_sens = np.zeros(l)
        self.AHUDemand_lat = np.zeros(l)
        self.T_supAHU = np.zeros(l)
        self.bd_name = bdName
        
    
    def AHUCalc(self,t,G_da_vent,AHUOnOff,AHUHUM,Sens_Recovery_eff,Lat_Recovery_eff,OutAirRatio,T_ext,x_ext,T_int,x_int,T_sup,x_sup):
        
        '''
        Solution for the single time step of the Air Handling Unit
        
        Parameters:
            t: timestep [-]
            G_da_vent: mass flow rate [kg/s]
            AHUOnOff: AHU availabilty [-]
            AHUHUM: Humidistat availability [-]
            SensRecoveryEff: Sensible heat recovery efficiency [-]
            LatRecoveryEff: Latent heat recovery efficiency [-]
            OutAirRatio: Outdoor air ratio [-]
            T_ext: external temperature [°C]
            x_ext: external specific humidity [kg_v/kg_da]
            T_int: zone internal temperature [°C]
            x_int: zone internal specific humidity [kg_v/kg_da]
            T_sup: supply temperature [°C]
            x_sup: supplyspecific humidity [kg_v/kg_da]
            
        Returns:
            None
        '''       
        
        # Check input data type 
        
        if not isinstance(t, int):
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input t is not an interger: t {t}')
        if (not isinstance(G_da_vent, float)) or G_da_vent < 0:
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input G_da_vent is not a positive float: G_da_vent {G_da_vent}')
        if not AHUOnOff in [0,1,-1]:
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input AHUOnOff must be 0,1,-1: AHUOnOff {AHUOnOff}')
        if not isinstance(AHUHUM, bool):
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input AHUHUM is not a boolean: AHUHUM {AHUHUM}')
        if (not isinstance(Sens_Recovery_eff, float)) or Sens_Recovery_eff < 0. or Sens_Recovery_eff > 1.:
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input Sens_Recovery_eff must be included in the range 0-1: Sens_Recovery_eff {Sens_Recovery_eff}')
        if (not isinstance(Lat_Recovery_eff, float)) or Lat_Recovery_eff < 0. or Lat_Recovery_eff > 1.:
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input Lat_Recovery_eff must be included in the range 0-1: Lat_Recovery_eff {Lat_Recovery_eff}')
        if (not isinstance(OutAirRatio, float)) or OutAirRatio < 0. or OutAirRatio > 1.:
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input OutAirRatio must be included in the range 0-1: OutAirRatio {OutAirRatio}')
        if not isinstance(T_ext, float):
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input T_ext is not an interger: T_ext {T_ext}')
        if not isinstance(T_int, float):
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input T_int is not an interger: T_int {T_int}')
        if not isinstance(T_sup, float):
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input T_sup is not an interger: T_sup {T_sup}')    
        if not isinstance(x_ext, float):
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input x_ext is not an interger: x_ext {x_ext}')
        if not isinstance(x_int, float):
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input x_int is not an interger: x_int {x_int}')
        if not isinstance(x_sup, float):
            raise TypeError(f'ERROR AHUCalc, bd {self.bd_name}, time step {t}, input x_sup is not an interger: x_sup {x_sup}')
        
        # Control input data quality
        
        if T_ext < -15 or T_ext > 60:
            wrn(f"WARNING AHU class, AHUCalc method, bd {self.bd_name}, time step {t}, input temperature outside limit boundary [-15,60]: T_ext {T_ext}")
        if T_int < -15 or T_int > 60:
            wrn(f"WARNING AHU class, AHUCalc method, bd {self.bd_name}, time step {t}, input temperature outside limit boundary [-15,60]: T_int {T_int}")
        if T_sup < -15 or T_sup > 60:
            wrn(f"WARNING AHU class, AHUCalc method, bd {self.bd_name}, time step {t}, input temperature outside limit boundary [-15,60]: T_sup {T_sup}")
        if x_ext < 0.0005 or x_ext > 0.040:
            wrn(f"WARNING AHU class, AHUCalc method, bd {self.bd_name}, time step {t}, input humidity outside limit boundary [0.0005,0.04]: x_ext {x_ext}")
        if x_int < 0.0005 or x_int > 0.040:
            wrn(f"WARNING AHU class, AHUCalc method, bd {self.bd_name}, time step {t}, input humidity outside limit boundary [0.0005,0.04]: x_int {x_int}")
        if x_sup < 0.0005 or x_sup > 0.040:
            wrn(f"WARNING AHU class, AHUCalc method, bd {self.bd_name}, time step {t}, input humidity outside limit boundary [0.0005,0.04]: x_sup {x_sup}")       
            
        # Set some input variables
        
        self.T_supAHU[t] = T_sup
        self.x_sup = x_sup
        
        
        # Saturation conditions check
        SatCond, psat = checkSatCond(T_ext,x_ext,self.p_atm)
        if SatCond == False:
            x_ext = 0.99*0.622*psat/(self.p_atm-0.99*psat)
            if AHUOnOff == 0:
                self.x_sup = x_ext
        SatCond, psat = checkSatCond(T_int,x_int,self.p_atm)
        if SatCond == False:
            wrn('ERROR  AHUCalc method, bd {self.bd_name}, time step {t}, Zone conditions outside saturation limit')
        SatCond, psat = checkSatCond(self.T_supAHU[t],self.x_sup,self.p_atm)
        if SatCond == False:
            wrn('ERROR:  AHUCalc method, bd {self.bd_name}, time step {t}, Supply conditions outside saturation limit')
        
        
        # Pre-processing on Heat Recovery and Mixer
        self.h_ext = self.cp_air*T_ext + (self.r_0 + self.cpv*T_ext)*x_ext
        self.h_z = self.cp_air*T_int + (self.r_0 + self.cpv*T_int)*x_int
        
        
        # Heat Recovery Bypass in condition of free heating or cooling
        if (T_int-T_ext)*AHUOnOff > 0:
            self.T_hr = T_ext + Sens_Recovery_eff*(T_int-T_ext)
            self.T_out = T_int - Sens_Recovery_eff*(T_int-T_ext)
        else:
            self.T_hr = T_ext
            self.T_out = T_int
            OutAirRatio = 1
        if (x_int-x_ext)*AHUOnOff > 0:
            self.x_hr = x_ext + Lat_Recovery_eff*(x_int-x_ext)
        else:
            self.x_hr = x_ext
        self.h_hr = self.cp_air*self.T_hr +(self.r_0 + self.cpv*self.T_hr)*self.x_hr
        
        
        # Mixer
        self.h_mix = OutAirRatio*self.h_hr + (1 - OutAirRatio)*self.h_z
        self.x_mix = OutAirRatio*self.x_hr + (1 - OutAirRatio)*x_int
        self.T_mix = (self.h_mix - self.r_0*self.x_mix)/(self.cp_air + self.cpv*self.x_mix)
        
        
        # BATTERIES DEMAND CALCULATION
        
        # SENSIBLE AND LATENT CONTROL MODE
        if AHUHUM == True:
        
            # Heating mode
            if AHUOnOff == 1:
                                
                
                # Adiabatic Saturator
                self.x_as = self.x_sup
                p_as = self.p_atm*self.x_as/(0.622+self.x_as)
                if p_as >= 610.5:
                    self.T_as = 237.3*np.log(p_as/610.5)/(17.269-np.log(p_as/610.5))
                else:
                    self.T_as = 265.5*np.log(p_as/610.5)/(21.875-np.log(p_as/610.5))
                self.h_as = self.cp_air*self.T_as + (self.r_0 + self.cpv*self.T_as)*self.x_as
                
                
                # Pre-Heater control
                if self.h_mix < self.h_as:
                    self.h_ph = self.h_as
                    self.x_ph = self.x_mix
                    self.T_ph = (self.h_ph-self.r_0*self.x_ph)/(1.006+self.cpv*self.x_ph)
                else:
                    # In this case the Pre-Heater is by-passed
                    self.h_ph = self.h_mix
                    self.x_ph = self.x_mix
                    self.T_ph = self.T_mix
                self.h_as = self.h_ph
                
                
                # Humidification need check
                if self.x_ph > self.x_sup:
                    self.x_as = self.x_ph
                    self.x_sup = self.x_as
                    self.T_as = (self.h_as - self.r_0*self.x_as)/(self.cp_air+self.cpv*self.x_as)
                
                
                # Controlling T_supAHU > T_as
                if self.T_supAHU[t] < self.T_as:
                    self.T_supAHU[t] = self.T_as
                self.h_sup = self.cp_air*self.T_supAHU[t]+(self.r_0+self.cpv*self.T_supAHU[t])*self.x_sup
                
                
                # Batteries Demand [kW]
                # Pre-Heater
                preh_Dem = G_da_vent*(self.h_ph-self.h_mix)
                preh_Dem_sens = G_da_vent*self.cp_air*(self.T_ph-self.T_mix)
                preh_Dem_lat = G_da_vent*self.r_0*(self.x_ph-self.x_mix)
                # Adiabatic Saturator
                sat_Dem = G_da_vent*(self.h_as-self.h_ph)
                sat_Dem_lat = G_da_vent*self.r_0*(self.x_as-self.x_ph)
                sat_Dem_sens = sat_Dem-sat_Dem_lat
                # Post-Heater
                posth_Dem = G_da_vent*(self.h_sup-self.h_as)
                posth_Dem_sens = G_da_vent*self.cp_air*(self.T_supAHU[t]-self.T_as)
                posth_Dem_lat = G_da_vent*self.r_0*(self.x_sup - self.x_as)
                # Total Demand
                self.AHUDemand[t] = preh_Dem + sat_Dem + posth_Dem
                self.AHUDemand_sens[t] = preh_Dem_sens + sat_Dem_sens + posth_Dem_sens
                self.AHUDemand_lat[t] = preh_Dem_lat + sat_Dem_lat + posth_Dem_lat
                
                
                # Cooling mode
            elif AHUOnOff == -1:
                
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
                    if self.T_supAHU[t] <= self.T_mix:
                        self.T_de = self.T_supAHU[t]
                    else:
                        self.T_de = self.T_mix
                        self.T_supAHU[t] = self.T_de
                self.h_de = self.cp_air*self.T_de+(self.r_0+self.cpv*self.T_de)*self.x_de
                self.h_sup = self.cp_air*self.T_supAHU[t]+(self.r_0+self.cpv*self.T_supAHU[t])*self.x_sup
                
                
                # Batteries Demand [kW]
                # Dehumidificator
                deu_Dem = G_da_vent*(self.h_de - self.h_mix)
                deu_Dem_sens = G_da_vent*self.cp_air*(self.T_de - self.T_mix)
                deu_Dem_lat = G_da_vent*self.r_0*(self.x_de - self.x_mix)
                # Post-Heater
                posth_Dem = G_da_vent*(self.h_sup - self.h_de)
                posth_Dem_sens = G_da_vent*self.cp_air*(self.T_supAHU[t] - self.T_de)
                posth_Dem_lat = G_da_vent*self.r_0*(self.x_sup - self.x_de)
                # Total Demand
                self.AHUDemand[t] = deu_Dem + posth_Dem
                self.AHUDemand_sens[t] = deu_Dem_sens + posth_Dem_sens
                self.AHUDemand_lat[t] = deu_Dem_lat + posth_Dem_lat
                

            elif AHUOnOff == 0:
                
                # Plant Off
                self.T_hr, self.T_mix, self.T_supAHU[t] = T_ext, T_ext, T_ext
                self.x_hr, self.x_mix, self.x_sup = x_ext, x_ext, x_ext
                self.h_hr, self.h_mix, self.h_sup = self.h_ext, self.h_ext, self.h_ext
                
                
                # Batteries Demand [kW]
                self.AHUDemand[t] = 0
                self.AHUDemand_sens[t] = 0
                self.AHUDemand_lat[t] = 0
            
            else:
                print('AHUOnOff value not allowed')
        
        
        # SENSIBLE CONTROL MODE
        if AHUHUM == False:
            
            
            # Heating mode
            if AHUOnOff == 1:
                
                
                # Pre-Heater and Adiabatic Saturator doesn't exist!!
                self.h_ph = self.h_mix
                self.x_ph = self.x_mix
                self.T_ph = self.T_mix
                self.h_as = self.h_ph
                self.x_as = self.x_ph
                self.T_as = self.T_ph
                
                
                # Post-Heater
                if self.T_supAHU[t]<self.T_as:
                    self.T_supAHU[t] = self.T_as
                
                self.x_sup = self.x_as
                self.h_sup = self.cp_air*self.T_supAHU[t]+(self.r_0+self.cpv*self.T_supAHU[t])*self.x_sup
                
                
                # Batteries Demand [kW]
                self.AHUDemand[t] = G_da_vent*self.cp_air*(self.T_supAHU[t]-self.T_as)
                self.AHUDemand_sens[t] = self.AHUDemand[t]
                self.AHUDemand_lat[t] = 0

    
            elif AHUOnOff == -1:
                # Cooling mode
                
                if self.T_supAHU[t] > self.T_mix:
                    self.T_supAHU[t] = self.T_mix
                    self.x_sup = self.x_mix
                    self.h_sup = self.h_mix
                else:
                    SatCond, psat = checkSatCond(self.T_supAHU[t],self.x_mix,self.p_atm)
                    if SatCond == False:
                        # print('Info: supply temperature is too low -- changed')

                        if self.T_supAHU[t] < 0:
                            p_sat = 610.5*np.exp((21.875*self.T_supAHU[t])/(265.5+self.T_supAHU[t]))
                        else:
                            p_sat = 610.5*np.exp((17.269*self.T_supAHU[t])/(237.3+self.T_supAHU[t]))
                        self.x_sup = 0.99*0.622*p_sat/(self.p_atm-0.99*p_sat)
                        self.h_sup = self.cp_air*self.T_supAHU[t]+(self.r_0+self.cpv*self.T_supAHU[t])*self.x_sup
                    else:
                        self.x_sup= self.x_mix
                        self.h_sup = self.cp_air*self.T_supAHU[t]+(self.r_0+self.cpv*self.T_supAHU[t])*self.x_sup
                
                
                #  Batteries Demand [kW]
                self.AHUDemand[t] = G_da_vent*(self.h_sup - self.h_mix)
                self.AHUDemand_sens[t] = self.AHUDemand[t]
                self.AHUDemand_lat[t] = 0
                
                
            elif AHUOnOff == 0:
                # Plant OFF
                
                self.T_hr, self.T_mix, self.T_supAHU[t] = T_ext, T_ext, T_ext
                self.x_hr, self.x_mix, self.x_sup = x_ext, x_ext, x_ext
                self.h_hr, self.h_mix, self.h_sup = self.h_ext, self.h_ext, self.h_ext
                
                # Batteries Demand [kW]
                self.AHUDemand[t] = 0
                self.AHUDemand_sens[t] = 0
                self.AHUDemand_lat[t] = 0
            
            else:
                sys.exit('AHUOnOff value not allowed at time step: '+str(t))
                
#%%--------------------------------------------------------------------------------------------------- 
#%%

if __name__ == '__main__':
    p_atm = 101325
    SatCond, psat = checkSatCond(12,0.009,p_atm)
    if SatCond == False:
        print('ERROR: Supply conditions outside saturation limit')
