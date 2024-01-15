# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 14:59:11 2024

@author: vivijac14771
"""

import numpy as np
from scipy.optimize import least_squares

class calibration():
       
    
    def __init__(self):
        '''This class contains all the functions needed to perform the calibration of the UBEM
        '''
        return
    
    
    
    def call_sim(self, x):
        '''Function to call EUReCA's ubem simulation with the new parameters. 
        It is called at each iteration of the calibration process. 
        

        Parameters
        ----------
        x : TYPE
            Solution array.

        Returns
        -------
        c_gas : TYPE
            Array with the gas consumption obtained with EUReCA's simulation.

        '''
        
        T_sp, f_w = self.unpack_solution(x)
        #
        ######################## To be replaced ##############################
        # chiama eureca qui
        # leggi risultati
        # c_gas è un vettore di N elementi costituito dal consumo annuale di gas di N edifici
        c_gas = T_sp + f_w # linea senza senso
        ######################################################################
        return c_gas
    
    
    def costfun(self, x, c_gas_meas):
        '''Function that calculates the vector of residuals between simulated and measured gas consumption.

        Parameters
        ----------
        x : TYPE
            Solution array.
        c_gas_meas : TYPE
            Array with the measured gas consumption.

        Returns
        -------
        res : TYPE
            Vector of residuals between measured and simulated gas consumption.
            
        Note
        ----
        The MSE is calculated inside the least squares from the vector of residuals

        '''
        
        c_gas = self.call_sim(x)
        res = c_gas - c_gas_meas
        return res
    
    
    def unpack_solution(self, x):
        '''Auxiliary function to read the solution array.

        Parameters
        ----------
        x : TYPE
            Solution array.

        Returns
        -------
        T_sp : TYPE
            Array with the indoor air temperature setpoints.
        f_w : TYPE
            Array with the multipliers of the external wall areas.

        '''
        N = int(len(x)/2.)
        T_sp = x[0:N-1]
        f_w  = x[N:2*N-1]
        
        return T_sp, f_w
    
    
    def calibrate(self, c_gas_meas, 
                  limits_setpoint = [18,22],
                  limits_fwalls = [0.75,1.25]
                  ):
        '''Function that performs the calibration of the buildings parameters 
           (temperature setpoints and external wall area multiplier) in order to 
           reduce the gap between measured and simulated gas consumption.
        

        Parameters
        ----------
        c_gas_meas : Numpy Array
            Array with the measured gas consumption.
        limits_setpoint : list, optional
            DESCRIPTION. The default is [18,22].
        limits_fwalls : list, optional
            DESCRIPTION. The default is [0.75,1.25].

        Returns
        -------
        x_opt : Numpy Array
            Optimal solution found by the calibration process.
        fmin : Float
            Value of the cost function at the solution.
        residuals : Numpy Array
            Vector of residuals at the solution.
        nfev : Int
            Number of function evaluations.

        '''

        N = len(c_gas_meas)
        # Set initial feasible solution
        T_sp_init = 20.*np.ones((N,))
        f_w_init = np.ones((N,))
        x0 = np.concatenate((T_sp_init, f_w_init))
        # Set domain
        T_sp_lb = limits_setpoint[0]*np.ones((N,))
        T_sp_ub = limits_setpoint[1]*np.ones((N,))
        f_w_lb = limits_fwalls[0]*np.ones((N,))
        f_w_ub = limits_fwalls[1]*np.ones((N,))
        lb = np.concatenate((T_sp_lb, f_w_lb))
        ub = np.concatenate((T_sp_ub, f_w_ub))
        # Call optimization to miminize objective function
        opt = least_squares(self.costfun, 
                            x0, 
                            bounds = (lb,ub), 
                            method = 'trf',
                            max_nfev = 50,
                            args = (c_gas_meas))
        # Read the solution
        x_opt = opt.x
        fmin = opt.cost
        residuals = opt.fun
        nfev = opt.nfev
        
        # Verbal description of the termination reason.
        print(opt.message)
        
        return x_opt, fmin, residuals, nfev
        
        
        