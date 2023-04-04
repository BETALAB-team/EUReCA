# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 11:42:23 2021

@author: ROBERTO
"""

import matplotlib.pyplot as plt
import pandas as pd
import xlrd
import matplotlib.mlab as mlab
import math
import time
import numpy as np

def distrEventi(n,x,pdf):
    y_guess=np.random.rand(int(n))
    cdf=np.cumsum(pdf)
    [cdf,index]=np.unique(cdf,1)

    x_event=np.interp(y_guess,cdf,x[index])
    time_event=np.round(x_event)
    return time_event

class DHW():
    
    
    total_days=365 #number of days for the simulation
    nuses=4 #number of uses
    
    vol_aver_drawoff=np.array([15/60,360/60,480/60,480/60],dtype=float) #Liters/min mu
    vol_desv_drawoff=np.array([60/60,120/60,120/60,120/60],dtype=float) #Liters/min sigma
    time_aver_drawoff=np.array([1,10,10,5],dtype=float) #minutes
    proportion_event=np.array([0.14,0.36,0.10,0.40],dtype=float) #probability in percentile of 1 for each event to happen from the total daily volume
    dist=np.array([[0.01,0.01,0.01,0.01,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.01],
                   [0.01,0.01,0.01,0.01,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.01],
                   [0,0.0000001,0.0000002,0.0000003,0.0000004,0.0000005,0.0103092783505155,0.0206185567010309,0.0257731958762887,0.0309278350515464,0.0463917525773196,0.0515463917525773,0.0515463917525773,0.0515463917525773,0.0515463917525773,0.0515463917525773,0.0515463917525773,0.134020618556701,0.22680412371134,0.134020618556701,0.0309278350515464,0.0206185567010309,0.0103092783505155,0],
                   [0,0.00001,0.00002,0.00003,0.0380952380952381,0.142857142857143,0.238095238095238,0.142857142857143,0.0380952380952381,0.019047619047619,0.019047619047619,0.019047619047619,0.019047619047619,0.019047619047619,0.019047619047619,0.019047619047619,0.019047619047619,0.0285714285714286,0.0761904761904762,0.0761904761904762,0.0285714285714286,0.019047619047619,0.019047619047619,0]],dtype=float).T    
   
    
    def __init__(self, volume_unit, numunits, time_step):
        
                
        vol_total_drawoff=self.vol_aver_drawoff*self.time_aver_drawoff
        
        vol_total_drawoff_use=np.zeros((self.nuses))
        
        
        for j in range(self.nuses):
            vol_total_drawoff_use[j]=volume_unit*self.proportion_event[j]
                
        draw_offs=np.zeros((self.nuses))
        draw_offs_dec=np.zeros((self.nuses)) 
        
        for j in range(self.nuses):
            draw_offs_dec[j]=np.abs(vol_total_drawoff_use[j]/self.vol_aver_drawoff[j]-int(vol_total_drawoff_use[j]/self.vol_aver_drawoff[j]))
            if draw_offs_dec[j]>=0.5:
                draw_offs[j]=math.ceil(vol_total_drawoff_use[j]/self.vol_aver_drawoff[j])
            else:
                draw_offs[j]=np.round(vol_total_drawoff_use[j]/self.vol_aver_drawoff[j])
                
        time_steps_hour=60/time_step
        time_steps_day=24*time_steps_hour              

        array_time_step=np.zeros(4)
        for i in range(4):
            array_time_step[i]=time_step
        
        n_max=array_time_step/self.time_aver_drawoff
        for i in range(len(n_max)):
            n_max[i]=math.ceil(n_max[i])
        
        Volume_use_daily_array=np.zeros((self.total_days,int(time_steps_day)))
        Volume_use_arrayb0=np.zeros((int(self.total_days*time_steps_day)))
        dist_use=np.zeros((24,12))
        Volume_use_array=np.zeros((int(time_steps_day*self.total_days),self.nuses))
        Volume_use_sum=np.zeros((1,self.nuses))
        Volume_aver_drawoff_final=np.zeros(self.nuses)
        Volume_desv_drawoff_final=np.zeros((1,self.nuses))
        Volume_use_time=np.zeros(np.int(time_steps_day))
        
        inicio=time.time()
        
        
        
        Volume_use_unit=np.zeros((int(self.total_days*time_steps_day),numunits))

        for units0 in range(numunits):
    
            for use in range(self.nuses):
    
                dist_use=np.repeat(self.dist[:,use], 12, axis=0)
                dist_use_t=np.reshape(dist_use,(288,1))/time_steps_hour
    
                vol_aver_drawoff_use=self.vol_aver_drawoff[use]
                vol_desv_drawoff_use=self.vol_desv_drawoff[use]
                time_aver_drawoff_use=self.time_aver_drawoff[use]
                draw_offs_use=int(draw_offs[use])
                n_max1=n_max[use]
    
                for day in range(self.total_days):
    
                    time_event=distrEventi(draw_offs_use,np.arange(0,time_steps_day),dist_use_t)
    
                    if use==0:
                        m=vol_aver_drawoff_use
                        v=vol_desv_drawoff_use
                        mu=np.log((m**2)/np.sqrt(v+m**2))
                        sigma=np.sqrt(np.log(v/(m**2)+1))
                        Flow_rated=np.random.lognormal(mu,sigma,int(draw_offs_use))
           
                    else:
                        Flow_rated=np.random.normal(vol_aver_drawoff_use,vol_desv_drawoff_use,int(draw_offs_use))
    
                    Volume_use=np.array(np.double(np.abs(Flow_rated*time_aver_drawoff_use)))
                
                        
                    for i in range(int(time_steps_day)):
    
                        index=np.where(time_event==i)
    
                        if len(index)>n_max1:
                            index1=index[0,int(n_max1)]
                        else:
                            index1=index
                            
                        if draw_offs_use==1:
                            if index1==0:
                                Volume_use_time=Volume_use
                        else:
                            Volume_use_time=Volume_use[index1]
                        
                        Volume_use_daily_array[(day,i)]=np.sum(Volume_use_time)/time_aver_drawoff_use
    
                Volume_use_resh1=np.reshape(Volume_use_daily_array,(1,int(self.total_days*time_steps_day)))
                Volume_use_array[:,use]=Volume_use_resh1
    
                Volume_use_sum[0,use]=np.sum(Volume_use_array[:,use])
    
            Volume_use_unit[:,units0]=np.sum(Volume_use_array,axis=1)
    
            # if numunits>1:
            Volume_use_arrayb0=np.sum(Volume_use_unit, axis = 1)
            # else:
            #     Volume_use_arrayb0=Volume_use_unit[:,0]
    
        Volume_totalb1 = np.sum(Volume_use_arrayb0)
        
        self.number_units = numunits
        
        self.volume_profile = Volume_use_arrayb0
        self.total_volume = Volume_totalb1
        
        self.Volume_meanb1=((Volume_totalb1)/numunits)/self.total_days
        
    def plot_prof(self):
        plt.plot(self.volume_profile)
        print(f'The number of units is: {self.number_units}')


if __name__ == '__main__':
    V = np.random.randint(50,200,3)
    n_units = 2*np.ones(20, dtype = int)
    ts = 5
    
    
    start = time.time()
    dhw_list = {}
    for i in range(3):
        dhw_list[i] = DHW(V[i], n_units[i], ts)   
    
    total = time.time() - start


#%%
