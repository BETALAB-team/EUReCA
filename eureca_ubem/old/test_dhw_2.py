import logging
import math
import time

import matplotlib.pyplot as plt
import numpy as np

from eureca_building.config import CONFIG
from eureca_building.schedule_properties import domestic_hot_water_prop
from eureca_building.fluids_properties import water_properties
from eureca_building.schedule import Schedule
from eureca_building.exceptions import InvalidScheduleType

def _event_distribution(number_of_daily_events,daily_vector_distibution, number_of_days):
    """DEPRECATED: Not working

    Parameters
    ----------
    n Quanti eventi al giorno
    x
    pdf vettore distribuzione

    Returns
    -------

    """
    y_guess = np.random.rand(number_of_days, int(number_of_daily_events))
    cdf = np.cumsum(daily_vector_distibution)
    x_event = np.interp(y_guess,cdf,np.arange(len(cdf)))
    # get randomly som time events from the cdf and y_guess selected randomly
    return np.round(x_event).astype(int)

sim = 100
start = time.time()

volume_unit = 300  # l/d??
time_step = 5  # min
time_steps_hour = int(60 / time_step)
time_steps_day = int(24 * time_steps_hour)
final_consumption = np.zeros([365*time_steps_hour*24,sim])

for prove in range(sim):

    total_days = domestic_hot_water_prop["total_days"]

    consumed_volume = np.zeros([total_days,time_steps_day,len(domestic_hot_water_prop["DHWcalc_uses"].items())])

    use_number = 0
    for use, values in domestic_hot_water_prop["DHWcalc_uses"].items():

        average_time_drawoff = values["average_time"]
        # l/d split by drawoff type
        daily_vol_total_drawoff_use = volume_unit * values["percentage_to_total_consumption"]
        daily_drawoff_on_event = np.round(daily_vol_total_drawoff_use / values["mean"]).astype(int)
        # Maximum number of on events in a timestep for each drawoff
        # n_max_time_step_on_events = np.ceil(time_step / average_time_drawoff)

        # This does a resampling with a linear interpolation method  to keep 1 as sum of each column
        temporal_dist = np.interp(np.arange(0, 24, 1/time_steps_hour), np.arange(0, 24), values["temporal_distribution"])/time_steps_hour

        # array with the time steps when the use is on [number of daily on events, number of days]
        temporal_dist_year = _event_distribution(daily_drawoff_on_event, temporal_dist, total_days)
        if values["volume_pdf"] == "lognormal":
            m = values["mean"]
            v = values["std"]
            mu = np.log((m ** 2) / np.sqrt(v + m ** 2))
            sigma = np.sqrt(np.log(v / (m ** 2) + 1))
            flow_rate = np.random.lognormal(mu, sigma, [total_days, daily_drawoff_on_event]) # l/min
        elif values["volume_pdf"] == "normal":
            flow_rate = np.random.normal(values["mean"], values["std"], [total_days, daily_drawoff_on_event])
        else:
            raise ValueError(f'DHW calculation. The volume flow rate probabilty distribution function is not allowed. PDF: {values["volume_pdf"]}. Aloowed PDFs: [lognormal, normal]')
        volume_use = np.abs(flow_rate) # l ad accensione
        # To avoid negative consumptions (possible with normal dist)
        volume_use[volume_use < 0] = 0.

        for ts in range(daily_drawoff_on_event):
            consumed_volume[np.arange(total_days), temporal_dist_year[:,ts],  use_number] = volume_use[:,ts]

        # The rescale needed because otherwise the rounding process provoke an underestimation
        consumed_volume[:,:,use_number] = consumed_volume[:,:,use_number] / (consumed_volume[:,:,use_number].sum()/365) * daily_vol_total_drawoff_use

        use_number += 1
    total = consumed_volume.sum(axis =2)

    final_consumption[:,prove] = total.reshape(365*24*time_steps_hour)

import pandas as pd
total_consumed_volume = final_consumption.sum(axis=1)
ts = 2
total_consumed_volume_rs = np.interp(np.arange(0, 365*24, 1/ts), np.arange(0, 365*24, 1/time_steps_hour), total_consumed_volume) * time_steps_hour/ts

total_consumed_volume_rs = pd.Series(total_consumed_volume, index = pd.date_range(start='1/1/2018 00:00', periods=8760*time_steps_hour, freq = f"{time_step}min")).resample(f"{60/ts*60}S").sum()

plt.plot( np.arange(0, 365*24, 1/time_steps_hour),total_consumed_volume)
plt.plot( np.arange(0, 365*24, 1/ts),total_consumed_volume_rs)
plt.show()
print(f"{total_consumed_volume.sum()}")
print(f"{total_consumed_volume_rs.sum()}")
#
# stop = time.time()
# print(f"\nRun {sim} simulations. ",
#       f"\nTotal time: {(stop-start):.2f} s",
#       f"\nTime per sim : {(stop-start)/sim:.2f} s")
#
#
# fig, [ax11,ax12] = plt.subplots(nrows=2)
# ax11.plot(total)
# ax12.plot(total.mean(axis=1))
# plt.show()
#
# fig, [[ax11,ax12],[ax21,ax22]] = plt.subplots(nrows=2,ncols=2)
#
# ax11.plot(consumed_volume[:,:,0].mean(axis=0))
# ax12.plot(consumed_volume[:,:,1].mean(axis=0))
# ax21.plot(consumed_volume[:,:,2].mean(axis=0))
# ax22.plot(consumed_volume[:,:,3].mean(axis=0))
# plt.show()

# Number of days = 365
# nuses = domestic_hot_water_prop["nuses"]
# # Number of drawoffs types
#
# dist = domestic_hot_water_prop["dist"]
# # Creation of arrays
# # Volume_use_daily_array = np.zeros((total_days, int(time_steps_day)))
# # Volume_use_arrayb0 = np.zeros((int(total_days * time_steps_day)))
# # dist_use = np.zeros((24, 12))
# # Volume_use_array = np.zeros((int(time_steps_day * total_days), nuses))
# # Volume_use_sum = np.zeros((1, nuses))
# # Volume_aver_drawoff_final = np.zeros(nuses)
# # Volume_desv_drawoff_final = np.zeros((1, nuses))
# # Volume_use_time = np.zeros(np.int(time_steps_day))
# # Volume_use_unit = np.zeros((int(total_days * time_steps_day), numunits))
#
# for use in range(4):
#
#
#
#     for d in range(total_days):
#         consumed_volume[d,time_event_year[:,d],use] = volume_use[:,d]


#
# ##########################################
# # START FOR CYCLO TO UNIT AND USE
# ###############################
#
# # for units0 in range(numunits):
#
# # for use in range(nuses):
# units0 = 0
# use = 1
#
#
# vol_aver_drawoff_use = vol_mean_drawoff[use]
# vol_desv_drawoff_use = vol_std_dev_drawoff[use]
# time_aver_drawoff_use = time_aver_drawoff[use]
# draw_offs_use = int(draw_offs[use])
# n_max1 = n_max[use]
# dist_use_t = dist_t[:,use]
#
# draw_offs_use = daily_drawoff_on_event[use]
#
# time_event_year = _event_distribution(draw_offs_use, dist_use_t)
# # Gets the daily time steps when some drawoff events occurs
#
# for day in range(total_days):
#
#     time_event = time_event_year[:,day]
#
#     for i in range(int(time_steps_day)):
#
#         index = np.where(time_event == i)
#
#         # index è il timestep del nessimo evento
#         # n_max1 numero di eventi massimi in un timestep
#         # Possono esserci casi in cui il timestep è lo stesso
#
#         if len(index) > n_max1:
#             index1 = index[0, int(n_max1)]
#         else:
#             index1 = index
#
#         if draw_offs_use == 1:
#             # Caso in cui c'è un solo evento al giorno
#             if index1 == 0:
#                 # Parte incomprensibile
#                 Volume_use_time = Volume_use
#         else:
#             Volume_use_time = Volume_use[index1]
#
#         Volume_use_daily_array[(day, i)] = np.sum(Volume_use_time) / time_aver_drawoff_use
#
# Volume_use_resh1 = np.reshape(Volume_use_daily_array, (1, int(total_days * time_steps_day)))
# Volume_use_array[:, use] = Volume_use_resh1
#
# Volume_use_sum[0, use] = np.sum(Volume_use_array[:, use])
#
# Volume_use_unit[:, units0] = np.sum(Volume_use_array, axis=1)
#
# # if numunits>1:
# Volume_use_arrayb0 = np.sum(Volume_use_unit, axis=1)
# # else:
# #     Volume_use_arrayb0=Volume_use_unit[:,0]
#
# Volume_totalb1 = np.sum(Volume_use_arrayb0)
# #
# # number_units = numunits
# #
# # volume_profile = Volume_use_arrayb0
# # total_volume = Volume_totalb1
# #
# # Volume_meanb1 = ((Volume_totalb1) / numunits) / total_days
# return volume_profile