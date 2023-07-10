"""
File with HVAC systems csvs
This is an internal class, where typical systems performances are stored in a disctionary, using csv parsed strings
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

systems_info = {
"CondensingBoiler":
"""
Size [kW]	eta_nom [%]	eta_int [%]	T_gn_w [°C]	f_cor_Pn [-]	P_int [%]	T_gn_Pint [°C]	f_cor_Pint [-]	location [str]
10	98.6	103.9	70	0.2	0.3	35	0.2	internal
30	98.3	104.5	70	0.2	0.3	35	0.2	internal
100	99.1	104.9	70	0.2	0.3	35	0.2	tech_room
300	99.2	105.0	70	0.2	0.3	35	0.2	tech_room
""",

"TraditionalBoiler":
"""
Size [kW]	eta_nom [%]	eta_int [%]	T_gn_w [°C]	f_cor_Pn [-]	P_int [%]	T_gn_Pint [°C]	f_cor_Pint [-]	location [str]
10	92.4	89.6	70	0.04	0.3	50	0.05	internal
30	95.0	91.0	70	0.04	0.3	50	0.05	internal
100	95.2	94.0	70	0.04	0.3	50	0.05	tech_room
300	96.0	98.0	70	0.04	0.3	50	0.05	tech_room
""",

"SplitAirCooler":
"""
Size [kW]	EER_100 [-]	EER_75 [-]	EER_50 [-]	EER_25 [-]
3	2.35	2.68	2.94	2.83
5	2.35	2.68	2.94	2.83
10	2.35	2.68	2.94	2.83
20	2.35	2.68	2.94	2.83
""",

"SplitAirConditioner":
"""
Size [kW]	EER_medio [-]	inverter [str]
3	4.22	yes
5	3.42	yes
10	3.21	yes
15	3.25	yes
""",

"ChillerAirtoWater":
"""
Size [kW]	EER_100 [-]	EER_75 [-]	EER_50 [-]	EER_25 [-]
3	2.35	2.68	2.94	2.83
5	2.35	2.68	2.94	2.83
10	2.35	2.68	2.94	2.83
20	2.35	2.68	2.94	2.83
""",

##################################
############ EN_15316 ############
##################################

"EN_15316_emission_control_heating_efficiency":
"""
Type	Efficiency [-]	Convective fraction [-]
High Temp Radiator	0.820	0.65
Low Temp Radiator	0.882	0.6
Fan coil	0.862	1
Radiant surface	0.857	0.35
""",
"EN_15316_distribution_heating_efficiency":
"""
Type	Efficiency [-]
Centralized	0.92
Single	0.97
""",
"EN_15316_generation_heating_efficiency":
"""
Type	5 kW	10 kW	15 kW	20 kW	30 kW	40 kW	50 kW	100 kW	200 kW	400 kW
Traditional Gas Boiler	0.90048455	0.905	0.907641369	0.90951545	0.912156819	0.9140309	0.91548455	0.92	0.92451545	0.9290309
Oil Boiler	0.90048455	0.905	0.907641369	0.90951545	0.912156819	0.9140309	0.91548455	0.92	0.92451545	0.9290309
Stove	0.7909691	0.8	0.805282738	0.8090309	0.814313638	0.8180618	0.8209691	0.83	0.8390309	0.8480618
Condensing Gas Boiler	0.9869897	0.99	0.991760913	0.9930103	0.994771213	0.9960206	0.9969897	1	1.0030103	1.0060206
A-W Heat Pump	3.2	3.2	3.2	3.2	3.2	3.2	3.2	3.2	3.2	3.2
""",
"EN_15316_generation_heating_auxiliary_electric_load":
"""
Type	5 kW	10 kW	15 kW	20 kW	30 kW	40 kW	50 kW	100 kW	200 kW	400 kW
Traditional Gas Boiler	0.4074	0.4148	0.4222	0.4296	0.4444	0.4592	0.474	0.548	0.696	0.992
Oil Boiler	0.4074	0.4148	0.4222	0.4296	0.4444	0.4592	0.474	0.548	0.696	0.992
Stove	0.15	0.15	0.15	0.15	0.15	0.15	0.15	0.15	0.15	0.15
Condensing Gas Boiler	0.4074	0.4148	0.4222	0.4296	0.4444	0.4592	0.474	0.548	0.696	0.992
A-W Heat Pump	0.	0.	0.	0.	0.	0.	0.	0.	0.	0.
""",
"EN_15316_emission_control_cooling_efficiency":
"""
Type	Efficiency [-]	Convective fraction [-]
Fan coil	0.748	1
Radiant surface	0.777	0.35
Split system	0.842	1
""",
"EN_15316_distribution_cooling_efficiency":
"""
Type	Efficiency [-]
Centralized	0.97
Single	0.97
""",
"EN_15316_generation_cooling_seasonal_performance_factor":
"""
Type	SPF [-]
A-A split	2.27
A-W chiller	2.7
""",
"EN_15316_emission_control_distribution_DHW":
"""
Type	Efficiency [-]
All	0.97
""",
}

# # import abc
# import os
# import io
# import logging
#
# import pandas as pd
# import numpy as np
# global systems_info_dict
# systems_info_dict = {}
# for k,v in systems_info.items():
#     systems_info_dict[k] = pd.read_csv(io.StringIO(v), sep = "\t")
#     if "Type" in systems_info_dict[k].columns:
#         systems_info_dict[k].set_index("Type", drop=True, inplace = True)
# for v, k in systems_info_dict.items():
#     print(k)
#     print("###############")

