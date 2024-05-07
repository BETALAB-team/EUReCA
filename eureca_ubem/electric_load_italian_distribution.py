"""
File with HVAC systems csvs
This is an internal class, where typical systems performances are stored in a disctionary, using csv parsed strings
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

from scipy.stats import lognorm, triang, weibull_min, gamma, norm, expon, uniform
#import matplotlib.pyplot as plt
import io
import pandas as pd
import numpy as np

electric_load_italian_distribution = {
"Appliances_penetration":
"""
Regione	Refrigerator	Freezer	Washing machine	Clothes dryer	Dish washer	Electric cooking	Electric oven	Television	Monitor	Light	Small appliances	Cooling split
Piemonte	0.996	0.248	0.973	0.016	0.442	0.032	0.736	0.954	0.504	0.985	0.984	0.063
ValleDAosta	0.989	0.395	0.972	0.035	0.427	0.053	0.816	0.964	0.515	0.982	0.994	0.004
Lombardia	0.998	0.238	0.963	0.055	0.490	0.025	0.764	0.976	0.566	0.988	0.993	0.155
Veneto	0.998	0.357	0.973	0.045	0.505	0.023	0.661	0.954	0.530	0.972	0.996	0.237
FriuliVeneziaGiulia	0.991	0.310	0.965	0.066	0.432	0.025	0.789	0.951	0.550	0.980	0.992	0.145
Liguria	0.993	0.243	0.949	0.006	0.379	0.026	0.670	0.959	0.494	0.978	0.989	0.084
EmiliaRomagna	0.992	0.302	0.965	0.072	0.427	0.021	0.706	0.968	0.534	0.980	0.992	0.242
Toscana	0.995	0.290	0.951	0.064	0.548	0.005	0.835	0.964	0.548	0.977	0.977	0.112
Umbria	0.989	0.390	0.960	0.014	0.425	0.017	0.786	0.952	0.508	0.973	0.990	0.056
Marche	0.995	0.336	0.956	0.042	0.497	0.014	0.751	0.961	0.605	0.985	0.997	0.121
Lazio	0.993	0.204	0.950	0.010	0.351	0.029	0.621	0.960	0.510	0.981	0.982	0.149
Abruzzo	0.994	0.387	0.964	0.016	0.425	0.029	0.634	0.980	0.461	0.985	0.970	0.066
Molise	0.998	0.258	0.981	0.030	0.359	0.017	0.719	0.973	0.430	0.984	0.987	0.051
Campania	0.999	0.203	0.981	0.008	0.220	0.037	0.697	0.975	0.546	0.982	0.988	0.144
Puglia	0.998	0.125	0.955	0.012	0.308	0.029	0.736	0.947	0.526	0.981	0.980	0.187
Basilicata	0.994	0.222	0.974	0.031	0.294	0.061	0.780	0.971	0.487	0.987	0.980	0.111
Calabria	0.998	0.260	0.954	0.023	0.294	0.025	0.725	0.971	0.481	0.976	0.966	0.150
Sicilia	0.997	0.141	0.953	0.004	0.179	0.016	0.697	0.952	0.490	0.981	0.982	0.225
Sardegna	0.998	0.404	0.968	0.027	0.252	0.022	0.631	0.967	0.530	0.989	0.984	0.228
TrentinoAltoAdige	0.999	0.371	0.983	0.063	0.533	0.288	0.853	0.958	0.566	0.996	0.993	0.025
""",
"National_distributions": {
    "Refrigerator":{
            "Name":"lognorm",
            "s" : 0.3471,
            "scale" : 350.75921818799463,
            "loc" : 0.,
        },
    "Freezer":{
            "Name":"expon",
            "scale" : 362.5329,
            "loc" : 0.,
        },
    "Washing machine":{
            "Name":"norm",
            "scale" : 131.7528,
            "loc" : 263.7713,
        },
    "Clothes dryer":{
            "Name":"expon",
            "scale" : 299.2194,
            "loc" : 0.,
        },
    "Dish washer": {
        "Name": "expon",
        "scale": 325.6673,
        "loc": 0.,
    },
    "Electric cooking": {
        "Name": "triang",
        "c": 1.,
        "scale": 1187.7350999999999,
        "loc": 69.8649,
    },
    "Electric oven": {
        "Name": "expon",
        "scale": 136.2937,
        "loc": 0.,
    },
    "Television": {
        "Name": "weibull_min",
        "c": 1.865,
        "scale": 192.6315,
    },
    "Monitor": {
        "Name": "weibull_min",
        "c": 1.8383,
        "scale": 151.6776,
    },
    "Light": {
        "Name": "lognorm",
        "s": 0.6698,
        "scale": 327.0784335215426,
        "loc": 0.,
    },
    "Small appliances": {
        "Name": "gamma",
        "a": 6.3351,
        "scale": 70.2625,
    },
    "Cooling split": {
        "Name": "lognorm",
        "s": 0.83,
        "scale": 326.6861748036186,
        "loc": 0.,
    },
}
}

global electric_load_italian_dict
electric_load_italian_dict = {}
v = electric_load_italian_distribution["Appliances_penetration"]
electric_load_italian_dict["Appliances_penetration"] = pd.read_csv(io.StringIO(v), sep = "\t")
electric_load_italian_dict["Appliances_penetration"].set_index("Regione", drop=True, inplace = True)

for k, el in electric_load_italian_distribution["National_distributions"].items():
    if el["Name"] == "lognorm":
        electric_load_italian_dict[f"PDF {k}"] = lognorm(el["s"], scale = el["scale"], loc = el["loc"])
    elif el["Name"] == "expon":
        electric_load_italian_dict[f"PDF {k}"] = expon(scale = el["scale"], loc = el["loc"])
    elif el["Name"] == "norm":
        electric_load_italian_dict[f"PDF {k}"] = norm(scale = el["scale"], loc = el["loc"])
    elif el["Name"] == "triang":
        electric_load_italian_dict[f"PDF {k}"] = triang(el["c"], scale = el["scale"], loc = el["loc"])
    elif el["Name"] == "weibull_min":
        electric_load_italian_dict[f"PDF {k}"] = weibull_min(el["c"], scale = el["scale"])
    elif el["Name"] == "gamma":
        electric_load_italian_dict[f"PDF {k}"] = gamma(el["a"], scale = el["scale"])
    elif el["Name"] == "uniform":
        raise NotImplementedError
    else:
        electric_load_italian_dict[f"PDF {k}"] = None

def get_italian_random_el_loads(number, region):
    """This function creates a df of stochastic annual appliances electric consumption from ISTAT dataset.
    The calculation is done for each type of appliance and a input number of dwellings

    Parameters
    ----------
    number : int
        number of dwellings
    region : str
        Italian region from which take the data. Available:
        ValleDAosta, Piemonte, Liguria, Lombardia, Veneto, TrentinoAltoAdige,
        FriuliVeneziaGiulia, EmiliaRomagna, Umbria, Toscana, Marche, Abruzzo,
        Lazio, Campania, Basilicata, Molise, Puglia, Calabria, Sicilia, Sardegna

    Returns
    -------
    pandas.DataFrame
        Dataframe with row number of units to get, columns appliances types

    Raises
    ------
    ValueError
        Error in case reagion is not an allowed value
    """
    if region not in ["ValleDAosta", "Piemonte", "Liguria", "Lombardia", "Veneto", "TrentinoAltoAdige",
        "FriuliVeneziaGiulia", "EmiliaRomagna", "Umbria", "Toscana", "Marche", "Abruzzo",
        "Lazio", "Campania", "Basilicata", "Molise", "Puglia", "Calabria", "Sicilia", "Sardegna"]:
        raise ValueError(f"""
Region not valid: {region}. Allowed regions: ValleDAosta, Piemonte, Liguria, Lombardia, Veneto, TrentinoAltoAdige,
FriuliVeneziaGiulia, EmiliaRomagna, Umbria, Toscana, Marche, Abruzzo,
Lazio, Campania, Basilicata, Molise, Puglia, Calabria, Sicilia, Sardegna
""")
    values = pd.DataFrame(np.random.random(size = (number,12)), columns = electric_load_italian_dict["Appliances_penetration"].columns)
    appl_presence = electric_load_italian_dict["Appliances_penetration"].loc[region] > values
    # appl = appl_presence[appl_presence == True]
    # appl = pd.Series(np.random.random(size = len(appl)), index = appl.index)
    loads = pd.DataFrame()
    for ap in appl_presence.columns:
        loads[ap] = electric_load_italian_dict[f"PDF {ap}"].rvs(size = number) * appl_presence[ap]

    loads["Tot"] = loads.sum(axis =1)

    return loads #, loads[loads!=0]


# import time
# import os
# st = time.time()
# loads, loads_wo_zeros = get_italian_random_el_loads(1000000,"Veneto")
#
# for l in loads:
#     print(f"{time.time()-st:.2f} s")
#
#     fig, ax = plt.subplots(figsize = (15,15))
#     df = loads_wo_zeros[l]
#     df.hist(ax = ax,bins= 100, density=True)
#     ax_ = ax.twinx()
#     x = np.linspace(0,1000, 100)
#
#     max_lim = 0.001
#     if l!= 'Tot':
#         ax_.plot(x, electric_load_italian_dict[f"PDF {l}"].pdf(x), lw=2, alpha=0.6, label=f'{k}', color = 'r')
#         max_lim = np.max(electric_load_italian_dict[f"PDF {l}"].pdf(x))
#     ax_.set_ylim(0,max_lim)
#     ax.set_ylim(0, max_lim)
#     fig.suptitle(f"{l}")
#     plt.tight_layout()
#     fig.savefig(os.path.join("el_load_dist",f"el_load_{l}.png"))
#     #plt.show()
#     plt.close()