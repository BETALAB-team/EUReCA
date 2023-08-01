"""
File with HVAC systems csvs
This is an internal class, where typical systems performances are stored in a disctionary, using csv parsed strings
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

from scipy.stats import lognorm
import matplotlib.pyplot as plt
import io
import pandas as pd
import numpy as np

electric_load_italian_distribution = {
"Appliances_penetration":
"""
Regione	Televisore	Frigorifero	Lavatrice	Asciugatrice
ValleDAosta	0.990	0.944	0.915	0.903
Piemonte	0.993	0.918	0.996	0.995
Liguria	0.997	0.992	0.906	0.958
Lombardia	0.976	0.989	0.961	0.995
Veneto	0.932	1.000	0.942	0.919
TrentinoAltoAdige	0.927	0.981	0.956	0.921
FriuliVeneziaGiulia	0.999	0.906	0.909	0.904
EmiliaRomagna	0.990	0.970	0.971	0.943
Umbria	0.950	0.915	0.988	0.955
Toscana	0.995	0.989	0.980	0.905
Marche	0.934	0.909	0.968	0.910
Abruzzo	0.927	0.971	0.984	0.904
Lazio	0.995	0.955	0.942	0.909
Campania	0.980	0.932	0.953	0.974
Basilicata	0.936	0.925	0.966	0.912
Molise	0.960	0.903	0.970	0.973
Puglia	0.959	0.976	0.979	0.965
Calabria	0.914	0.927	0.959	0.916
Sicilia	0.917	0.957	0.921	0.938
Sardegna	0.953	0.952	0.998	0.945
""",
"National_distributions": {
    "Televisore":{
            "Name":"lognorm",
            "s" : 0.9,
        },
    "Frigorifero":{
            "Name":"lognorm",
            "s" : 0.9,
        },
    "Lavatrice":{
            "Name":"lognorm",
            "s" : 0.9,
        },
    "Asciugatrice":{
            "Name":"lognorm",
            "s" : 0.9,
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
        electric_load_italian_dict[f"PDF {k}"] = lognorm(el["s"])
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
    values = pd.DataFrame(np.random.random(size = (number,4)), columns = electric_load_italian_dict["Appliances_penetration"].columns)
    appl_presence = electric_load_italian_dict["Appliances_penetration"].loc[region] > values
    # appl = appl_presence[appl_presence == True]
    # appl = pd.Series(np.random.random(size = len(appl)), index = appl.index)
    loads = pd.DataFrame()
    for ap in appl_presence.columns:
        loads[ap] = electric_load_italian_dict[f"PDF {ap}"].rvs(size = number) * appl_presence[ap]

    loads["Tot"] = loads.sum(axis =1)

    return loads

# import time
# st = time.time()
# loads = get_italian_random_el_loads(10000,"Veneto")
#
# print(f"{time.time()-st:.2f} s")
#
# df = loads["Televisore"]
# df.hist(bins= 25)
# ax = plt.gca()
# ax_ = ax.twinx()
# x = np.linspace(0,
#                     20, 100)
#
# ax_.plot(x, electric_load_italian_dict[f"PDF {k}"].pdf(x), lw=2, alpha=0.6, label=f'{k}')
# plt.show()