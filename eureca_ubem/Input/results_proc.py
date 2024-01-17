import os
import glob
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

from eureca_building.config import load_config
load_config("config.json")

snappy_f = glob.glob( os.path.join('belzoni_new_hvac','**.snappy') )

info_bui_gs = gpd.read_file(os.path.join('belzoni_new_hvac','Buildings_summary.geojson'), index = "Index")
res_bui = info_bui_gs[info_bui_gs["End Use"] == "Residential"]["index"]

info_bui = pd.read_csv(os.path.join('belzoni_new_hvac','Buildings_summary.csv'), delimiter = ";", index_col = [0])
info_bui = info_bui.loc[res_bui.values]

results = pd.DataFrame(index=info_bui.index,columns=['DHW_dem_ISO [kWh]', 'DHW_dem_DHW_Calc [kWh]'])
ts_cons_results = pd.DataFrame(columns=pd.MultiIndex.from_product([info_bui.index,['DHW_dem_ISO [W]', 'DHW_dem_DHW_Calc [W]']]),
                                                          index = pd.date_range(start = "02/03/2019 00:00", end = "11/27/2019 23:00", freq = "30min"))


for file in info_bui.index:
    df_ISO = pd.read_parquet(os.path.join('belzoni_new_hvac',f'Results Bd {file}.parquet.snappy'))
    cons_ISO = df_ISO['TZ DHW demand [W]'].sum().sum()/1000
    df_DHW_calc = pd.read_parquet(os.path.join('belzoni_new_hvac_DHW_calc',f'Results Bd {file}.parquet.snappy'))
    cons_DHW_calc = df_DHW_calc['TZ DHW demand [W]'].sum().sum()/1000
    results['DHW_dem_ISO [kWh]'].loc[file] = cons_ISO
    results['DHW_dem_DHW_Calc [kWh]'].loc[file] = cons_DHW_calc

    ts_cons_results[file,'DHW_dem_ISO [W]'] = df_ISO['TZ DHW demand [W]'].values
    ts_cons_results[file,'DHW_dem_DHW_Calc [W]'] = df_DHW_calc['TZ DHW demand [W]'].values

filt = info_bui["Zone net floor area [m2]"] > 5.

demands = ts_cons_results.sum().unstack()
spec_demands = (demands.divide(info_bui["Zone net floor area [m2]"], axis = 0)/1000)

(spec_demands.loc[filt].iloc[400:]).plot(kind = "bar")

#
# plt.close()
#
# fig, [[ax1, ax2], [ax11,ax21]] = plt.subplots(ncols = 2, nrows=2, figsize=(10,10))
#
# ax1.set_axisbelow(True)
# ax2.set_axisbelow(True)
# results['El. Cons [kWh]'].hist(ax = ax1, bins = 50)
# results['El. Cons dw [kWh]'].hist(ax = ax2, bins = 100)
# ax1.set_title('El. Cons [kWh]')
# ax2.set_title('El. Cons dw [kWh]')
# ax1.set_xlabel('El. Cons [kWh]')
# ax2.set_xlabel('El. Cons [kWh]')
# ax1.set_ylabel('Num. [-]')
# ax2.set_ylabel('Num. [-]')
# ax1.set_xlim(0,150000)
#
# results.sort_values('El. Cons dw [kWh]',ascending=True)['El. Cons dw [kWh]'].plot(ax=ax21, grid = True)
# results.sort_values('El. Cons [kWh]',ascending=True)['El. Cons [kWh]'].plot(ax=ax11, grid = True)
# # ax11.set_yscale('log')
# # ax21.set_yscale('log')
# ax11.set_ylabel('El. Cons [kWh]')
# ax21.set_ylabel('El. Cons [kWh]')
# ax11.set_ylabel('El. Cons [kWh]')
# ax21.set_ylabel('El. Cons [kWh]')
# ax11.set_xticks([])
# ax21.set_xticks([])
# plt.tight_layout()
#
# fig.savefig(os.path.join('belzoni_new_hvac','el_loads.svg'))
# plt.show()
#
# s_results = results.sort_values('El. Cons dw [kWh]',ascending=True)
