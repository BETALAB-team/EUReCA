import pandas as pd
import geopandas as gpd 
import numpy as np 
import json
from eureca_pubem import scenario_process as sc
import itertools

ENV_ORDER = ["none", "shallow", "medium", "deep"]

ENV_ALLOWED = {
    "none": ["none", "shallow", "medium", "deep"],
    "shallow": ["shallow", "medium", "deep"],
    "medium": ["medium", "deep"],
    "deep": ["deep"]
}

HEAT_ALLOWED = {
    "boiler": ["dhn", "hp_le", "hp_me", "hp_he"],
    "dhn": ["dhn", "hp_le", "hp_me", "hp_he"],
    "hp_le": ["dhn", "hp_le", "hp_me", "hp_he"],
    "hp_me": ["dhn", "hp_me", "hp_he"],
    "hp_he": ["hp_he", "dhn"]
}

DHW_ALLOWED = HEAT_ALLOWED.copy()

PV_STEPS = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0]


HP_TYPES = {"hp_le", "hp_me", "hp_he"}

def get_heat_key(value):
    if "boiler" in value:
        return "boiler"
    if "dhn" in value:
        return "dhn"
    if "hp_le" in value:
        return "hp_le"
    if "hp_me" in value:
        return "hp_me"
    if "hp_he" in value:
        return "hp_he"
    raise ValueError(value)

def compare_config_with_base(base_dictionary,
                               baseline_gdf_path,
                               configuration, 
                               ee_measure_path,
                               pv_type_path,
                               hp_catalog_path,
                               grid_pricing_path,
                               dhn_pricing_path,
                               fuels_path,
                               spot_price_path,
                               weatherfile_path,
                               mycity,
                               baseline_scenario,
                               r=0.04,
                               T = 20
                               ):
    interv_dict, building_info = sc.make_dictionary(baseline_geojson = baseline_gdf_path,
                                                   city=mycity,
                                                   baseline_scenario=baseline_scenario,
                                                   weatherfile_path=weatherfile_path,
                                                   intervention_dictionary = configuration)
    
    new_dict = build_dict_gen(building_info,
                       interv_dict,
                       baseline_gdf_path ,
                       configuration, 
                       ee_measure_path,
                       pv_type_path,
                       hp_catalog_path,
                       grid_pricing_path,
                       dhn_pricing_path,
                       fuels_path,
                       spot_price_path

                       )
    
    for idx, building in new_dict.items():
        building["costs"]["financial"] = {}
        baseline_building = base_dictionary[idx]
        costs_ante_hourly = building["costs"]["operational cost"]["hourly_total"]
        costs_ex_hourly = baseline_building["costs"]["operational cost"]["hourly_total"]
        building["costs"]["financial"]["savings_hourly"] = costs_ex_hourly - costs_ante_hourly 
        building["costs"]["financial"]["savings_yearly"] = np.sum(building["costs"]["financial"]["savings_hourly"])
        S = building["costs"]["financial"]["savings_yearly"]
        C = building["costs"]["capital cost"]["total"]
        building["costs"]["financial"]["NPV"] = -C + S * (1 - (1 + r)**(-T)) / r
    return new_dict



def generate_configs_for_building(b0,i=0):
    env_options = ENV_ALLOWED[b0["env"]]
    
    heat_options = HEAT_ALLOWED[get_heat_key(b0["heat"])].copy()
    if b0["heat"] not in heat_options:
        heat_options.append(b0["heat"])
    
    dhw_options = DHW_ALLOWED[get_heat_key(b0["dhw"])].copy()
    if b0["dhw"] not in dhw_options:
        dhw_options.append(b0["dhw"])
    
    pv_min = b0["pv_percentage"]
    pv_options = [p for p in PV_STEPS if p >= pv_min]
    
    configs = []
    
    for env, heat, dhw, pv in itertools.product(env_options, heat_options, dhw_options, pv_options):
        configs.append({
            "env": env,
            "heat": heat,
            "dhw": dhw,
            "fuel": b0["fuel"],
            "pv_percentage": pv
        })
        
        
        # if i ==1 :
        #     print({
        #         "env": env,
        #         "heat": heat,
        #         "dhw": dhw,
        #         "fuel": b0["fuel"],
        #         "pv_percentage": pv
        #     })

    return configs


def optimize_configuration_per_building(config0, 
                                        base_dictionary,
                                        baseline_gdf_path,
                                        ee_measure_path,
                                        pv_type_path,
                                        hp_catalog_path,
                                        grid_pricing_path,
                                        dhn_pricing_path,
                                        fuels_path,
                                        spot_price_path,
                                        weatherfile_path,
                                        mycity,
                                        baseline_scenario,
                                        r=0.04,
                                        T=20):

    candidate_lists = {}
    indices = {}
    best_configs = {}
    best_npvs = {}

    for bid, b0 in config0.items():
        if bid == 48: 
            cands = [b0] + list(generate_configs_for_building(b0, i = 1))
        else:
            cands = [b0] + list(generate_configs_for_building(b0))
        candidate_lists[bid] = cands
        # print(len(cands))
        indices[bid] = 0
        best_configs[bid] = b0
        best_npvs[bid] = 0.0

    active = True

    while active:
        active = False

        current_config = {}

        for bid in config0:
            idx = indices[bid]
            cands = candidate_lists[bid]

            if idx < len(cands):
                current_config[bid] = cands[idx]
                active = True
            else:
                current_config[bid] = best_configs[bid]
                
        

        if not active:
            break

        result = compare_config_with_base(base_dictionary=base_dictionary,
                                          baseline_gdf_path=baseline_gdf_path,
                                          configuration=current_config, 
                                          ee_measure_path=ee_measure_path,
                                          pv_type_path=pv_type_path,
                                          hp_catalog_path=hp_catalog_path,
                                          grid_pricing_path=grid_pricing_path,
                                          dhn_pricing_path=dhn_pricing_path,
                                          fuels_path=fuels_path,
                                          spot_price_path=spot_price_path,
                                          weatherfile_path=weatherfile_path,
                                          mycity=mycity,
                                          baseline_scenario=baseline_scenario,
                                          r=r,
                                          T=T)

        for bid in config0:
            idx = indices[bid]
            cands = candidate_lists[bid]
 
            if idx < len(cands):
                conf = cands[idx]
                npv = result[bid]["costs"]["financial"]["NPV"]


                if npv > best_npvs[bid]:
                    best_npvs[bid] = npv
                    best_configs[bid] = conf
                    
        # print(best_npvs.values())
                    
        for bid in config0:
            if indices[bid] < len(candidate_lists[bid]):
                indices[bid] += 1

    return best_configs, compare_config_with_base(base_dictionary=base_dictionary,
                                      baseline_gdf_path=baseline_gdf_path,
                                      configuration=best_configs, 
                                      ee_measure_path=ee_measure_path,
                                      pv_type_path=pv_type_path,
                                      hp_catalog_path=hp_catalog_path,
                                      grid_pricing_path=grid_pricing_path,
                                      dhn_pricing_path=dhn_pricing_path,
                                      fuels_path=fuels_path,
                                      spot_price_path=spot_price_path,
                                      weatherfile_path=weatherfile_path,
                                      mycity=mycity,
                                      baseline_scenario=baseline_scenario,
                                      r=r,
                                      T=T)
# def optimize_configuration_per_building(config0, 
#                                         base_dictionary,
#                                         baseline_gdf_path,
#                                         ee_measure_path,
#                                         pv_type_path,
#                                         hp_catalog_path,
#                                         grid_pricing_path,
#                                         dhn_pricing_path,
#                                         fuels_path,
#                                         spot_price_path,
#                                         weatherfile_path,
#                                         mycity,
#                                         baseline_scenario,
#                                         r=0.04,
#                                         T = 20
#                                         ):
#     best_configs = {}
    
#     for bid, b0 in config0.items():
#         print(bid)
#         best_npv = 0
#         best_conf = b0
        
#         for conf in generate_configs_for_building(b0):
#             test_config = {bid: conf}

#             result = compare_config_with_base(base_dictionary=base_dictionary,
#                                            baseline_gdf_path=baseline_gdf_path,
#                                            configuration = test_config, 
#                                            ee_measure_path=ee_measure_path,
#                                            pv_type_path=pv_type_path,
#                                            hp_catalog_path=hp_catalog_path,
#                                            grid_pricing_path=grid_pricing_path,
#                                            dhn_pricing_path=dhn_pricing_path,
#                                            fuels_path=fuels_path,
#                                            spot_price_path=spot_price_path,
#                                            weatherfile_path=weatherfile_path,
#                                            mycity=mycity,
#                                            baseline_scenario=baseline_scenario,
#                                            r=r,
#                                            T=T)
            
            
#             npv = result[bid]["costs"]["financial"]["NPV"]
            
#             if npv > best_npv:
#                 best_npv = npv
#                 best_conf = conf
            
#             print(conf, npv, best_npv)
        
#         best_configs[bid] = best_conf
    
#     return best_configs



def extract_config(gdf):
    def get_fuel(val):
        if isinstance(val, str) and "_" in val and "hp" not in val:
            return val.split("_", 1)[1]
        return None

    return {
        row["id"]: {
            "env": row["EEdepth"],
            "heat": row["SHSource"],
            "dhw": row["DHWsource"],
            "fuel": get_fuel(row["SHSource"]) or get_fuel(row["DHWsource"]),
            "pv_percentage": row["PVpercentage"]
        }
        for _, row in gdf.iterrows()
    }

def load_buildings(input_data):

    if isinstance(input_data, gpd.GeoDataFrame):
        gdf = input_data.copy()
    elif isinstance(input_data, str):
        gdf = gpd.read_file(input_data)
    else:
        raise ValueError("invalid_input")

    required_columns = [
        "Name",
        "EEdepth",
        "SHSource",
        "DHWsource",
        "PVType",
        "PVpercentage"
    ]

    missing = [c for c in required_columns if c not in gdf.columns]

    if missing:
        raise ValueError(f"{missing}")

    return gdf

def p95(x):
    return np.percentile(x, 95)

def pick_hp(capacity, catalog):
    for hp in sorted(catalog, key=lambda x: x["max_heating_capacity_kW"]):
        if hp["max_heating_capacity_kW"] >= capacity:
            return hp
    return max(catalog, key=lambda x: x["max_heating_capacity_kW"])

def hp_cost(building, interv, hp_catalog, factor=1.0):
    base = building["base"]
    meta = building["meta"]

    sh_arr = np.array(base["space_heating"])
    dhw_arr = np.array(base["dhw_demand"])

    sh_after = interv.get("sh_source", (None, None))[1]
    dhw_after = interv.get("dhw_source", (None, None))[1]

    sh_final = sh_after if sh_after is not None else meta["SHSource"]
    dhw_final = dhw_after if dhw_after is not None else meta["DHWsource"]

    sh_hp = sh_final in HP_TYPES
    dhw_hp = dhw_final in HP_TYPES

    if not sh_hp and not dhw_hp:
        return 0

    if sh_hp and dhw_hp:
        if sh_final == dhw_final:
            cap = p95(sh_arr + dhw_arr)
            hp = pick_hp(cap, hp_catalog[sh_final.split("_")[1]])
            return hp["cost_eur"] * factor
        else:
            cap_sh = p95(sh_arr)
            cap_dhw = p95(dhw_arr)
            hp1 = pick_hp(cap_sh, hp_catalog[sh_final.split("_")[1]])
            hp2 = pick_hp(cap_dhw, hp_catalog[dhw_final.split("_")[1]])
            return (hp1["cost_eur"] + hp2["cost_eur"]) * factor

    if dhw_hp and not sh_hp:
        if not sh_hp:
            cap = p95(dhw_arr)
            hp = pick_hp(cap, hp_catalog[dhw_final.split("_")[1]])
            return hp["cost_eur"] * factor
        if sh_final != dhw_final:
            cap = p95(dhw_arr)
            hp = pick_hp(cap, hp_catalog[dhw_final.split("_")[1]])
            return hp["cost_eur"] * factor
        cap_sh = p95(sh_arr)
        cap_both = p95(sh_arr + dhw_arr)
        if cap_sh >= cap_both:
            return 0
        hp = pick_hp(cap_both, hp_catalog[sh_final.split("_")[1]])
        return hp["cost_eur"] * factor

    if sh_hp and not dhw_hp:
        if not dhw_hp:
            cap = p95(sh_arr)
            hp = pick_hp(cap, hp_catalog[sh_final.split("_")[1]])
            return hp["cost_eur"] * factor
        if sh_final != dhw_final:
            cap = p95(sh_arr)
            hp = pick_hp(cap, hp_catalog[sh_final.split("_")[1]])
            return hp["cost_eur"] * factor
        cap_dhw = p95(dhw_arr)
        cap_both = p95(sh_arr + dhw_arr)
        if cap_dhw >= cap_both:
            return 0
        hp = pick_hp(cap_both, hp_catalog[sh_final.split("_")[1]])
        return hp["cost_eur"] * factor

    return 0

def dhn_cost(building, interv, connection_cost):
    meta = building["meta"]

    sh_before = interv.get("sh_source", (meta["SHSource"], None))[0]
    dhw_before = interv.get("dhw_source", (meta["DHWsource"], None))[0]

    sh_after = interv.get("sh_source", (None, None))[1]
    dhw_after = interv.get("dhw_source", (None, None))[1]

    after_has_dhn = (sh_after == "dhn") or (dhw_after == "dhn")
    before_has_dhn = (sh_before == "dhn") or (dhw_before == "dhn")

    if after_has_dhn and not before_has_dhn:
        return connection_cost

    return 0

def compute_peak_cost(bought, time_index, price_per_kw):
    df = pd.DataFrame({
        "bought": bought,
        "time": time_index
    })

    df = df.set_index("time")

    monthly_peak = df["bought"].resample("M").max()

    monthly_cost = monthly_peak * price_per_kw

    hourly_cost = monthly_cost.reindex(df.index, method="ffill") / df.resample("M").size().reindex(df.index, method="ffill")

    return hourly_cost.values

def fuel_cost_array(fuel_name, demand_wh, fuels):
    f = fuels[fuel_name]
    var_price = f["price"]["variable"]["value"]
    fixed = f["price"]["fixed"]["value"]

    fixed_per_step = fixed / 8760
    
    return (demand_wh / 1000.0) * var_price + fixed_per_step

def hourly_dhn_cost(operation, p, area):
    demand= operation["dhn_bought"]/1_000_000
    area_fee = area * p["area fee per m2"]
    fixed_heat = p["fixed heat price per MWh"]
    var_heat = p["variable heat price per MWh"]
    admin = p["admin fee yearly"]
    sub_fixed = p["subscription fixed yearly per unit"]
    sub_var = p["subscription variable price per MWh"]
    vat = p["VAT"]

    hourly = demand * (var_heat + sub_var)
    hourly += fixed_heat / 8760
    total = hourly + (area_fee.values[0] + admin + sub_fixed) / 8760
    total *= (1 + vat)
    return total
    
def hourly_grid_cost(operation, prices, pricing):
    df = prices.copy()
    df["time"] = pd.to_datetime(df["time"], utc=True)
    df = df.sort_values("time")

    mask = (df["time"] >= "2025-01-01") & (df["time"] < "2026-01-01")
    df = df.loc[mask]

    time = df["time"]
    price = df["price"].astype(float).values

    bought = np.array(operation["electricity_bought"], dtype=float) / 1000.0
    sold = np.array(operation["electricity_sold"], dtype=float) / 1000.0

    n = min(len(price), len(bought))
    price = price[:n]
    time = time.iloc[:n]
    bought = bought[:n]
    sold = sold[:n]

    buy_cfg = pricing["buy"]
    sell_cfg = pricing["sell"]

    energy_tax = float(buy_cfg["energy tax per kWh"])
    cert = float(buy_cfg["electricity certificate cost per kWh"])
    grid_var = float(buy_cfg["grid local distribution cost monthly per kWh usage"])
    vat = float(buy_cfg["VAT"])

    grid_comp = float(sell_cfg["grid compensation"])
    tax_credit = float(sell_cfg["tax credit"])

    peak_price = float(buy_cfg["grid local distribution cost monthly per kW peak"])
    fixed_monthly = float(buy_cfg["grid local distribution cost monthly fix"])

    buy_price = (price + energy_tax + cert + grid_var) * (1 + vat)
    sell_price = price + grid_comp + tax_credit

    hourly_cost = bought * buy_price - sold * sell_price

    months = time.dt.to_period("M")
    monthly_add = np.zeros(len(hourly_cost))

    for m in months.unique():
        idx = (months == m).values
        peak_kw = np.max(bought[idx])
        hours = np.sum(idx)
        monthly_cost = peak_kw * peak_price + fixed_monthly
        monthly_add[idx] = monthly_cost / hours

    return hourly_cost + monthly_add

def build_dict_gen(building_info,
                   intervention_dict,
                   baseline_gdf_path,
                   configuration, 
                   ee_measure_path,
                   pv_type_path,
                   hp_catalog_path,
                   grid_pricing_path,
                   dhn_pricing_path,
                   fuels_path,
                   spot_price_path
                   ):
    
    
    baseline_gdf = load_buildings(baseline_gdf_path)
    with open(ee_measure_path, "r") as f:
        ee_measure = json.load(f)
    with open(pv_type_path, "r") as f:
        pv_installs = json.load(f)
    with open(hp_catalog_path, "r") as f:
        hp_catalog = json.load(f)
    with open(grid_pricing_path, "r") as f:
        pricing = json.load(f)
    with open(dhn_pricing_path, "r") as f:
        dnn_pricing = json.load(f)
    with open(fuels_path, "r") as f:
        fuels = json.load(f)
    prices = pd.read_csv(spot_price_path) 
    prices["price"] = prices["price"].astype(float)
    Buildings_dict = {}    
    for idx, building in building_info.items():
        EUR_to_SEK = 10.86
        meta={}
        interventions = {}
        operation = {}
        costs = {}
        idx = int(idx)
        if idx in intervention_dict.keys():
            interventions = intervention_dict[idx]
        meta["PV_type"] = baseline_gdf.loc[baseline_gdf["id"] == idx, "PVType"].values[0]
        pv_production = building ["pv_production"]
        electricity_need = building["hp_electricity"] + building["base"]["appliance_electricity"]
        operation["electricity_bought"] = np.maximum(electricity_need - pv_production, 0)
        operation["electricity_sold"] = np.maximum(pv_production - electricity_need , 0)
        operation["dhn_bought"] = building["thermal"]["dhn_demand"]
        isfuel = any(x in building["meta"]["DHWsource"] for x in ["gas", "oil", "bio"]) or any(x in building["meta"]["SHSource"] for x in ["gas", "oil", "bio"])
        if isfuel:
            fuel = next(
                (
                    f for val in [building["meta"]["DHWsource"], building["meta"]["SHSource"]]
                    if isinstance(val, str) and "_" in val and "hp" not in val
                    for f in [val.split("_", 1)[1]]
                    if f in {"gas", "oil", "bio"}
                ),
                None
            )
            
            operation["fuel_bought"] = building["thermal"][f"{fuel}_demand"]
            meta["fuel"]=fuel
            
        area = building["base"]["pv_available_area"]
        floors = baseline_gdf.loc[baseline_gdf["id"] == idx, "Floors"]
        used_area = area*floors

        efficiency_measure_cost = 0
        pv_install_cost = 0
        hp_install_cost = 0 
        dhn_connection_cost = 0
        envelope_type = baseline_gdf.loc[baseline_gdf["id"] == idx, "Envelope"].values[0]
        building_type = envelope_type[:3]
        costs["capital cost"]={}
        costs["operational cost"] = {}
        
        if 'envelope' in intervention_dict[idx].keys():
            efficiency_measure = intervention_dict[idx]['envelope']
            before_measure = efficiency_measure[0]
            if before_measure in ["medium", "deep"]:
                before_measure = before_measure + " " + building_type
            after_measure = efficiency_measure[1]
            if after_measure in ["medium", "deep"]:
                after_measure = after_measure + " " + building_type
            wall_area = building["base"]["opaque_exposed_area"] - building["base"]["pv_available_area"]
            roof_area = building["base"]["pv_available_area"]
            window_area = building ["base"]["glazing_area"]
            wall_cost_fix=ee_measure[after_measure]["wall cost constant"] - ee_measure[before_measure]["wall cost constant"]
            roof_cost_fix=ee_measure[after_measure]["roof cost constant"] - ee_measure[before_measure]["roof cost constant"]
            window_cost_fix=ee_measure[after_measure]["window cost constant"] - ee_measure[before_measure]["window cost constant"]
            wall_cost_persqm=ee_measure[after_measure]["roof cost per square meter"] - ee_measure[before_measure]["roof cost per square meter"]
            roof_cost_persqm=ee_measure[after_measure]["wall cost per square meter"] - ee_measure[before_measure]["wall cost per square meter"]
            window_cost_persqm=ee_measure[after_measure]["window cost per square meter"] - ee_measure[before_measure]["window cost per square meter"]
            efficiency_measure_cost = wall_cost_fix + roof_cost_fix + window_cost_fix \
                                    + wall_cost_persqm * wall_area \
                                    + window_cost_persqm * window_area\
                                    + roof_cost_persqm * roof_area
                                    
        if "PVpercentage" in intervention_dict[idx].keys():
            pv_install = intervention_dict[idx]["PVpercentage"]
            before_measure = pv_install[0]
            after_measure = pv_install[1]
            pv_type = meta["PV_type"]
            pv_install_area = (after_measure - before_measure) * building["base"]["pv_available_area"]/100
            pv_install_fix_cost = pv_installs[pv_type]["cost_fixed"]
            pv_install_persqm_cost = pv_installs[pv_type]["cost_per_m2"]
            pv_install_cost = pv_install_persqm_cost * pv_install_area + pv_install_fix_cost
            
        if any(x in intervention_dict[idx].keys() for x in ["dhw_source","sh_source"]):

            hp_install_cost = hp_cost(building, intervention_dict[idx], hp_catalog, factor=1.0) * EUR_to_SEK 
            dhn_connection_cost =  dhn_cost(building, intervention_dict[idx], connection_cost = 24000)
            

        costs["capital cost"]["efficiency_measure_cost"]  =   efficiency_measure_cost     
        costs["capital cost"]["pv_install_cost"]  =   pv_install_cost      
        costs["capital cost"]["heating_systems"] = hp_install_cost + dhn_connection_cost 
        costs["capital cost"]["total"] = sum(list(costs["capital cost"].values()))
        costs["operational cost"]["electricity"] = np.pad(hourly_grid_cost(operation, prices, pricing["1"]), (8760 - 8759, 0), mode='edge')
        costs["operational cost"]["electricity"] = np.pad(costs["operational cost"]["electricity"], (8760 - len(costs["operational cost"]["electricity"]), 0), mode='edge')
        costs["operational cost"]["district_heating"] = hourly_dhn_cost(operation, dnn_pricing["1"]["buy"], used_area)
        if isfuel:
            costs["operational cost"]["fuel"] = fuel_cost_array(meta["fuel"], operation["fuel_bought"], fuels)*EUR_to_SEK
        costs["operational cost"]["hourly_total"] = np.sum(list(costs["operational cost"].values()), axis=0)
        costs["operational cost"]["yearly_total"] = np.sum(costs["operational cost"]["hourly_total"])
            
        Buildings_dict[idx] = {"meta": meta,
                               "interventions" : interventions,
                               "operation" : operation,
                               "costs": costs}
    
    return Buildings_dict