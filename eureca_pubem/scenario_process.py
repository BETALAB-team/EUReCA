import numpy as np 
import geopandas as gpd
import json
import pandas as pd
from eureca_pubem.scenario.orchestrator import  create_scenario, apply_dhn_changes
from eureca_pubem import scenario_gui as app

def load_hp_catalogue(path):
    with open(path, "r") as f:
        return json.load(f)
    

def choose_heat_pump(hp_demand_array, efficiency_class, hp_catalogue):
    """
    hp_demand_array : thermal kWh array (hourly)
    efficiency_class : 'le' | 'me' | 'he'
    hp_catalogue : loaded json dict
    """

    peak_load = np.max(hp_demand_array)  # kW if hourly kWh

    candidates = sorted(
        hp_catalogue[efficiency_class],
        key=lambda x: x["max_heating_capacity_kW"]
    )

    for hp in candidates:
        if hp["max_heating_capacity_kW"] >= peak_load:
            return hp

    raise ValueError(
        f"No heat pump in class {efficiency_class} "
        f"can cover peak load {peak_load:.2f} kW"
    )

def load_geojson_with_envelope_prefix(path: str, depth: int):
    if depth < 0:
        raise ValueError("n must be 0 or positive")

    gdf = gpd.read_file(path)

    if "Envelope" not in gdf.columns:
        raise KeyError("Column 'Envelope' not found in GeoJSON")

    if depth > 0:
        gdf["Envelope"] = gdf["Envelope"].astype(str) + "+" * depth 

    return gdf


def create_gdf_dictionary (city_object):

    rows = []            
    values_dict = {}    

    for ident, bd in city_object.buildings_objects.items():

        space_heating = sum(
            np.maximum(obj.sensible_load_hourly + obj.latent_load, 0)
            for obj in bd._thermal_zones_list
        )

        space_cooling = -sum(
            np.minimum(obj.sensible_load_hourly + obj.latent_load, 0)
            for obj in bd._thermal_zones_list
        )

        dhw_demand = sum(
            obj.domestic_hot_water_demand
            for obj in bd._thermal_zones_list
        )

        appliance_electricity = sum(
            obj.electric_load
            for obj in bd._thermal_zones_list
        )

        opaque_exposed_area = sum(
            obj.ext_roof_area + obj.ext_wall_opaque_area
            for obj in bd._thermal_zones_list
        )

        glazing_area = sum(
            obj.ext_wall_glazed_area
            for obj in bd._thermal_zones_list
        )

        pv_available_area = sum(
            obj.ext_roof_area
            for obj in bd._thermal_zones_list
        )

        geometry = city_object.buildings_info[ident]["geometry"]
        centroid = geometry.centroid

        values_dict[ident] = {
            "space_heating": space_heating,
            "space_cooling": space_cooling,
            "dhw_demand": dhw_demand,
            "appliance_electricity": appliance_electricity,
            "opaque_exposed_area": opaque_exposed_area,
            "glazing_area": glazing_area,
            "pv_available_area": pv_available_area,
        }

        rows.append({
            "building_id": ident,
            "geometry": centroid
        })

    gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs=city_object.cityjson.crs)

    return gdf, values_dict

def build_thermal_tech_arrays(geo, mycity):
    """
    Creates per-building thermal demand arrays:
        - dhn_demand
        - hp_demand
        - split_demand
        - {fuel}_demand (dynamic keys)

    All values remain thermal kWh (same as original arrays).
    """

    gdf = gpd.read_file(geo) if isinstance(geo, str) else geo
    results = {}

    for _, row in gdf.iterrows():
        bid = str(row["id"])
        level = row["EEdepth"]

        data = mycity[level][bid]

        sh = data["space_heating"]
        dhw = data["dhw_demand"]
        sc = data["space_cooling"]

        # initialize
        tech_arrays = {
            "dhn_demand": np.zeros_like(sh),
            "hp_demand": np.zeros_like(sh),
            "split_demand": np.zeros_like(sh),
        }


        sh_source = row["SHSource"]

        if sh_source == "dhn":
            tech_arrays["dhn_demand"] += sh
        elif sh_source.startswith("hp"):
            tech_arrays["hp_demand"] += sh
        elif sh_source.startswith("boiler_"):
            fuel = sh_source.split("_", 1)[1]
            tech_arrays.setdefault(f"{fuel}_demand", np.zeros_like(sh))
            tech_arrays[f"{fuel}_demand"] += sh


        dhw_source = row["DHWsource"]

        if dhw_source == "dhn":
            tech_arrays["dhn_demand"] += dhw
        elif dhw_source.startswith("hp"):
            tech_arrays["hp_demand"] += dhw
        elif dhw_source.startswith("boiler_"):
            fuel = dhw_source.split("_", 1)[1]
            tech_arrays.setdefault(f"{fuel}_demand", np.zeros_like(sh))
            tech_arrays[f"{fuel}_demand"] += dhw


        sc_source = row["SCsource"]

        if sc_source == "split":
            tech_arrays["split_demand"] += sc
        elif sc_source.startswith("hp"):
            tech_arrays["hp_demand"] += sc

        results[bid] = tech_arrays
    return results

def demand_aggregation(geo, mycity):
    with open("C:/Works/EUReCA/EUReCA/eureca_pubem/hp_config.json", "r") as f:
        hp_catalogue = json.load(f)
    gdf = gpd.read_file(geo) if isinstance(geo, str) else geo
    epw = "C:/Works/EUReCA/EUReCA/eureca_ubem/Input/SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw"
    buildings = initialize_buildings(gdf, mycity)
    T_out = load_epw_drybulb(epw)
    epw_df = pd.read_csv(epw, skiprows=8, header=None)
    
    GHI = epw_df.iloc[:, 13].astype(float).to_numpy()  
    T_air = epw_df.iloc[:, 6].astype(float).to_numpy() 
    calculate_thermal_demands(buildings)
    choose_heat_pumps(buildings, hp_catalogue)
    calculate_pv_production(gdf, buildings, GHI, T_air)
    calculate_hp_electricity(buildings, T_out)
    # calculate_boiler_consumption(buildings)
    # calculate_split_consumption(buildings)
    # calculate_total_electricity(buildings)

    return buildings

def initialize_buildings(geo, mycity):
    """
    Initializes building dictionary.

    Returns:
        dict[building_id] = {
            meta: {...},
            base: {...raw arrays from mycity...}
        }
    """

    gdf = gpd.read_file(geo) if isinstance(geo, str) else geo
    buildings = {}

    for _, row in gdf.iterrows():
        bid = str(row["id"])
        level = row["EEdepth"]

        buildings[bid] = {
            "meta": {
                "level": level,
                "SHSource": row.get("SHSource"),
                "DHWsource": row.get("DHWsource"),
                "SCsource": row.get("SCsource"),
            },
            "base": mycity[level][bid],
        }

    return buildings

import numpy as np


def calculate_thermal_demands(buildings):
    """
    Updates buildings dict in place.

    Adds:
        buildings[bid]["thermal"] = {
            dhn_demand,
            hp_demand,
            split_demand,
            {fuel}_demand
        }
    """

    for bid, b in buildings.items():

        base = b["base"]
        meta = b["meta"]

        sh = base["space_heating"]
        dhw = base["dhw_demand"]
        sc = base["space_cooling"]

        tech = {
            "dhn_demand": np.zeros_like(sh),
            "hp_demand": np.zeros_like(sh),
            "split_demand": np.zeros_like(sh),
        }

        def assign(source, arr):
            if not isinstance(source, str):
                return

            if source == "dhn":
                tech["dhn_demand"] += arr

            elif source.startswith("hp"):
                tech["hp_demand"] += arr

            elif source.startswith("boiler_"):
                fuel = source.split("_", 1)[1]
                key = f"{fuel}_demand"
                tech.setdefault(key, np.zeros_like(sh))
                tech[key] += arr

            elif source == "split":
                tech["split_demand"] += arr

        assign(meta.get("SHSource"), sh)
        assign(meta.get("DHWsource"), dhw)
        assign(meta.get("SCsource"), sc)

        b["thermal"] = tech
        
import numpy as np


def choose_heat_pumps(buildings, hp_catalogue, oversize_factor=1.0):
    """
    Updates buildings in place.

    Adds:
        buildings[bid]["hp"] = {
            model,
            max_heating_capacity_kW,
            SCOP,
            cost_eur,
            design_peak_kW
        }
        OR None
    """
    for bid, b in buildings.items():

        thermal = b.get("thermal", {})
        hp_array = thermal.get("hp_demand")/1000
        if hp_array is None or np.max(hp_array) <= 0:
            b["hp"] = None
            continue

        efficiency_class = None
        for field in ["SHSource", "DHWsource", "SCsource"]:
            source = b["meta"].get(field)
            if isinstance(source, str) and source.startswith("hp_"):
                efficiency_class = source.split("_", 1)[1]

                break

        if efficiency_class is None:
            b["hp"] = None
            continue

        peak = np.max(hp_array) * oversize_factor
        candidates = sorted(
            hp_catalogue[efficiency_class],
            key=lambda x: x["max_heating_capacity_kW"]
        )

        chosen = None
        for hp in candidates:
            if hp["max_heating_capacity_kW"] >= peak:
                chosen = {
                    **hp,
                    "design_peak_kW": peak
                }
                break

        if chosen is None:
            raise ValueError(
                f"No HP in class {efficiency_class} "
                f"covers {peak:.2f} kW (building {bid})"
            )

        b["hp"] = chosen
        
        
import numpy as np
def load_epw_drybulb(epw_path):
    df = pd.read_csv(
        epw_path,
        skiprows=8,
        header=None
    )

    # EPW drybulb column index = 6
    # (year,month,day,hour,min,...,drybulb)
    drybulb = df.iloc[:, 6].astype(float).to_numpy()

    return drybulb

def calculate_hp_electricity(
    buildings,
    epw,
    supply_temp_C=45.0,
    delta_T=5.0,
    cop_max=6.0,
    cop_min=1.2,
):
    """
    Carnot-based COP model with calibration to rated SCOP at 7°C.

    Updates buildings in place:
        - hp_cop
        - hp_electricity
    """
    if isinstance(epw, dict):
        T_out_C = np.asarray(epw["drybulb"], dtype=float)
    else:
        T_out_C = np.asarray(epw, dtype=float)
    T_out_K = T_out_C + 273.15
    T_hot_K = supply_temp_C + 273.15

    T_ref_C = 7.0
    T_ref_K = T_ref_C + 273.15

    for bid, b in buildings.items():

        hp = b.get("hp")
        if not hp:
            b["hp_cop"] = np.zeros_like(T_out_K)
            b["hp_electricity"] = np.zeros_like(T_out_K)
            continue

        scop = hp["SCOP"]
        thermal = b["thermal"]["hp_demand"]

        cop_carnot_ref = T_hot_K / (T_hot_K - T_ref_K + delta_T)

        eta_c = scop / cop_carnot_ref

        cop_carnot = T_hot_K / (T_hot_K - T_out_K + delta_T)

        COP = eta_c * cop_carnot

        COP = np.clip(COP, cop_min, cop_max)

        electricity = np.divide(
            thermal,
            COP,
            out=np.zeros_like(thermal),
            where=COP > 0
        )

        b["hp_cop"] = COP
        b["hp_electricity"] = electricity
        
        
def electricity_df(data):
    return pd.DataFrame({
        k: v["base"]["appliance_electricity"] +
           (v["hp_electricity"] if v["hp_electricity"] is not None else 0)
        for k, v in data.items()
    })

def dhn_df(data):
    df = pd.DataFrame({
        k: v["thermal"]["dhn_demand"] 
        for k, v in data.items()
    })
    df=df/1000
    return df

def pv_df(data):
    return pd.DataFrame({
        k: v["pv_production"] if v["pv_production"] is not None else 0
        for k, v in data.items()
    })

def dataframes(data):
    
    return{"electricity_consumption": electricity_df(data),
           "pv_generation": pv_df(data),
           "dhn_demand":dhn_df(data)}

def calculate_pv_production(gdf, buildings, GHI, T_air):

    with open("C:/Works/EUReCA/EUReCA/eureca_pubem/pv_config.json", "r") as f:
        pv_catalogue = json.load(f)

    for _, row in gdf.iterrows():

        bid = str(row["id"])

        pv_type = row.get("PVType")
        pv_percentage = row.get("PVpercentage", 0)

        if pv_percentage == 0 or pv_type not in pv_catalogue:
            buildings[bid]["pv_production"] = np.zeros_like(GHI)
            continue

        pv_data = pv_catalogue[pv_type]

        module_area = pv_data["module_area"]
        p_stc = pv_data["max_power"]
        noct = pv_data["noct"]
        beta = pv_data["temperature_coefficient"]

        base = buildings[bid]["base"]
        available_area = base["pv_available_area"]

        installed_area = available_area * pv_percentage / 100
        n_modules = installed_area / module_area

        irradiance = GHI

        T_cell = T_air + irradiance * (noct - 20) / 800
        eta_temp = 1 + beta * (T_cell - 25)

        P_module = p_stc * (irradiance / 1000) * eta_temp
        P_total = P_module * n_modules

        buildings[bid]["pv_production"] = P_total
        
        
def create_baseline(input_gdf, 
                    input_city, 
                    baseline_gdf, 
                    distributions_geojson,
                    weather_path, 
                    street_path):
    input_gdf["id"] = input_gdf["building_id"].astype(int)
    def create_network_gdfs(input_gdf, geojson_path):

        geo = gpd.read_file(geojson_path)
        geo = geo.to_crs(input_gdf.crs)
        electrical_geo = geo[geo["type"] == "electrical"].copy()
        dhn_geo = geo[geo["type"] == "dhn"].copy()



        def next_ids(n, used_ids):
            ids = []
            i = 0
            while len(ids) < n:
                if i not in used_ids:
                    ids.append(i)
                    used_ids.add(i)
                i += 1
            return ids



        used_ids = set(input_gdf["id"])

        electrical_geo["id"] = next_ids(len(electrical_geo), used_ids)

        electrical_input = input_gdf[["id", "geometry"]].copy()
        electrical_input["type"] = -1

        electrical_geo = electrical_geo[["id", "geometry"]]
        electrical_geo["type"] = +1

        electrical_gdf = gpd.GeoDataFrame(
            pd.concat([electrical_input, electrical_geo], ignore_index=True),
            geometry="geometry",
            crs=input_gdf.crs
        )


        used_ids = set(input_gdf["id"]) | set(electrical_gdf["id"])

        dhn_geo["id"] = next_ids(len(dhn_geo), used_ids)

        dhn_input = input_gdf[["id", "geometry"]].copy()
        dhn_input["type"] = -1
        dhn_input["capacity"] = None

        dhn_geo = dhn_geo[["id", "geometry", "capacity"]]
        dhn_geo["type"] = +1

        dhn_gdf = gpd.GeoDataFrame(
            pd.concat([dhn_input, dhn_geo], ignore_index=True),
            geometry="geometry",
            crs=input_gdf.crs
        )

        dhn_gdf = dhn_gdf[["id", "type", "capacity", "geometry"]]
        electrical_gdf = electrical_gdf[["id", "type", "geometry"]]

        return electrical_gdf, dhn_gdf

    elec_gdf, dhn_gdf = create_network_gdfs(input_gdf=input_gdf, geojson_path=distributions_geojson)
    a = demand_aggregation(baseline_gdf, input_city)
    b = dataframes(a)
    streets = gpd.read_file(street_path)
    streets =streets.to_crs(input_gdf.crs)

    DESIGN = create_scenario (  dhn_demand_path = b["dhn_demand"],
                                dhn_nodes_path = dhn_gdf,
                                elec_consumption_path =  b["electricity_consumption"],
                                elec_production_path =  b["pv_generation"],
                                elec_nodes_path = elec_gdf,
                                climate_file_path = weather_path, 
                                streets_path = streets)
    
    BASELINE = create_scenario ( 
            dhn_demand_path = b["dhn_demand"],
            climate_file_path = weather_path, 
            mode = "BASELINE",
            base = DESIGN
            )
    
    return BASELINE

def analyze_intervention(baseline_geojson, city, baseline_scenario, weatherfile_path, mode="dictionary", intervention_dictionary = None):
    if mode == "dictionary" :
        # interv_dict = intervention_dictionary
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
        def build_intervention_dict(baseline, gdf):

            interventions = {}
        
            for idx in gdf.index:
        
                base = baseline.loc[idx]
                new  = gdf.loc[idx]
        
                changes = {}
        
                if base["EEdepth"] != new["EEdepth"]:
                    changes["envelope"] = (base["EEdepth"], new["EEdepth"])
        
                if base["SHSource"] != new["SHSource"]:
                    changes["sh_source"] = (base["SHSource"], new["SHSource"])
        
                if base["DHWsource"] != new["DHWsource"]:
                    changes["dhw_source"] = (base["DHWsource"], new["DHWsource"])
        
                if base["PVType"] != new["PVType"]:
                    changes["PVtype"] = (base["PVType"], new["PVType"])
        
                if float(base["PVpercentage"]) != float(new["PVpercentage"]):
                    changes["PVpercentage"] = (
                        float(base["PVpercentage"]),
                        float(new["PVpercentage"])
                    )
        
                interventions[idx] = changes
        
            return interventions
        def building_editor_from_dict(geojson_path, config):
        
            gdf = load_buildings(baseline_geojson)
            baseline = gdf.copy()
        
            _validate_config(gdf, config)
        
            gdf = _apply_config(gdf, config)
        
            interventions = build_intervention_dict(baseline, gdf)
        
            return {"gdf": gdf}, interventions
        
        
        def _validate_config(gdf, config):
            pass
        
        
        def _apply_config(gdf, config):
        
            for idx, row in config.items():
                if idx not in gdf["id"].values:
                    raise ValueError(f"{idx}")
        
                env  = row.get("env")
                heat = row.get("heat")
                dhw  = row.get("dhw")
                fuel = row.get("fuel")
                pv_type = row.get("pv_type")
                pv_perc = row.get("pv_percentage")
        
                if env is not None:
                    gdf.loc[gdf["id"] == idx, "EEdepth"] = env
        
                if heat is not None:
                    if heat == "boiler":
                        if not fuel:
                            raise ValueError(f"{idx}")
                        gdf.loc[gdf["id"] ==idx, "SHSource"] = f"boiler_{fuel}"
                    else:
                        gdf.loc[gdf["id"] ==idx, "SHSource"] = heat
        
                if dhw is not None:
                    if dhw == "boiler":
                        if not fuel:
                            raise ValueError(f"{idx}")
                        gdf.loc[gdf["id"] ==idx, "DHWsource"] = f"boiler_{fuel}"
                    else:
                        gdf.loc[gdf["id"] ==idx, "DHWsource"] = dhw
        
                if pv_type is not None:
                    gdf.loc[gdf["id"] ==idx, "PVType"] = pv_type
        
                if pv_perc is not None:
                    gdf.loc[gdf["id"] ==idx, "PVpercentage"] = max(
                        float(gdf.loc[gdf["id"] ==idx, "PVpercentage"]),
                        float(pv_perc)
                    )
            return gdf
        
        q , interv_dict = building_editor_from_dict(baseline_geojson, intervention_dictionary)
        a = demand_aggregation(q["gdf"], city)
        b = dataframes (a)
        #something to make the intervention dataframe
    if mode == "app":
        q, interv_dict = app.building_editor(baseline_geojson)
        a = demand_aggregation(q["gdf"], city)
        b = dataframes(a) 
    NEW = create_scenario ( 
            dhn_demand_path =  b["dhn_demand"],
            elec_consumption_path = b["electricity_consumption"],
            elec_production_path =  b["pv_generation"],
            climate_file_path = weatherfile_path, 
            mode = "RETROFIT",
            base = baseline_scenario
            )
    
    return NEW, a, interv_dict


def make_dictionary(baseline_geojson, city, baseline_scenario, weatherfile_path, mode="dictionary", intervention_dictionary = None):
    if mode == "dictionary" :
        # interv_dict = intervention_dictionary
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
        def build_intervention_dict(baseline, gdf):

            interventions = {}
            gdf = gdf.set_index("id")
            baseline = baseline.set_index("id")
        
            for idx in gdf.index:
        
                base = baseline.loc[idx]
                new  = gdf.loc[idx]
        
                changes = {}
        
                if base["EEdepth"] != new["EEdepth"]:
                    changes["envelope"] = (base["EEdepth"], new["EEdepth"])
        
                if base["SHSource"] != new["SHSource"]:
                    changes["sh_source"] = (base["SHSource"], new["SHSource"])
        
                if base["DHWsource"] != new["DHWsource"]:
                    changes["dhw_source"] = (base["DHWsource"], new["DHWsource"])
        
                if base["PVType"] != new["PVType"]:
                    changes["PVtype"] = (base["PVType"], new["PVType"])
        
                if float(base["PVpercentage"]) != float(new["PVpercentage"]):
                    changes["PVpercentage"] = (
                        float(base["PVpercentage"]),
                        float(new["PVpercentage"])
                    )
        
                interventions[idx] = changes
        
            return interventions
        def building_editor_from_dict(geojson_path, config):
        
            gdf = load_buildings(baseline_geojson)
            baseline = gdf.copy()
        
            _validate_config(gdf, config)
        
            gdf = _apply_config(gdf, config)
        
            interventions = build_intervention_dict(baseline, gdf)
        
            return {"gdf": gdf}, interventions
        
        
        def _validate_config(gdf, config):
            pass
        
        
        def _apply_config(gdf, config):
        
            for idx, row in config.items():
                if idx not in gdf["id"].values:
                    raise ValueError(f"{idx}")
        
                env  = row.get("env")
                heat = row.get("heat")
                dhw  = row.get("dhw")
                fuel = row.get("fuel")
                pv_type = row.get("pv_type")
                pv_perc = row.get("pv_percentage")
        
                if env is not None:
                    gdf.loc[gdf["id"] == idx, "EEdepth"] = env
        
                if heat is not None:
                    if heat == "boiler":
                        if not fuel:
                            raise ValueError(f"{idx}")
                        gdf.loc[gdf["id"] ==idx, "SHSource"] = f"boiler_{fuel}"
                    else:
                        gdf.loc[gdf["id"] ==idx, "SHSource"] = heat
        
                if dhw is not None:
                    if dhw == "boiler":
                        if not fuel:
                            raise ValueError(f"{idx}")
                        gdf.loc[gdf["id"] ==idx, "DHWsource"] = f"boiler_{fuel}"
                    else:
                        gdf.loc[gdf["id"] ==idx, "DHWsource"] = dhw
        
                if pv_type is not None:
                    gdf.loc[gdf["id"] ==idx, "PVType"] = pv_type
        
                if pv_perc is not None:
                    gdf.loc[gdf["id"] ==idx, "PVpercentage"] = max(
                        float(gdf.loc[gdf["id"] ==idx, "PVpercentage"]),
                        float(pv_perc)
                    )
            return gdf
        
        q , interv_dict = building_editor_from_dict(baseline_geojson, intervention_dictionary)
        a = demand_aggregation(q["gdf"], city)
    if mode == "app":
        q, interv_dict = app.building_editor(baseline_geojson)


    
    return  interv_dict, a


