import os
import copy
import warnings
from time import time
from pathlib import Path
#%%
import numpy as np
import geopandas as gpd

from eureca_building.config import load_config
from eureca_ubem.city import City

from eureca_pubem import scenario_process as sc
from eureca_pubem import buildings_cost as bc
from eureca_pubem import dhn_costs as dc
from eureca_pubem import grid_costs as gc
from eureca_pubem import market

from visualization import write_map, compute_total_dhn_generation

warnings.filterwarnings("ignore")

_NO_CHANGE = "__NO_CHANGE__"

 
PATHS = {
    "config": r".\Example_District_Config.json",
    "weather": r".\SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
    "schedules": r"TransitionMatrixSchedule.csv",
    "materials": r".\tabula_sverige.xlsx",
    "city_model": r".\soderman_limited_reproject_2.geojson",
    "systems": r".\systems.xlsx",
    "output_folder": r".\bigcheck",

    "distribution": r"C:\Works\EUReCA\EUReCA\eureca_ubem\Input\distribution.geojson",
    "baseline_geojson": r"C:\Works\EUReCA\EUReCA\eureca_ubem\Input\soderman_limited_reproject_2.geojson",
    "roads": r"C:\Works\EUReCA\EUReCA\eureca_ubem\Input\roads.geojson",
    "weather_abs": r"C:\Works\EUReCA\EUReCA\eureca_ubem\Input\SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",

    "dhn_pipes": r"C:\Works\EUReCA\EUReCA\eureca_pubem\dhn\_pipe_diameters.json",
    "cables": r"C:\Works\EUReCA\EUReCA\eureca_pubem\grid\swedish_cable.json",

    "ee_measures": r"C:\Works\EUReCA\EUReCA\eureca_pubem\EE_measures_catalog.json",
    "pv_config": r"C:\Works\EUReCA\EUReCA\eureca_pubem\pv_config.json",
    "hp_config": r"C:\Works\EUReCA\EUReCA\eureca_pubem\hp_config.json",
    "fuels": r"C:\Works\EUReCA\EUReCA\eureca_pubem\fuels.json",
    "spot_price": r"C:\Works\EUReCA\EUReCA\eureca_pubem\spot_price_se3.csv",
}


ASSUMPTIONS = {
    "discount_rate": 0.04,
    "horizon_years": 20,

    "hp_scop": 3.0,
    "boiler_efficiency": 0.9,

    "pv_yield_kwh_per_m2_year": 160.0,
    "pv_self_consumption_ratio": 0.6,

    "dhn_connection_cost": 30_000.0,

    "pv_capex_per_m2": 2500.0,
    "hp_capex_sh": 90_000.0,
    "hp_capex_dhw": 40_000.0,
    "boiler_capex_sh": 35_000.0,
    "boiler_capex_dhw": 20_000.0,
    "generator_capex": 50_000.0,

    "boiler_fuel_cost_per_kwh": 1.2,
    "generator_cost_per_kwh": 2.5,
    "generator_kwh_per_year": 0.0,

    "electricity_spot_price_per_kwh": 1.0,

    "operation_unit": "Wh",

    "primary_energy_factors": {
        "electricity_bought": 1.8,
        "dhn_bought": 0.7,
        "bio": 0.2,
        "biomass": 0.2,
        "pellet": 0.2,
        "wood": 0.2,
        "gas": 1.1,
        "natural_gas": 1.1,
        "oil": 1.1,
        "fuel_default": 1.0,
    },
}


SUPPLIER_COSTS = {
    "grids": {
        0: {
            "electricity_purchase_cost_per_kwh": 1.0,
            "grid_fixed_cost_yearly": 0.0,
        },
    },
    "dhns": {
        "District Heating Supply 8": {
            "heat_supply_cost_per_mwh": 500.0,
            "fixed_cost_yearly": 0.0,
        },
    },
}


OPTIMIZATION_SETTINGS = {
    "max_rounds": 1000,
    "min_rounds": 2,
    "tolerance": 1e-2,

    "local_search": {
        "grid": {
            "active_levers": [
                "pricing.buy.grid local distribution cost monthly fix",
                "pricing.buy.grid local distribution cost monthly per kW peak",
                "pricing.buy.grid local distribution cost monthly per kWh usage",
            ],
            "step_sizes": {
                "pricing.buy.grid local distribution cost monthly fix": 25.0,
                "pricing.buy.grid local distribution cost monthly per kW peak": 2.5,
                "pricing.buy.grid local distribution cost monthly per kWh usage": 0.005,
            },
            "bounds": {
                "pricing.buy.grid local distribution cost monthly fix": (0.0, 3000.0),
                "pricing.buy.grid local distribution cost monthly per kW peak": (0.0, 500.0),
                "pricing.buy.grid local distribution cost monthly per kWh usage": (0.0, 2.0),
            },
        },

        "dhn": {
            "active_levers": [
                "pricing.area fee per m2",
                "pricing.fixed heat price per MWh",
                "pricing.variable heat price per MWh",
                "pricing.admin fee yearly",
                "pricing.subscription fixed yearly per unit",
                "pricing.subscription variable price per MWh",
                "connection_cost",
            ],
            "step_sizes": {
                "pricing.area fee per m2": 1.0,
                "pricing.fixed heat price per MWh": 5.0,
                "pricing.variable heat price per MWh": 10.0,
                "pricing.admin fee yearly": 50.0,
                "pricing.subscription fixed yearly per unit": 100.0,
                "pricing.subscription variable price per MWh": 5.0,
                "connection_cost": 1000.0,
            },
            "bounds": {
                "pricing.area fee per m2": (0.0, 200.0),
                "pricing.fixed heat price per MWh": (0.0, 1000.0),
                "pricing.variable heat price per MWh": (0.0, 5000.0),
                "pricing.admin fee yearly": (0.0, 10000.0),
                "pricing.subscription fixed yearly per unit": (0.0, 50000.0),
                "pricing.subscription variable price per MWh": (0.0, 1000.0),
                "connection_cost": (0.0, 50000.0),
            },
        },
    },
}


GRID_REGULATION = {
    "allowed_revenue_yearly": 8_000_000.0,
    "background_similarity_factor": 0.30,

    "background_baseline_customer_count": 8_000.0,
    "background_baseline_electricity_bought_kwh_year": 32_000_000.0,
    "background_baseline_electricity_sold_kwh_year": 2_000_000.0,
    "background_baseline_peak_kw": 12_000.0,

    "revenue_gap_weight": 1.0,
    "tariff_change_weight": 1_000.0,
}


ZERO_STUDY_AREA_BASELINE_GRID_TOTALS = {
    "customer_count": 0.0,
    "electricity_bought_kwh_year": 0.0,
    "electricity_sold_kwh_year": 0.0,
    "peak_kw": 0.0,
}


def wait_for_user(message, interactive=True):
    if interactive:
        input(f"\n{message}\nPress Enter to continue...")


def load_buildings(input_data):
    if isinstance(input_data, gpd.GeoDataFrame):
        gdf = input_data.copy()
    elif isinstance(input_data, str):
        gdf = gpd.read_file(input_data)
    else:
        raise ValueError("invalid_input")

    required_columns = [
        "id",
        "Name",
        "EEdepth",
        "SHSource",
        "DHWsource",
        "PVType",
        "PVpercentage",
    ]

    missing = [c for c in required_columns if c not in gdf.columns]

    if missing:
        raise ValueError(f"missing_columns: {missing}")

    return gdf


def simulate_city_envelopes():
    output_folder = PATHS["output_folder"]
    os.makedirs(output_folder, exist_ok=True)

    load_config(PATHS["config"])

    cities_envelopes = {}

    for i, envelope in enumerate(["none", "shallow", "medium", "deep"]):
        print(f"\nSimulating envelope scenario: {envelope}")

        gdf = sc.load_geojson_with_envelope_prefix(PATHS["city_model"], i)

        scenario_city = City(
            city_model=gdf,
            epw_weather_file=PATHS["weather"],
            end_uses_types_file=PATHS["schedules"],
            envelope_types_file=PATHS["materials"],
            systems_templates_file=PATHS["systems"],
            shading_calculation=True,
            building_model="2C",
            output_folder=output_folder,
        )

        scenario_city.simulate()
        gdf, value_store = sc.create_gdf_dictionary(scenario_city)
        cities_envelopes[envelope] = value_store

    return gdf, cities_envelopes


def create_initial_baseline(mycity_gdf, mycity):
    baseline = sc.create_baseline(
        input_gdf=mycity_gdf,
        input_city=mycity,
        distributions_geojson=PATHS["distribution"],
        baseline_gdf=PATHS["baseline_geojson"],
        weather_path=PATHS["weather_abs"],
        street_path=PATHS["roads"],
    )

    baseline_gdf = load_buildings(PATHS["baseline_geojson"])
    config0 = bc.extract_config(baseline_gdf)

    interv_dict, building_info = sc.make_dictionary(
        baseline_geojson=PATHS["baseline_geojson"],
        city=mycity,
        baseline_scenario=baseline,
        weatherfile_path=PATHS["weather_abs"],
        intervention_dictionary=config0,
    )

    retrofit, buildings, _ = sc.analyze_intervention(
        baseline_geojson=PATHS["baseline_geojson"],
        city=mycity,
        baseline_scenario=baseline,
        weatherfile_path=PATHS["weather_abs"],
        intervention_dictionary=interv_dict,
    )

    return baseline, baseline_gdf, config0, interv_dict, building_info, retrofit, buildings


def make_intervention_from_config(config, mycity, baseline):
    interv_dict, building_info = sc.make_dictionary(
        baseline_geojson=PATHS["baseline_geojson"],
        city=mycity,
        baseline_scenario=baseline,
        weatherfile_path=PATHS["weather_abs"],
        intervention_dictionary=config,
    )

    retrofit, buildings, _ = sc.analyze_intervention(
        baseline_geojson=PATHS["baseline_geojson"],
        city=mycity,
        baseline_scenario=baseline,
        weatherfile_path=PATHS["weather_abs"],
        intervention_dictionary=interv_dict,
    )

    return interv_dict, building_info, retrofit, buildings


def compute_system_capex(retrofit):
    dhn_pipe_costs = dc.build_cost_table(
        pipe_data=PATHS["dhn_pipes"]
    )

    capital_costs_dhn = dc.compute_dhn_cost(
        dhn_pipe_changes=retrofit.dhn_pipe_changes,
        area_type="urban",
        assumptions={
            "replacement_factor_ground": 0.5,
            "removal_factor": 0.3,
        },
        pipe_json=dhn_pipe_costs,
    )

    capital_costs_grid = gc.compute_grid_cost(
        grid_line_changes=retrofit.grid_line_changes,
        area_type="town",
        assumptions={
            "replacement_factor_ground": 1.0,
            "removal_factor": 0.2,
            "ground_share": 0.55,
            "rest_share": 0.45,
        },
        cable_json=PATHS["cables"],
        cable_key="name",
    )

    return capital_costs_dhn, capital_costs_grid


def run_market_optimization(buildings, retrofit, capital_costs_grid):
    buildings = market.reachables(
        buildings,
        retrofit.District_Heating_Systems,
        retrofit.Electrical_Network,
    )

    buildings = market.build_network_maps(
        buildings,
        retrofit.District_Heating_Systems,
        retrofit.Electrical_Network,
    )

    dhn_levers, grid_levers = market.initialize_market_levers(
        dhns=retrofit.District_Heating_Systems,
        grids=retrofit.Electrical_Network,
    )

    building_option_cache = market.build_building_option_cache(
        buildings=buildings,
        grids=retrofit.Electrical_Network,
        dhns=retrofit.District_Heating_Systems,
        grid_levers=grid_levers,
        dhn_levers=dhn_levers,
        assumptions=ASSUMPTIONS,
    )

    tech_data = market.prepare_global_assumptions(ASSUMPTIONS)

    avg_n_occ = sum(d["meta"]["n_occ"] for d in buildings.values()) / len(buildings)

    market_results = market.optimize_market_levers(
        building_option_cache=building_option_cache,
        initial_grid_levers=grid_levers,
        initial_dhn_levers=dhn_levers,
        tech_data=tech_data,
        supplier_costs=SUPPLIER_COSTS,
        optimization_settings=OPTIMIZATION_SETTINGS,
        grid_regulation=GRID_REGULATION,
        zero_study_area_baseline_grid_totals=ZERO_STUDY_AREA_BASELINE_GRID_TOTALS,
        avg_n_occ=avg_n_occ,
        study_area_grid_capex=capital_costs_grid,
    )

    return buildings, market_results


def prices_from_market_results(market_results):
    grid_levers = market_results["grid_levers"]
    dhn_levers = market_results["dhn_levers"]

    if isinstance(grid_levers, dict):
        first_grid_key = next(iter(grid_levers))
        first_grid_lever = grid_levers[first_grid_key]
    else:
        first_grid_lever = grid_levers[0]

    first_dhn_key = next(iter(dhn_levers))
    first_dhn_lever = dhn_levers[first_dhn_key]

    grid_pricing = {
        "1": first_grid_lever["pricing"]
    }

    dhn_pricing = {
        "1": {
            "buy": first_dhn_lever["pricing"]
        }
    }

    return grid_pricing, dhn_pricing


def build_building_dictionary(
    building_info,
    interv_dict,
    configuration,
    grid_pricing,
    dhn_pricing,
):
    return bc.build_dict_gen(
        building_info,
        interv_dict,
        baseline_gdf_path=PATHS["baseline_geojson"],
        configuration=configuration,
        ee_measure_path=PATHS["ee_measures"],
        pv_type_path=PATHS["pv_config"],
        hp_catalog_path=PATHS["hp_config"],
        grid_pricing_path=grid_pricing,
        dhn_pricing_path=dhn_pricing,
        fuels_path=PATHS["fuels"],
        spot_price_path=PATHS["spot_price"],
    )


def optimize_building_response(
    config_current,
    current_dictionary,
    baseline_dictionary,
    grid_pricing,
    dhn_pricing,
    mycity,
    baseline,
    move_fraction=0.10,
    r=0.04,
    T=25,
):
    return bc.optimize_configuration_per_building_batch_step(
        config_current=config_current,
        current_dictionary=current_dictionary,
        baseline_dictionary=baseline_dictionary,
        baseline_gdf_path=PATHS["baseline_geojson"],
        ee_measure_path=PATHS["ee_measures"],
        pv_type_path=PATHS["pv_config"],
        hp_catalog_path=PATHS["hp_config"],
        grid_pricing_path=grid_pricing,
        dhn_pricing_path=dhn_pricing,
        fuels_path=PATHS["fuels"],
        spot_price_path=PATHS["spot_price"],
        weatherfile_path=PATHS["weather_abs"],
        mycity=mycity,
        baseline_scenario=baseline,
        move_fraction=move_fraction,
        r=r,
        T=T,
    )


def values_are_different(old, new):
    if isinstance(old, np.ndarray) or isinstance(new, np.ndarray):
        try:
            return not np.array_equal(np.asarray(old), np.asarray(new), equal_nan=True)
        except Exception:
            return True

    try:
        return old != new
    except Exception:
        return True


def diff_dicts(old, new, path=""):
    changes = []

    if isinstance(old, dict) and isinstance(new, dict):
        old_keys = set(old.keys())
        new_keys = set(new.keys())

        for key in old_keys - new_keys:
            changes.append((f"{path}{key}", old[key], "__MISSING__"))

        for key in new_keys - old_keys:
            changes.append((f"{path}{key}", "__MISSING__", new[key]))

        for key in old_keys & new_keys:
            changes.extend(diff_dicts(old[key], new[key], path=f"{path}{key}."))

        return changes

    if isinstance(old, list) and isinstance(new, list):
        max_len = max(len(old), len(new))

        for i in range(max_len):
            if i >= len(old):
                changes.append((f"{path}[{i}]", "__MISSING__", new[i]))
            elif i >= len(new):
                changes.append((f"{path}[{i}]", old[i], "__MISSING__"))
            else:
                changes.extend(diff_dicts(old[i], new[i], path=f"{path}[{i}]."))

        return changes

    if values_are_different(old, new):
        clean_path = path[:-1] if path.endswith(".") else path
        changes.append((clean_path, old, new))

    return changes


def count_changed_buildings(changes):
    building_ids = set()

    for key_path, _, _ in changes:
        first = str(key_path).split(".")[0]
        building_ids.add(first)

    return len(building_ids)


def summarize_changed_buildings(changes):
    out = {}

    for key_path, old_value, new_value in changes:
        parts = str(key_path).split(".")
        bid = parts[0]
        field = ".".join(parts[1:]) if len(parts) > 1 else ""

        if bid not in out:
            out[bid] = []

        out[bid].append((field, old_value, new_value))

    return out


def print_market_summary(
    game_step,
    grid_pricing,
    dhn_pricing,
    total_dhn_generation_kwh=None,
    n_changed_buildings=None,
    n_changed_fields=None,
):
    print(f"\nGame step {game_step}")

    if total_dhn_generation_kwh is not None:
        print(f"Total DHN generation: {total_dhn_generation_kwh:,.2f} kWh")

    if n_changed_buildings is not None:
        print(f"Changed buildings: {n_changed_buildings}")

    if n_changed_fields is not None:
        print(f"Changed fields: {n_changed_fields}")

    print("Grid pricing:")
    print(grid_pricing)

    print("DHN pricing:")
    print(dhn_pricing)


def print_config_changes(game_step, changes, limit=80):
    n_buildings = count_changed_buildings(changes)

    if not changes:
        print(f"Stopping at game_step={game_step}: no building changed configuration.")
        return

    print(f"Changes at game_step={game_step}:")
    print(f"  Changed buildings: {n_buildings}")
    print(f"  Changed fields: {len(changes)}")

    for key_path, old_value, new_value in changes[:limit]:
        print(f"  {key_path}: {old_value} -> {new_value}")

    if len(changes) > limit:
        print(f"  Showing first {limit} of {len(changes)} changed fields.")


def make_history_record(
    game_step,
    config,
    optimal_set,
    changes,
    market_results,
    grid_pricing,
    dhn_pricing,
    capital_costs_grid,
    capital_costs_dhn,
    total_dhn_generation_kwh,
):
    return {
        "game_step": game_step,
        "config": copy.deepcopy(config),
        "optimal_set": copy.deepcopy(optimal_set),
        "changes": copy.deepcopy(changes),
        "n_changed_fields": len(changes),
        "n_changed_buildings": count_changed_buildings(changes),
        "changed_buildings_summary": copy.deepcopy(summarize_changed_buildings(changes)),
        "market_results": copy.deepcopy(market_results),
        "grid_pricing": copy.deepcopy(grid_pricing),
        "dhn_pricing": copy.deepcopy(dhn_pricing),
        "capital_costs_grid": copy.deepcopy(capital_costs_grid),
        "capital_costs_dhn": copy.deepcopy(capital_costs_dhn),
        "total_dhn_generation_kwh": total_dhn_generation_kwh,
    }


def run_initial_step(mycity_gdf, mycity, move_fraction=0.10, open_browser=True):
    baseline, baseline_gdf, config0, interv_dict, building_info, retrofit, buildings = create_initial_baseline(
        mycity_gdf=mycity_gdf,
        mycity=mycity,
    )
    print(1)
    capital_costs_dhn, capital_costs_grid = compute_system_capex(retrofit)
    print(2)
    buildings, market_results = run_market_optimization(
        buildings=buildings,
        retrofit=retrofit,
        capital_costs_grid=capital_costs_grid,
    )
    print(3)
    grid_pricing, dhn_pricing = prices_from_market_results(market_results)

    baseline_dictionary = build_building_dictionary(
        building_info=building_info,
        interv_dict=interv_dict,
        configuration=config0,
        grid_pricing=grid_pricing,
        dhn_pricing=dhn_pricing,
    )
    print(4)
    current_optimal_config, current_optimal_set, hist = optimize_building_response(
        config_current=config0,
        current_dictionary=baseline_dictionary,
        baseline_dictionary=baseline_dictionary,
        grid_pricing=grid_pricing,
        dhn_pricing=dhn_pricing,
        mycity=mycity,
        baseline=baseline,
        move_fraction=move_fraction,
        r=0.04,
        T=25,
    )
    print(5)
    initial_changes = diff_dicts(config0, current_optimal_config)

    total_dhn_generation_kwh = compute_total_dhn_generation(
        retrofit.District_Heating_Systems,
        unit="kWh",
    )
    print(6)
    write_map(
        game_step=0,
        buildings_gdf=baseline_gdf,
        retrofit=retrofit,
        config=current_optimal_config,
        optimal_set=current_optimal_set,
        grid_pricing=grid_pricing,
        dhn_pricing=dhn_pricing,
        capital_costs_grid=capital_costs_grid,
        capital_costs_dhn=capital_costs_dhn,
        changes=initial_changes,
        assumptions=ASSUMPTIONS,
        output_folder=PATHS["output_folder"],
        open_browser=open_browser,
    )

    return {
        "baseline": baseline,
        "baseline_gdf": baseline_gdf,
        "config0": config0,
        "baseline_dictionary": baseline_dictionary,
        "current_optimal_config": current_optimal_config,
        "current_optimal_set": current_optimal_set,
        "retrofit": retrofit,
        "buildings": buildings,
        "market_results": market_results,
        "grid_pricing": grid_pricing,
        "dhn_pricing": dhn_pricing,
        "capital_costs_grid": capital_costs_grid,
        "capital_costs_dhn": capital_costs_dhn,
        "changes": initial_changes,
        "total_dhn_generation_kwh": total_dhn_generation_kwh,
    }


def run_one_market_building_step(
    game_step,
    config_current,
    baseline_dictionary,
    mycity,
    baseline,
    buildings_gdf,
    move_fraction=0.10,
    open_browser=False,
):
    interv_dict, building_info, retrofit, buildings = make_intervention_from_config(
        config=config_current,
        mycity=mycity,
        baseline=baseline,
    )

    capital_costs_dhn, capital_costs_grid = compute_system_capex(retrofit)

    buildings, market_results = run_market_optimization(
        buildings=buildings,
        retrofit=retrofit,
        capital_costs_grid=capital_costs_grid,
    )

    grid_pricing, dhn_pricing = prices_from_market_results(market_results)

    current_dictionary = build_building_dictionary(
        building_info=building_info,
        interv_dict=interv_dict,
        configuration=config_current,
        grid_pricing=grid_pricing,
        dhn_pricing=dhn_pricing,
    )

    new_config, optimal_set, hist = optimize_building_response(
        config_current=config_current,
        current_dictionary=current_dictionary,
        baseline_dictionary=baseline_dictionary,
        grid_pricing=grid_pricing,
        dhn_pricing=dhn_pricing,
        mycity=mycity,
        baseline=baseline,
        move_fraction=move_fraction,
        r=0.04,
        T=25,
    )

    changes = diff_dicts(config_current, new_config)

    total_dhn_generation_kwh = compute_total_dhn_generation(
        retrofit.District_Heating_Systems,
        unit="kWh",
    )

    write_map(
        game_step=game_step,
        buildings_gdf=buildings_gdf,
        retrofit=retrofit,
        config=new_config,
        optimal_set=optimal_set,
        grid_pricing=grid_pricing,
        dhn_pricing=dhn_pricing,
        capital_costs_grid=capital_costs_grid,
        capital_costs_dhn=capital_costs_dhn,
        changes=changes,
        assumptions=ASSUMPTIONS,
        output_folder=PATHS["output_folder"],
        open_browser=open_browser,
    )

    return {
        "interv_dict": interv_dict,
        "building_info": building_info,
        "retrofit": retrofit,
        "buildings": buildings,
        "capital_costs_dhn": capital_costs_dhn,
        "capital_costs_grid": capital_costs_grid,
        "market_results": market_results,
        "grid_pricing": grid_pricing,
        "dhn_pricing": dhn_pricing,
        "current_dictionary": current_dictionary,
        "new_config": new_config,
        "optimal_set": optimal_set,
        "changes": changes,
        "total_dhn_generation_kwh": total_dhn_generation_kwh,
    }


def open_final_map():
    final_map_path = Path(PATHS["output_folder"]) / "maps" / "current_map.html"

    print("\nFinal map:")
    print(final_map_path.resolve())

    try:
        import webbrowser
        webbrowser.open(final_map_path.resolve().as_uri())
    except Exception as e:
        print(f"Could not open browser automatically: {e}")
        print("Open this file manually:")
        print(final_map_path.resolve())


def run_visual_game(
    max_game_steps=20,
    move_fraction=0.1,
    interactive=True,
    open_browser_step0=True,
    open_browser_each_step=False,
):
    print("\nStarting envelope simulations.")
    mycity_gdf, mycity = simulate_city_envelopes()
    print("\nEnvelope simulations completed.")

    wait_for_user(
        "Envelope simulations are completed. Next stage: create baseline and run game step 0.",
        interactive=interactive,
    )

    initial = run_initial_step(
        mycity_gdf=mycity_gdf,
        mycity=mycity,
        move_fraction=move_fraction,
        open_browser=open_browser_step0,
    )

    baseline = initial["baseline"]
    baseline_gdf = initial["baseline_gdf"]
    baseline_dictionary = initial["baseline_dictionary"]

    current_optimal_config = copy.deepcopy(initial["current_optimal_config"])
    current_optimal_set = copy.deepcopy(initial["current_optimal_set"])

    history = []

    history.append(
        make_history_record(
            game_step=0,
            config=current_optimal_config,
            optimal_set=current_optimal_set,
            changes=initial["changes"],
            market_results=initial["market_results"],
            grid_pricing=initial["grid_pricing"],
            dhn_pricing=initial["dhn_pricing"],
            capital_costs_grid=initial["capital_costs_grid"],
            capital_costs_dhn=initial["capital_costs_dhn"],
            total_dhn_generation_kwh=initial["total_dhn_generation_kwh"],
        )
    )

    print_market_summary(
        game_step=0,
        grid_pricing=initial["grid_pricing"],
        dhn_pricing=initial["dhn_pricing"],
        total_dhn_generation_kwh=initial["total_dhn_generation_kwh"],
        n_changed_buildings=count_changed_buildings(initial["changes"]),
        n_changed_fields=len(initial["changes"]),
    )

    print_config_changes(
        game_step=0,
        changes=initial["changes"],
    )

    if count_changed_buildings(initial["changes"]) == 0:
        print("\nStopping after step 0: no building changed configuration.")
        open_final_map()

        return {
            "mycity_gdf": mycity_gdf,
            "mycity": mycity,
            "baseline": baseline,
            "baseline_gdf": baseline_gdf,
            "baseline_dictionary": baseline_dictionary,
            "final_config": current_optimal_config,
            "final_optimal_set": current_optimal_set,
            "history": history,
        }

    for game_step in range(1, max_game_steps + 1):
        wait_for_user(
            f"Game step {game_step - 1} completed. Next stage: run game step {game_step}.",
            interactive=interactive,
        )

        previous_config = copy.deepcopy(current_optimal_config)

        result = run_one_market_building_step(
            game_step=game_step,
            config_current=previous_config,
            baseline_dictionary=baseline_dictionary,
            mycity=mycity,
            baseline=baseline,
            buildings_gdf=baseline_gdf,
            open_browser=open_browser_each_step,
            move_fraction=move_fraction,
        )

        new_config = copy.deepcopy(result["new_config"])
        changes = result["changes"]

        print_market_summary(
            game_step=game_step,
            grid_pricing=result["grid_pricing"],
            dhn_pricing=result["dhn_pricing"],
            total_dhn_generation_kwh=result["total_dhn_generation_kwh"],
            n_changed_buildings=count_changed_buildings(changes),
            n_changed_fields=len(changes),
        )

        print_config_changes(
            game_step=game_step,
            changes=changes,
        )

        current_optimal_config = new_config
        current_optimal_set = copy.deepcopy(result["optimal_set"])

        history.append(
            make_history_record(
                game_step=game_step,
                config=current_optimal_config,
                optimal_set=current_optimal_set,
                changes=changes,
                market_results=result["market_results"],
                grid_pricing=result["grid_pricing"],
                dhn_pricing=result["dhn_pricing"],
                capital_costs_grid=result["capital_costs_grid"],
                capital_costs_dhn=result["capital_costs_dhn"],
                total_dhn_generation_kwh=result["total_dhn_generation_kwh"],
            )
        )

        if count_changed_buildings(changes) == 0:
            print(f"\nStopping at game_step={game_step}: no building changed configuration.")
            break

        repeated_changes = diff_dicts(previous_config, current_optimal_config)

        if count_changed_buildings(repeated_changes) == 0:
            print(f"\nStopping at game_step={game_step}: configuration did not change.")
            break

    open_final_map()

    return {
        "mycity_gdf": mycity_gdf,
        "mycity": mycity,
        "baseline": baseline,
        "baseline_gdf": baseline_gdf,
        "baseline_dictionary": baseline_dictionary,
        "final_config": current_optimal_config,
        "final_optimal_set": current_optimal_set,
        "history": history,
    }


# if __name__ == "__main__":
#     a = time()

#     results = run_visual_game(
#         max_game_steps=20,
#         move_fraction=0.10,
#         interactive=False,
#         open_browser_step0=True,
#         open_browser_each_step=True,
#     )

#     b = time()
#     print(f"Finished in {(b - a) / 60:.2f} minutes")



mycity_gdf, mycity = simulate_city_envelopes()
print("\nEnvelope simulations completed.")

wait_for_user(
    "Envelope simulations are completed. Next stage: create baseline and run game step 0.",
    interactive=False,
)

from eureca_pubem import buildings_cost as bc
initial = run_initial_step(
    mycity_gdf=mycity_gdf,
    mycity=mycity,
    move_fraction=0.1,
    open_browser=True,
)

baseline = initial["baseline"]
baseline_gdf = initial["baseline_gdf"]
baseline_dictionary = initial["baseline_dictionary"]

current_optimal_config = copy.deepcopy(initial["current_optimal_config"])
current_optimal_set = copy.deepcopy(initial["current_optimal_set"])

history = []

history.append(
    make_history_record(
        game_step=0,
        config=current_optimal_config,
        optimal_set=current_optimal_set,
        changes=initial["changes"],
        market_results=initial["market_results"],
        grid_pricing=initial["grid_pricing"],
        dhn_pricing=initial["dhn_pricing"],
        capital_costs_grid=initial["capital_costs_grid"],
        capital_costs_dhn=initial["capital_costs_dhn"],
        total_dhn_generation_kwh=initial["total_dhn_generation_kwh"],
    )
)

print_market_summary(
    game_step=0,
    grid_pricing=initial["grid_pricing"],
    dhn_pricing=initial["dhn_pricing"],
    total_dhn_generation_kwh=initial["total_dhn_generation_kwh"],
    n_changed_buildings=count_changed_buildings(initial["changes"]),
    n_changed_fields=len(initial["changes"]),
)

print_config_changes(
    game_step=0,
    changes=initial["changes"],
)

if count_changed_buildings(initial["changes"]) == 0:
    print("\nStopping after step 0: no building changed configuration.")
    open_final_map()

max_game_steps=20
for game_step in range(1, max_game_steps + 1):
    wait_for_user(
        f"Game step {game_step - 1} completed. Next stage: run game step {game_step}.",
        interactive=False,
    )

    previous_config = copy.deepcopy(current_optimal_config)

    result = run_one_market_building_step(
        game_step=game_step,
        config_current=previous_config,
        baseline_dictionary=baseline_dictionary,
        mycity=mycity,
        baseline=baseline,
        buildings_gdf=baseline_gdf,
        open_browser=True,
        move_fraction=0.1,
    )

    new_config = copy.deepcopy(result["new_config"])
    changes = result["changes"]

    print_market_summary(
        game_step=game_step,
        grid_pricing=result["grid_pricing"],
        dhn_pricing=result["dhn_pricing"],
        total_dhn_generation_kwh=result["total_dhn_generation_kwh"],
        n_changed_buildings=count_changed_buildings(changes),
        n_changed_fields=len(changes),
    )

    print_config_changes(
        game_step=game_step,
        changes=changes,
    )

    current_optimal_config = new_config
    current_optimal_set = copy.deepcopy(result["optimal_set"])

    history.append(
        make_history_record(
            game_step=game_step,
            config=current_optimal_config,
            optimal_set=current_optimal_set,
            changes=changes,
            market_results=result["market_results"],
            grid_pricing=result["grid_pricing"],
            dhn_pricing=result["dhn_pricing"],
            capital_costs_grid=result["capital_costs_grid"],
            capital_costs_dhn=result["capital_costs_dhn"],
            total_dhn_generation_kwh=result["total_dhn_generation_kwh"],
        )
    )

    if count_changed_buildings(changes) == 0:
        print(f"\nStopping at game_step={game_step}: no building changed configuration.")
        break

    repeated_changes = diff_dicts(previous_config, current_optimal_config)

    if count_changed_buildings(repeated_changes) == 0:
        print(f"\nStopping at game_step={game_step}: configuration did not change.")
        break

open_final_map()